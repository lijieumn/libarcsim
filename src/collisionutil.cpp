/*
  Copyright ©2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "collisionutil.hpp"

#include "simulation.hpp"
#include "magic.hpp"
#include <omp.h>
using namespace std;

namespace arcsim {

    void collect_leaves(BVHNode *node, vector<BVHNode *> &leaves);

    void assign_features(vector<BVHNode *> &leaves, const Mesh &mesh);

    AccelStruct::AccelStruct(const Mesh &mesh, bool ccd) :
            tree((Mesh &) mesh, ccd), root(tree._root), leaves(mesh.faces.size()) {
        if (root) {
            collect_leaves(root, leaves);
            if (arcsim::magic.use_representative_triangles) {
                assign_features(leaves, mesh);
            }
        }
    }

    void collect_leaves(BVHNode *node, vector<BVHNode *> &leaves) {
        if (node->isLeaf()) {
            int f = node->getFace()->index;
            if (f >= leaves.size())
                leaves.resize(f + 1);
            leaves[f] = node;
        } else {
            collect_leaves(node->getLeftChild(), leaves);
            collect_leaves(node->getRightChild(), leaves);
        }
    }

    void assign_features(vector<BVHNode *> &leaves, const Mesh &mesh) {
        for (int n = 0; n < mesh.nodes.size(); n++) {
            mesh.nodes[n]->assigned = false;
        }
        for (int e = 0; e < mesh.edges.size(); e++) {
            mesh.edges[e]->assigned = false;
        }
        for (int l = 0; l < leaves.size(); l++) {
            Face *face = leaves[l]->_face;
            for (int v = 0; v < 3; v++) {
                Node *node = face->v[v]->node;
                if (node->assigned) {
                    continue;
                }
                leaves[l]->_face->rNode.push_back(node);
                node->assigned = true;
            }
            for (int e = 0; e < 3; e++) {
                Edge *edge = face->adje[e];
                if (edge->assigned) {
                    continue;
                }
                leaves[l]->_face->rEdge.push_back(edge);
                edge->assigned = true;
            }

        }
    }

    void update_accel_struct(AccelStruct &acc) {
        if (acc.root)
            acc.tree.refit();
    }

    void mark_descendants(BVHNode *node, bool active);

    void mark_ancestors(BVHNode *node, bool active);

    void mark_all_inactive(AccelStruct &acc) {
        if (acc.root)
            mark_descendants(acc.root, false);
    }

    void mark_active(AccelStruct &acc, const Face *face) {
        if (acc.root)
            mark_ancestors(acc.leaves[face->index], true);
    }

    void mark_descendants(BVHNode *node, bool active) {
        node->_active = active;
        if (!node->isLeaf()) {
            mark_descendants(node->_left, active);
            mark_descendants(node->_right, active);
        }
    }

    void mark_ancestors(BVHNode *node, bool active) {
        node->_active = active;
        if (!node->isRoot())
            mark_ancestors(node->_parent, active);
    }

    void for_overlapping_faces(BVHNode *node, float thickness,
                               BVHCallback callback) {
        if (node->isLeaf() || !node->_active)
            return;
        for_overlapping_faces(node->getLeftChild(), thickness, callback);
        for_overlapping_faces(node->getRightChild(), thickness, callback);
        for_overlapping_faces(node->getLeftChild(), node->getRightChild(),
                              thickness, callback);
    }

    void find_overlapping_faces_single(BVHNode *node, float thickness, BVHCallback callback) {
        vector<BVHNode *> nodeStack;
        nodeStack.clear();
        nodeStack.push_back(node);
        while (!nodeStack.empty()) {
            BVHNode *node = nodeStack.back();
            nodeStack.pop_back();
            if (node->isLeaf() || !node->_active)
                continue;
            nodeStack.push_back(node->getLeftChild());
            nodeStack.push_back(node->getRightChild());
            find_overlapping_faces_double(node->getLeftChild(), node->getRightChild(), thickness, callback);
        }
    }

    void for_overlapping_faces(BVHNode *node0, BVHNode *node1, float thickness,
                               BVHCallback callback) {
        if (!node0->_active && !node1->_active)
            return;
        if (!overlap(node0->_box, node1->_box, thickness))
            return;
        if (node0->isLeaf() && node1->isLeaf()) {
            Face *face0 = node0->getFace(),
                    *face1 = node1->getFace();
            callback(face0, face1);
        } else if (node0->isLeaf()) {
            for_overlapping_faces(node0, node1->getLeftChild(), thickness, callback);
            for_overlapping_faces(node0, node1->getRightChild(), thickness, callback);
        } else {
            for_overlapping_faces(node0->getLeftChild(), node1, thickness, callback);
            for_overlapping_faces(node0->getRightChild(), node1, thickness, callback);
        }
    }

    void find_overlapping_faces_double(BVHNode *node0, BVHNode *node1, float thickness, BVHCallback callback) {
        vector<pair<BVHNode *, BVHNode *> > nodePairStack;
        nodePairStack.clear();
        nodePairStack.push_back(make_pair(node0, node1));
        while (!nodePairStack.empty()) {
            pair<BVHNode *, BVHNode *> nodePair = nodePairStack.back();
            nodePairStack.pop_back();
            BVHNode *n0 = nodePair.first;
            BVHNode *n1 = nodePair.second;
            if (!n0->_active && !n1->_active)
                continue;
            if (!overlap(n0->_box, n1->_box, thickness))
                continue;
            if (n0->isLeaf() && n1->isLeaf()) {
                Face *face0 = n0->getFace(),
                        *face1 = n1->getFace();
                callback(face0, face1);
            } else if (n0->isLeaf()) {
                nodePairStack.push_back(make_pair(n0, n1->getLeftChild()));
                nodePairStack.push_back(make_pair(n0, n1->getRightChild()));
            } else {
                nodePairStack.push_back(make_pair(n0->getLeftChild(), n1));
                nodePairStack.push_back(make_pair(n0->getRightChild(), n1));
            }
        }

    }

    void for_overlapping_faces(Face *face, BVHNode *node, float thickness,
                               BVHCallback callback) {
        if (!node->_active)
            return;
        BOX box;
        box += face->v[0]->node->x;
        box += face->v[1]->node->x;
        box += face->v[2]->node->x;
        if (!overlap(box, node->_box, thickness)) {
            return;
        }
        if (node->isLeaf()) {
            callback(face, node->getFace());
        } else {
            for_overlapping_faces(face, node->getLeftChild(), thickness, callback);
            for_overlapping_faces(face, node->getRightChild(), thickness, callback);
        }
    }

    vector<BVHNode *> collect_upper_nodes(const vector<AccelStruct *> &accs, int n);

    void for_overlapping_faces(const vector<AccelStruct *> &accs,
                               const vector<AccelStruct *> &obs_accs,
                               double thickness, BVHCallback callback,
                               bool parallel) {
        int nnodes = (int) ceil(sqrt(2 * omp_get_max_threads()));
        vector<BVHNode *> nodes = collect_upper_nodes(accs, nnodes);
        int nthreads = omp_get_max_threads();
        omp_set_num_threads(parallel ? omp_get_max_threads() : 1);
#pragma omp parallel for
        for (int n = 0; n < nodes.size(); n++) {
            if (arcsim::magic.self_contact) {
                for_overlapping_faces(nodes[n], thickness, callback);
            }
            for (int m = 0; m < n; m++) {
                if (arcsim::magic.self_contact) {
                    for_overlapping_faces(nodes[n], nodes[m], thickness, callback);
                }
            }
            for (int o = 0; o < obs_accs.size(); o++)
                if (obs_accs[o]->root)
                    for_overlapping_faces(nodes[n], obs_accs[o]->root, thickness,
                                          callback);
        }
        omp_set_num_threads(nthreads);
    }

    void find_overlapping_faces(const std::vector<AccelStruct *> &accs,
                                const std::vector<AccelStruct *> &obs_accs,
                                double thickness, BVHCallback callback,
                                bool parallel) {
        int nnodes = (int) ceil(sqrt(2 * omp_get_max_threads()));
        vector<BVHNode *> nodes = collect_upper_nodes(accs, nnodes);
        int nthreads = omp_get_max_threads();
        omp_set_num_threads(parallel ? omp_get_max_threads() : 1);
#pragma omp parallel for
        for (int n = 0; n < nodes.size(); n++) {
            if (arcsim::magic.self_contact) {
                find_overlapping_faces_single(nodes[n], thickness, callback);
            }
            for (int m = 0; m < n; m++) {
                if (arcsim::magic.self_contact) {
                    find_overlapping_faces_double(nodes[n], nodes[m], thickness, callback);
                }
            }
            for (int o = 0; o < obs_accs.size(); o++)
                if (obs_accs[o]->root)
                    find_overlapping_faces_double(nodes[n], obs_accs[o]->root, thickness, callback);
        }
        omp_set_num_threads(nthreads);
    }

    void for_overlapping_faces(const std::vector<Face *> &faces,
                               const std::vector<AccelStruct *> &obs_accs,
                               double thickness, BVHCallback callback,
                               bool parallel) {
        // int nthreads = omp_get_max_threads();
        // omp_set_num_threads(parallel ? omp_get_max_threads() : 1);
        for (int f = 0; f < faces.size(); f++) {
// #pragma omp parallel for
            for (int o = 0; o < obs_accs.size(); o++) {
                if (obs_accs[o]->root) {
                    for_overlapping_faces(faces[f], obs_accs[o]->root, thickness, callback);
                }
            }
        }
        // omp_set_num_threads(nthreads);
    }

    void for_faces_overlapping_obstacles(const vector<AccelStruct *> &accs,
                                         const vector<AccelStruct *> &obs_accs,
                                         double thickness, BVHCallback callback,
                                         bool parallel) {
        int nnodes = omp_get_max_threads();
        vector<BVHNode *> nodes = collect_upper_nodes(accs, nnodes);
        int nthreads = omp_get_max_threads();
        omp_set_num_threads(parallel ? omp_get_max_threads() : 1);
#pragma omp parallel for
        for (int n = 0; n < nodes.size(); n++)
            for (int o = 0; o < obs_accs.size(); o++)
                if (obs_accs[o]->root)
                    for_overlapping_faces(nodes[n], obs_accs[o]->root, thickness,
                                          callback);
        omp_set_num_threads(nthreads);
    }

    vector<BVHNode *> collect_upper_nodes(const vector<AccelStruct *> &accs,
                                          int nnodes) {
        vector<BVHNode *> nodes;
        for (int a = 0; a < accs.size(); a++)
            if (accs[a]->root)
                nodes.push_back(accs[a]->root);
        while (nodes.size() < nnodes) {
            vector<BVHNode *> children;
            for (int n = 0; n < nodes.size(); n++)
                if (nodes[n]->isLeaf())
                    children.push_back(nodes[n]);
                else {
                    children.push_back(nodes[n]->_left);
                    children.push_back(nodes[n]->_right);
                }
            if (children.size() == nodes.size())
                break;
            nodes = children;
        }
        return nodes;
    }

    vector<AccelStruct *> create_accel_structs(const vector<Mesh *> &meshes,
                                               bool ccd) {
        vector<AccelStruct *> accs(meshes.size());
        for (int m = 0; m < meshes.size(); m++)
            accs[m] = new AccelStruct(*meshes[m], ccd);
        return accs;
    }

    void destroy_accel_structs(vector<AccelStruct *> &accs) {
        for (int a = 0; a < accs.size(); a++) {
            for (int l = 0; l < accs[a]->leaves.size(); l++) {
                accs[a]->leaves[l]->_face->rEdge.clear();
                accs[a]->leaves[l]->_face->rNode.clear();
            }
            delete accs[a];
        }
    }

    template<typename Prim>
    int find_mesh(const Prim *p, const vector<Mesh *> &meshes) {
        for (int m = 0; m < meshes.size(); m++) {
            const vector<Prim *> &ps = get<Prim>(*meshes[m]);
            if (p->index < ps.size() && p == ps[p->index])
                return m;
        }
        return -1;
    }

    template int find_mesh(const Vert *, const vector<Mesh *> &);

    template int find_mesh(const Node *, const vector<Mesh *> &);

    template int find_mesh(const Edge *, const vector<Mesh *> &);

    template int find_mesh(const Face *, const vector<Mesh *> &);

    const vector<Mesh *> *meshes, *obs_meshes;

}