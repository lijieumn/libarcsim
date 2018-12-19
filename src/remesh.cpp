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

#include "remesh.hpp"
#include "blockvectors.hpp"
#include "cloth.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "localopt.hpp"
#include "magic.hpp"
#include "util.hpp"
#include <assert.h>
#include <cstdlib>
#include <cstdio>
using namespace std;

namespace arcsim {
    RemeshOp RemeshOp::inverse() const {
        RemeshOp iop;
        iop.added_verts = removed_verts;
        iop.removed_verts = added_verts;
        iop.added_nodes = removed_nodes;
        iop.removed_nodes = added_nodes;
        iop.added_edges = removed_edges;
        iop.removed_edges = added_edges;
        iop.added_faces = removed_faces;
        iop.removed_faces = added_faces;
        return iop;
    }

    void RemeshOp::updateImpacts() const {
        vector<Vec2> impactPoints;
        for (int i = 0; i < removed_faces.size(); i++) {
            append(impactPoints, removed_faces[i]->impactPoints);
            removed_faces[i]->impactPoints.clear();
        }
        for (int i = 0; i < removed_edges.size(); i++) {
            append(impactPoints, removed_edges[i]->impactPoints);
            removed_edges[i]->impactPoints.clear();
        }
        for (int i = 0; i < removed_nodes.size(); i++) {
            append(impactPoints, removed_nodes[i]->impactPoints);
            removed_nodes[i]->impactPoints.clear();
        }
        double eps = 1e-6;
        for (int i = 0; i < impactPoints.size(); i++) {
            Vec2 &p = impactPoints[i];
            for (int f = 0; f < added_faces.size(); f++) {
                Face *face = added_faces[f];
                // assert(face->impactPoints.empty());
                Vec3 b = get_barycentric_coords(p, face);
                if (b[0] < -eps || b[1] < -eps || b[2] < -eps) {    // not inside
                    continue;
                }
                int count = 0;  // number of bary that is equal to 0
                int index[2];
                for (int v = 0; v < 3; v++) {
                    if (b[v] < eps) {
                        index[count] = v;
                        count++;
                    }
                }
                if (count == 2) {   // on node
                    Node *n = face->v[3 - index[0] - index[1]]->node;
                    n->impactPoints.push_back(p);
                } else if (count == 1) {    // on edge
                    Node *n1 = face->v[(index[0] + 1) % 3]->node;
                    Node *n2 = face->v[(index[0] + 2) % 3]->node;
                    Edge *edge = getCommonEdge(n1, n2);
                    edge->impactPoints.push_back(p);
                } else if (count == 0) {    // inside face
                    face->impactPoints.push_back(p);
                }
            }
        }
    }

    template<class T>
    static void delete_all(const vector<T> &a) { for (size_t i = 0; i < a.size(); i++) delete a[i]; }

    void RemeshOp::cancel() {
        delete_all(added_verts);
        delete_all(added_nodes);
        delete_all(added_edges);
        delete_all(added_faces);
        added_edges.clear();
        added_faces.clear();
        added_nodes.clear();
        added_verts.clear();
        removed_edges.clear();
        removed_faces.clear();
        removed_nodes.clear();
        removed_verts.clear();
    }

    void RemeshOp::apply(Mesh &mesh) const {
        // cout << "removing " << removed_faces << ", " << removed_edges << ", " << removed_verts << " and adding " << added_verts << ", " << added_edges << ", " << added_faces << endl;
        for (int i = 0; i < removed_faces.size(); i++)
            mesh.remove(removed_faces[i]);
        for (int i = 0; i < removed_edges.size(); i++)
            mesh.remove(removed_edges[i]);
        for (int i = 0; i < removed_nodes.size(); i++)
            mesh.remove(removed_nodes[i]);
        for (int i = 0; i < removed_verts.size(); i++) {
            mesh.remove(removed_verts[i]);
            if (mesh.kdTree) {
                mesh.kdTree->removeVert(removed_verts[i]);
            }
        }
        for (int i = 0; i < added_verts.size(); i++) {
            mesh.add(added_verts[i]);
            if (mesh.kdTree) {
                mesh.kdTree->addVert(added_verts[i]);
            }
        }
        for (int i = 0; i < added_nodes.size(); i++)
            mesh.add(added_nodes[i]);
        for (int i = 0; i < added_edges.size(); i++)
            mesh.add(added_edges[i]);
        for (int i = 0; i < added_faces.size(); i++)
            mesh.add(added_faces[i]);
        // updateImpacts();

        for (size_t f = 0; f < added_faces.size(); f++) {
            compute_ms_data(added_faces[f]);
            for (int i = 0; i < 3; i++)
                compute_ms_data(added_faces[f]->adje[i]);
            for (int i = 0; i < 3; i++)
                compute_ms_data(added_faces[f]->v[i]->node);
        }
    }

    void RemeshOp::done() const {
        for (int i = 0; i < removed_verts.size(); i++) {
            if (removed_verts[i]->sizing) {
                delete removed_verts[i]->sizing;
            }
            delete removed_verts[i];
        }
        for (int i = 0; i < removed_nodes.size(); i++)
            delete removed_nodes[i];
        for (int i = 0; i < removed_edges.size(); i++)
            delete removed_edges[i];
        for (int i = 0; i < removed_faces.size(); i++)
            delete removed_faces[i];
    }

    ostream &operator<<(ostream &out, const RemeshOp &op) {
        out << "removed " << op.removed_verts << ", " << op.removed_nodes << ", "
            << op.removed_edges << ", " << op.removed_faces << ", added "
            << op.added_verts << ", " << op.added_nodes << ", " << op.added_edges
            << ", " << op.added_faces;
        return out;
    }

    template<typename T>
    void compose_removal(T *t, vector<T *> &added, vector<T *> &removed) {
        int i = find(t, added);
        if (i != -1) {
            remove(i, added);
            delete t;
        } else
            removed.push_back(t);
    }

    RemeshOp compose(const RemeshOp &op1, const RemeshOp &op2) {
        RemeshOp op = op1;
        for (int i = 0; i < op2.removed_verts.size(); i++)
            compose_removal(op2.removed_verts[i], op.added_verts, op.removed_verts);
        for (int i = 0; i < op2.removed_nodes.size(); i++)
            compose_removal(op2.removed_nodes[i], op.added_nodes, op.removed_nodes);
        for (int i = 0; i < op2.removed_edges.size(); i++)
            compose_removal(op2.removed_edges[i], op.added_edges, op.removed_edges);
        for (int i = 0; i < op2.removed_faces.size(); i++)
            compose_removal(op2.removed_faces[i], op.added_faces, op.removed_faces);
        for (int i = 0; i < op2.added_verts.size(); i++)
            op.added_verts.push_back(op2.added_verts[i]);
        for (int i = 0; i < op2.added_nodes.size(); i++)
            op.added_nodes.push_back(op2.added_nodes[i]);
        for (int i = 0; i < op2.added_faces.size(); i++)
            op.added_faces.push_back(op2.added_faces[i]);
        for (int i = 0; i < op2.added_edges.size(); i++)
            op.added_edges.push_back(op2.added_edges[i]);
        return op;
    }

// Local pop filter

    void optimize_node(Node *node) {
        vector<Node *> nodes(1, node);
        vector<Face *> faces;
        for (size_t v = 0; v < node->verts.size(); v++)
            append(faces, node->verts[v]->adjf);
        vector<Edge *> edges;
        for (size_t f = 0; f < faces.size(); f++) {
            const Face *face = faces[f];
            for (int i = 0; i < 3; i++)
                include(face->adje[i], edges);
        }
        // vector<Constraint*> no_cons;
        // PlasticityStash stash(faces, edges);
        // local_opt<PS>(nodes, faces, edges, no_cons);
        // stash.apply(faces, edges);
        local_opt<WS>(nodes, faces, edges);
    }

    void local_pop_filter(const vector<Face *> &fs) {
        vector<Node *> nodes;
        vector<Face *> faces;
        vector<Edge *> edges;
        for (size_t f = 0; f < fs.size(); f++)
            for (int i = 0; i < 3; i++)
                include(fs[f]->v[i]->node, nodes);
        for (size_t n = 0; n < nodes.size(); n++) {
            for (size_t v = 0; v < nodes[n]->verts.size(); v++) {
                const Vert *vert = nodes[n]->verts[v];
                for (size_t f = 0; f < vert->adjf.size(); f++)
                    include(vert->adjf[f], faces);
            }
        }
        for (size_t f = 0; f < faces.size(); f++)
            for (int i = 0; i < 3; i++)
                include(faces[f]->adje[i], edges);
        // vector<Constraint*> cons;
        // cons = proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes,
        //                              sim.friction, sim.obs_friction, false, true);
        // for (int h = 0; h < (int)sim.handles.size(); h++)
        //     append(cons, sim.handles[h]->get_constraints(sim.time));

        /*if (sim.frame > 0) {
        for (int i=0; i<nodes.size(); i++)
          Annotation::add(nodes[i]);
        for (int i=0; i<edges.size(); i++)
          Annotation::add(edges[i], Vec3(0,0,1));
        wait_key();}*/
        local_opt<WS>(nodes, faces, edges);
        //if (sim.frame>0) wait_key();
    }

// The actual operations

    int combine_label(int l0, int l1) { return (l0 == l1) ? l0 : 0; }

    int combine_label(int l0, int l1, int l2) {

        if (l0 == l1 && l0 == l2) {
            return l0;
        } else {
            return 0;
        }

    }

    RemeshOp split_edge(Edge *edge) {
        RemeshOp op;
        Node *node0 = edge->n[0],
                *node1 = edge->n[1],
                *node = new Node((node0->y + node1->y) / 2.,
                                 (node0->x + node1->x) / 2.,
                                 (node0->v + node1->v) / 2.,
                                 combine_label(node0->label, node1->label));
        node->acceleration = (node0->acceleration + node1->acceleration) / 2.;
        node->internal_force_density = (node0->internal_force_density + node1->internal_force_density) / 2.;
        node->inMesh = node0->inMesh;
        op.added_nodes.push_back(node);
        op.removed_edges.push_back(edge);
        op.added_edges.push_back(new Edge(node0, node, edge->theta_ideal,
                                          edge->label));
        op.added_edges.push_back(new Edge(node, node1, edge->theta_ideal,
                                          edge->label));
        Vert *vnew[2] = {NULL, NULL};
        for (int s = 0; s < 2; s++) {
            if (!edge->adjf[s])
                continue;
            Vert *v0 = edge_vert(edge, s, s),
                    *v1 = edge_vert(edge, s, 1 - s),
                    *v2 = edge_opp_vert(edge, s);
            if (s == 0 || is_seam_or_boundary(edge)) {
                vnew[s] = new Vert((v0->u + v1->u) / 2.,
                                   combine_label(v0->label, v1->label), v0->component);
                connect(vnew[s], node);
                op.added_verts.push_back(vnew[s]);
            } else
                vnew[s] = vnew[0];

            op.added_edges.push_back(new Edge(v2->node, node));

            Face *f = edge->adjf[s];
            op.removed_faces.push_back(f);
            Face *f0 = new Face(v0, vnew[s], v2, f->label, f->component);
            Face *f1 = new Face(vnew[s], v1, v2, f->label, f->component);
            if (std::min(aspect(f0), aspect(f1)) < 1e-3) {
                op.cancel();
                return op;
            }
            op.added_faces.push_back(f0);
            op.added_faces.push_back(f1);
        }

        if (arcsim::magic.enable_localopt) {
            // embedding_from_plasticity(op.removed_faces);
            op.apply(*edge->n[0]->mesh);
            node->y = (node0->y + node1->y) / 2.;
            optimize_node(node);
            // plasticity_from_embedding(op.added_faces);
            local_pop_filter(op.added_faces);
            op.inverse().apply(*edge->n[0]->mesh);
        }

        return op;
    }


    RemeshOp split_edge(Edge *edge, double b0) {

        double b1 = 1.0 - b0;

        RemeshOp op;
        Node *node0 = edge->n[0],
                *node1 = edge->n[1],
                *node = new Node((b0 * node0->y + b1 * node1->y),
                                 (b0 * node0->x + b1 * node1->x),
                                 (b0 * node0->v + b1 * node1->v),
                                 combine_label(node0->label, node1->label));
        node->acceleration = (b0 * node0->acceleration + b1 * node1->acceleration);
        node->internal_force_density = (b0 * node0->internal_force_density + b1 * node1->internal_force_density);
        node->x0 = (b0 * node0->x0 + b1 * node1->x0); // inserted to test
        node->inMesh = node0->inMesh;
        op.added_nodes.push_back(node);
        op.removed_edges.push_back(edge);
        op.added_edges.push_back(new Edge(node0, node, edge->theta_ideal,
                                          edge->label));
        op.added_edges.push_back(new Edge(node, node1, edge->theta_ideal,
                                          edge->label));
        Vert *vnew[2] = {NULL, NULL};


        for (int s = 0; s < 2; s++) {
            if (!edge->adjf[s])
                continue;
            Vert *v0 = edge_vert(edge, s, s),
                    *v1 = edge_vert(edge, s, 1 - s),
                    *v2 = edge_opp_vert(edge, s);
            if (s == 0 || is_seam_or_boundary(edge)) {
                if (s == 0) {
                    vnew[s] = new Vert((b0 * v0->u + b1 * v1->u),
                                       combine_label(v0->label, v1->label), v0->component);
                } else {
                    vnew[s] = new Vert((b1 * v0->u + b0 * v1->u),
                                       combine_label(v0->label, v1->label), v0->component);
                }

                //if( s == 0 ){
                //	std::cout << "vnew[0] = " << vnew[0]->u << std::endl;
                //} else {
                //	std::cout << "vnew[1] = " << vnew[1]->u << std::endl;
                //}

                connect(vnew[s], node);
                op.added_verts.push_back(vnew[s]);
            } else
                vnew[s] = vnew[0];
            Edge *newEdge = new Edge(v2->node, node);
            newEdge->temp2 = true;
            op.added_edges.push_back(newEdge);
            Face *f = edge->adjf[s];
            op.removed_faces.push_back(f);
            Face *f0 = new Face(v0, vnew[s], v2, f->label, f->component);
            Face *f1 = new Face(vnew[s], v1, v2, f->label, f->component);
            if (std::min(aspect(f0), aspect(f1)) < 1e-3) {
                op.cancel();
                return op;
            }
            op.added_faces.push_back(f0);
            op.added_faces.push_back(f1);
        }

        if (arcsim::magic.enable_localopt) {
            // embedding_from_plasticity(op.removed_faces);
            op.apply(*edge->n[0]->mesh);
            node->y = b0 * node0->y + b1 * node1->y;
            optimize_node(node);
            // plasticity_from_embedding(op.added_faces);
            local_pop_filter(op.added_faces);
            op.inverse().apply(*edge->n[0]->mesh);
        }

        return op;

    }


    RemeshOp collapse_edge(Edge *edge, int i) {
        RemeshOp op;
        Node *node0 = edge->n[i], *node1 = edge->n[1 - i];
        op.removed_nodes.push_back(node0);
        for (int e = 0; e < node0->adje.size(); e++) {
            Edge *edge1 = node0->adje[e];
            op.removed_edges.push_back(edge1);
            Node *node2 = (edge1->n[0] != node0) ? edge1->n[0] : edge1->n[1];
            if (node2 != node1 && !get_edge(node1, node2))
                op.added_edges.push_back(new Edge(node1, node2, edge1->theta_ideal,
                                                  edge1->label));
        }
        for (int s = 0; s < 2; s++) {
            Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
            if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i)))
                continue;
            op.removed_verts.push_back(vert0);
            for (int f = 0; f < vert0->adjf.size(); f++) {
                Face *face = vert0->adjf[f];
                op.removed_faces.push_back(face);
                if (!is_in(vert1, face->v)) {
                    Vert *verts[3] = {face->v[0], face->v[1], face->v[2]};
                    replace(vert0, vert1, verts);
                    Face *newFace = new Face(verts[0], verts[1], verts[2], face->label, face->component);
                    op.added_faces.push_back(newFace);
                    // degenerate test
                    double asp_old = aspect(face), asp_new = aspect(newFace);
                    double asp_min = node0->mesh->parent->remeshing.aspect_min;
                    if (asp_new < asp_min / 4 && asp_old >= asp_min / 4) {
                        op.cancel();
                        return RemeshOp();
                    }
                }
            }
        }

        if (arcsim::magic.enable_localopt) {
            // embedding_from_plasticity(op.removed_faces);
            op.apply(*edge->n[0]->mesh);
            // plasticity_from_embedding(op.added_faces);
            local_pop_filter(op.added_faces);
            op.inverse().apply(*edge->n[0]->mesh);
        }

        return op;
    }


    RemeshOp flip_edge(Edge *edge) {
        RemeshOp op;
        Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 1, 1),
                *vert2 = edge_opp_vert(edge, 0), *vert3 = edge_opp_vert(edge, 1);
        Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
        op.removed_edges.push_back(edge);
        op.added_edges.push_back(new Edge(vert2->node, vert3->node,
                                          -edge->theta_ideal, edge->label));
        op.removed_faces.push_back(face0);
        op.removed_faces.push_back(face1);
        op.added_faces.push_back(new Face(vert0, vert3, vert2, face0->label, face0->component));
        op.added_faces.push_back(new Face(vert1, vert2, vert3, face1->label, face1->component));

        if (arcsim::magic.enable_localopt) {
            // embedding_from_plasticity(op.removed_faces);
            op.apply(*edge->n[0]->mesh);
            // plasticity_from_embedding(op.added_faces);
            local_pop_filter(op.added_faces);
            op.inverse().apply(*edge->n[0]->mesh);
        }

        return op;
    }


    RemeshOp split_face(Face *face, double b0, double b1, double b2) {
        RemeshOp op;

        //std::cout << "split face start\n";

        // Three verts of the input face
        Vert *v0 = face->v[0];
        Vert *v1 = face->v[1];
        Vert *v2 = face->v[2];

        // New vert to be placed inside the input face
        Vec2 baryU = b0 * v0->u + b1 * v1->u + b2 * v2->u;
        int vertLabel = combine_label(v0->label, v1->label, v2->label);
        Vert *vert = new Vert(baryU, vertLabel, face->component);
        op.added_verts.push_back(vert);

        // Creating the three new faces and removing the original input face
        Face *f01n = new Face(v0, v1, vert, face->label, face->component);
        Face *f12n = new Face(v1, v2, vert, face->label, face->component);
        Face *f21n = new Face(v2, v0, vert, face->label, face->component);
        op.removed_faces.push_back(face);
        op.added_faces.push_back(f01n);
        op.added_faces.push_back(f12n);
        op.added_faces.push_back(f21n);

        // Corresponding node for each vert of input face
        Node *n0 = v0->node;
        Node *n1 = v1->node;
        Node *n2 = v2->node;

        // New node corresponding to the new vert
        Vec3 baryY = b0 * n0->y + b1 * n1->y + b2 * n2->y;
        Vec3 baryX = b0 * n0->x + b1 * n1->x + b2 * n2->x;
        Vec3 baryX0 = b0 * n0->x0 + b1 * n1->x0 + b2 * n2->x0;
        Vec3 baryV = b0 * n0->v + b1 * n1->v + b2 * n2->v;
        Vec3 baryAcc = b0 * n0->acceleration + b1 * n1->acceleration + b2 * n2->acceleration;
        int nodeLabel = combine_label(n0->label, n1->label, n2->label);
        Node *node = new Node(baryY, baryX, baryV, nodeLabel);
        node->temp = true;
        node->acceleration = baryAcc;
        node->x0 = baryX0;  // not sure if this should even be included.
        node->inMesh = n0->inMesh;
        op.added_nodes.push_back(node);

        // Connecting new vert to new node
        connect(vert, node);

        // Creating the three new edges linking the face nodes to the new central one
        Edge *e0n = new Edge(n0, node, 0, 0);
        Edge *e1n = new Edge(n1, node, 0, 0);
        Edge *e2n = new Edge(n2, node, 0, 0);
        op.added_edges.push_back(e0n);
        op.added_edges.push_back(e1n);
        op.added_edges.push_back(e2n);

        //std::cout << "split face end\n";

        return op;

    }

    RemeshOp collapse_face_vert(Vert *vert) {


        RemeshOp op;

        // Marking the vert for deletion
        op.removed_verts.push_back(vert);

        // The three faces attached to the vert
        Face *face0 = vert->adjf[0];
        Face *face1 = vert->adjf[1];
        Face *face2 = vert->adjf[2];

        // Marking the faces for deletion
        op.removed_faces.push_back(face0);
        op.removed_faces.push_back(face1);
        op.removed_faces.push_back(face2);

        // Node attached to the vert
        Node *node = vert->node;

        // The three edges attached to the node (by definition we always assume exactly 3)
        Edge *edge0 = node->adje[0];
        Edge *edge1 = node->adje[1];
        Edge *edge2 = node->adje[2];

        // The three unique nodes attached to these edges
        Node *node0 = (edge0->n[1] == node) ? edge0->n[0] : edge0->n[1];
        Node *node1 = (edge1->n[1] == node) ? edge1->n[0] : edge1->n[1];
        Node *node2 = (edge2->n[1] == node) ? edge2->n[0] : edge2->n[1];

        // The verts corresponding to the nodes
        Vert *vert0 = node0->verts[0];
        Vert *vert1 = node1->verts[0];
        Vert *vert2 = node2->verts[0];

        // Marking the node and edges for deletion
        op.removed_nodes.push_back(node);
        op.removed_edges.push_back(edge0);
        op.removed_edges.push_back(edge1);
        op.removed_edges.push_back(edge2);

        // New face which will be shared by node
        Face *face = new Face(vert0, vert1, vert2, face0->label, face0->component);
        //face->label = 0;
        op.added_faces.push_back(face);

        return op;

    }

/*RemeshOp collapse_edge_node(Node* node){
	
	RemeshOp op;
	
	// The vert corresponding to the node
	Vert* vert = node->verts[0];
	
	// Marking the vert for deletion
	op.removed_verts.push_back(vert);
	
	// Marking the node for deletion
	op.removed_nodes.push_back(node);
	
	// The edges attached to the node
	std::vector<Edge*> adje = node -> adje;
	
	// Marking the edges for deletion
	for(int i = 0; i < adje.size(); i++){
		op.removed_edges.push_back( adje[i] );
	}
	
	// The faces adjacent to the vert
	std::vector<Face*> adjf = vert -> adjf;
	
	// Marking the faces for deletion
	for(int i = 0; i < adjf.size(); i++){
		op.removed_faces.push_back( adjf[i] );
	}
	
	// Looping through the edges and connecting non-temp edges
	for(int i = 0; i < adje.size(); i++){
		
		if( adje[i]->temp2 ){
			continue;
		}
		
		for(int j = i+1; j < adje.size(); j++){
		
			if( adje[j]->temp2 ){
				continue;
			}
			
			Node* node0 = (adje[i]->n[1] == node) ? adje[i]->n[0] : adje[i]->n[1];
			Node* node1 = (adje[j]->n[1] == node) ? adje[j]->n[0] : adje[j]->n[1];
			Edge* newEdge = new Edge(node0,node1,0,0);
			op.added_edges.push_back( newEdge );
			
			for(int k = 0; k < adje.size(); k++){
				
				if( !adje[k]->temp2 ){
					continue;
				}
				
				if( adje[k]->n[1] == node ){
					Node* node2 = adje[k]->n[0];
					Face* newFace = new Face(node0->verts[0], node1->verts[0], node2->verts[0]);
					
					newFace->Dm = Mat2x2(newFace->v[1]->u - newFace->v[0]->u,
										 newFace->v[2]->u - newFace->v[0]->u); 
					newFace->invDm = newFace->Dm.inv();
					newFace->a = det(newFace->Dm)/2;
					
					if( newFace->a < 0 ){
						newFace = new Face(node0->verts[0], node2->verts[0], node1->verts[0]);
					} 
					op.added_faces.push_back( newFace );
				} 
					
				
				//Node* node2 = (adje[k]->n[1] == node) ? adje[k]->n[0] : adje[k]->n[1];
				
				//Face* newFace = new Face(node0->verts[0], node1->verts[0], node2->verts[0]);
				//op.added_faces.push_back( newFace );
				
			}
			
		}
	}
	
	//exit(0);

	

	
	return op;
	
		

*/



    // Each time an edge is deleted, we will want to merge the two faces it bisects

    // Edges which each contain the node we want to delete, but themselves do not share an edge,
    // will be made such that they do share an edge, with the deletion of the node in question
    // All other edges we simply delete without question, and we merge the faces bisected by them



    /*
    for(int i = 0; i < cloth.mesh.nodes.size(); i++){

        if( cloth.mesh.nodes[i] -> temp2 ){

            Node* node = cloth.mesh.nodes[i];

            for(int j = 0; j < node->adje.size(); j++){

                Edge* edge = node->adje[j];

                for(int k = 0; k < edge->adjf.size)(); k++){

                    Face*

                }

                op.removed_edges.push_back( edge );

            }

            op.removed_nodes.push_back( cloth.mesh.nodes[i] );

        }
    }

    op.apply(cloth.mesh);  */

//}

}