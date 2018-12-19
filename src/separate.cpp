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

#include "separate.hpp"

#include "collisionutil.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "magic.hpp"
#include "optimization.hpp"
#include "simulation.hpp"
#include "util.hpp"
#include "log.hpp"
#include <omp.h>
#include <set>
using namespace std;

namespace arcsim {
    static const int max_iter = 100;
    static const double &thickness = arcsim::magic.projection_thickness;

    static const vector<Mesh *> *old_meshes;
    static vector<Vec3> xold;
    static vector<Vec3> xold_obs;

    typedef Vec3 Bary; // barycentric coordinates
    typedef std::pair<Vec3, Vec3> Line3;

    struct EdgeIxn {
        Edge *edge;
        Face *face;
        mutable Vec3 G;

        EdgeIxn() :
                edge(NULL), face(NULL), G(0) {}

        EdgeIxn(Edge *edge, Face *face) :
                edge(edge), face(face), G(0) {}

        bool operator<(const EdgeIxn &rhs) const {
            // Don't use any variable with `mutable` qualifer here.
            return (this->edge < rhs.edge
                    || (this->edge == rhs.edge && this->face < rhs.face));
        }
    };

    struct Ixn {// intersection
        Face *f0, *f1;
        Bary b0, b1;
        Vec3 n;

        Ixn() {}

        Ixn(const Face *f0, const Bary &b0, const Face *f1, const Bary &b1,
            const Vec3 &n) : f0((Face *) f0), f1((Face *) f1), b0(b0), b1(b1), n(n) {}
    };

    struct cmpOneAxis {
        size_t axis;

        cmpOneAxis(size_t axis)
                : axis(axis) {}

        bool operator()(const Vec3 &lhs, const Vec3 &rhs) const {
            return lhs[axis] < rhs[axis];
        }
    };

    ostream &operator<<(ostream &out, const Ixn &ixn) {
        out << ixn.f0 << "@" << ixn.b0 << " " << ixn.f1 << "@" << ixn.b1 << " " << ixn.n;
        return out;
    }

    void update_active(const vector<AccelStruct *> &accs, const vector<Ixn> &ixns);

    int major_axis(const Vec3 &v);

    vector<Ixn> find_intersections(const vector<AccelStruct *> &accs,
                                   const vector<AccelStruct *> &obs_accs);

    vector<Ixn> find_intersections(const vector<Face *> &faces, const vector<AccelStruct *> &obs_accs);

    bool edge_face_intersection(const Edge *edge, const Face *face, Vec3 &pt);

    bool edge_face_intersection(const Vec3 &e0, const Vec3 &e1,
                                const Vec3 &f0, const Vec3 &f1, const Vec3 &f2, const Vec3 &fn, Vec3 &pt);

    bool get_line_of_intersection(const Face *face0, const Face *face1, Line3 &line);

    void solve_ixns(const vector<Ixn> &ixns);

    void init_meshes(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &old_meshes,
                     const std::vector<Mesh *> &obs_meshes) {
        arcsim::meshes = &meshes;
        arcsim::old_meshes = &old_meshes;
        arcsim::obs_meshes = &obs_meshes;
    }

    bool
    has_intersection(const std::vector<Face *> &added_faces, Mesh *mesh, const std::vector<AccelStruct *> &obs_accs) {
        vector<Ixn> new_ixns = find_intersections(added_faces, obs_accs);
        //vector<Ixn> new_ixns = find_intersections(accs, obs_accs);
        if (!new_ixns.empty()) return true;

        vector<Mesh *> meshes;
        meshes.push_back(mesh);
        vector<AccelStruct *> accs = create_accel_structs(meshes, false);
        new_ixns = find_intersections(added_faces, accs);
        destroy_accel_structs(accs);

        return !new_ixns.empty();
    }

    bool edge_face_intersection(const Vec3 &e0, const Vec3 &e1,
                                const Vec3 &x0, const Vec3 &x1, const Vec3 &x2, const Vec3 &n, Vec3 &pt) {
        const double numer = dot(x0 - e0, n);
        const double denom = dot(e1 - e0, n);
        if (std::abs(numer) < 1e-12 || std::abs(denom) < 1e-12)
            return false;

        // Find a parameter t in [0, 1].
        const double t = numer / denom;
        if (t < 0.0 || t > 1.0)
            return false;

        // Find a point of intersection.
        pt = e0 + t * (e1 - e0);

        // Check CCW (insdie) or CW (outside).
        const Vec3 &c0 = cross(x1 - x0, pt - x0);
        const Vec3 &c1 = cross(x2 - x1, pt - x1);
        const Vec3 &c2 = cross(x0 - x2, pt - x2);

        bool is_inside = (dot(n, c0) >= -1e-12
                          && dot(n, c1) >= -1e-12
                          && dot(n, c2) >= -1e-12);
        return is_inside;
    }

    bool edge_face_intersection(const Edge *edge, const Face *face, Vec3 &pt) {
        // edge
        const Vec3 &e0 = edge->n[0]->x;
        const Vec3 &e1 = edge->n[1]->x;

        // face
        const Vec3 &f0 = face->v[0]->node->x;
        const Vec3 &f1 = face->v[1]->node->x;
        const Vec3 &f2 = face->v[2]->node->x;
        const Vec3 &fn = nor<WS>(face);

        return edge_face_intersection(e0, e1, f0, f1, f2, fn, pt);
    }

    void ICM_separate(vector<Mesh *> &meshes, const vector<Mesh *> &old_meshes,
                      const vector<Mesh *> &obs_meshes, bool log) {

        arcsim::meshes = &meshes;
        arcsim::old_meshes = &old_meshes;
        arcsim::obs_meshes = &obs_meshes;
        arcsim::xold = node_positions(meshes);
        arcsim::xold_obs = node_positions(obs_meshes);
        for (int m = 0; m < meshes.size(); m++) {
            compute_ws_data(*meshes[m]);
        }
        vector<AccelStruct *> accs = create_accel_structs(meshes, false),
                obs_accs = create_accel_structs(obs_meshes, false);

        // double D = 1e-3;  // step size
        double D = arcsim::magic.projection_thickness * .1;

#ifndef SILENCE_ARGUS
        std::cout << "Separation with step_size: " << D << std::endl;
#endif

        std::set<EdgeIxn> edge_ixns;
        std::map<Node *, Vec3> node_to_G;

        double x_old = 1.0;  // Sum of contour length
        double x_new = 0.0;  // Sum of contour length
        static int frame = 0;
        Log *logger = Log::getInstance();
        VisualDebugger *vd = VisualDebugger::getInstance();
        int maxIter = 200;
        for (size_t iter = 0; iter < maxIter; ++iter) {
            const bool print_logs = (iter % 16 == 0);

            x_old = x_new;
            x_new = 0.0;

            const vector<Ixn> &new_ixns = find_intersections(accs, obs_accs);
            if (new_ixns.empty()) {
                if (log) {
                    if (iter == 0) {
                        logger->setIntersectionNumber(new_ixns.size());
                    }
                    logger->setIntersectionNumberAfter(new_ixns.size());
                }
                break;
            }

            // Find the edge-face intersections.
            // (All I need are just two Face* variables that might be intersecting.)
            edge_ixns.clear();
            for (size_t ixn_idx = 0; ixn_idx < new_ixns.size(); ++ixn_idx) {
                const Ixn &ixn = new_ixns[ixn_idx];

                for (size_t edge_idx = 0; edge_idx < 3; ++edge_idx) {
                    Vec3 pt_ixn;

                    // Check each face's three edges and find all intersections.
                    Edge *edge_ixn = ixn.f0->adje[edge_idx];
                    bool is_ixn = edge_face_intersection(edge_ixn, ixn.f1, pt_ixn);
                    if (is_ixn) {
                        edge_ixns.insert(EdgeIxn(edge_ixn, ixn.f1));
                        if (iter == 0) {
                            Vec3 e0 = edge_ixn->n[0]->x;
                            Vec3 e1 = edge_ixn->n[1]->x;

                            // vd->addVisualPoint3(pt_ixn, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(e0, e1, Vec3(1, 0, 0), 'A');
                            Vec3 f0 = ixn.f1->v[0]->node->x;
                            Vec3 f1 = ixn.f1->v[1]->node->x;
                            Vec3 f2 = ixn.f1->v[2]->node->x;
                            vd->addVisualLine3(f0, f1, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(f1, f2, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(f2, f0, Vec3(1, 0, 0), 'A');
                        }
                    }

                    // Check each face's three edges and find all intersections.
                    edge_ixn = ixn.f1->adje[edge_idx];
                    is_ixn = edge_face_intersection(edge_ixn, ixn.f0, pt_ixn);
                    if (is_ixn) {
                        edge_ixns.insert(EdgeIxn(edge_ixn, ixn.f0));
                        if (iter == 0) {
                            Vec3 e0 = edge_ixn->n[0]->x;
                            Vec3 e1 = edge_ixn->n[1]->x;

                            // vd->addVisualPoint3(pt_ixn, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(e0, e1, Vec3(1, 0, 0), 'A');
                            Vec3 f0 = ixn.f0->v[0]->node->x;
                            Vec3 f1 = ixn.f0->v[1]->node->x;
                            Vec3 f2 = ixn.f0->v[2]->node->x;
                            vd->addVisualLine3(f0, f1, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(f1, f2, Vec3(1, 0, 0), 'A');
                            vd->addVisualLine3(f2, f0, Vec3(1, 0, 0), 'A');
                        }
                    }
                }
            }
            if (iter == 0 && log) {
                logger->setIntersectionNumber(edge_ixns.size());
            }
            // out << endl;
            // break;
            if (edge_ixns.empty()) {
                if (log) logger->setIntersectionNumberAfter(edge_ixns.size());
                break;
            }
            if (iter == maxIter - 1 && log) {
                logger->setIntersectionNumber(edge_ixns.size());
            }
            // A sort of gradient descent method.
            // Compute gradient vector G and displace edge E and face A.
            for (std::set<EdgeIxn>::iterator it = edge_ixns.begin(); it != edge_ixns.end(); ++it) {
                Edge *E = it->edge;  // Edge E
                Face *A = it->face;  // Face A
                const Vec3 &N = nor<WS>(A);

                Vec3 G;  // Gradient vector.
                for (size_t face_idx = 0; face_idx < 2; ++face_idx) {
                    Face *Bi = E->adjf[face_idx];  // The edge's adjacent faces
                    if (Bi == NULL)
                        continue;
                    const Vec3 &Mi = nor<WS>(Bi);  // The surface normal of face B_i

                    Line3 l;
                    bool is_intersect = get_line_of_intersection(Bi, A, l);
                    // Q. What will happen if the function returns a point instead of a line?
                    if (!is_intersect) {
#ifndef SILENCE_ARGUS
                        std::cout << "One of the faces is not intersecting. Ignore it..."
                                  << std::endl;
#endif
                        continue;
                    }
                    x_new += 2.0 * norm(l.second - l.first);

                    // Get the direction vector of edge E.
                    // Vec3 E_vec = E->direction();
                    // to-do determine the direction.
                    Vec3 E_vec = normalize(E->n[1]->x - E->n[0]->x);
                    if (face_idx == 1)
                        E_vec = -E_vec; // Make sure that the edge is CCW w.r.t. Bi.
                    // Add an assertion here.

                    // Get the line of intersection between faces A and B_i.
                    Vec3 Ri_vec = cross(N, Mi);
                    // Check wheter faces A and Bi are parallel or not.
                    if (norm2(Ri_vec) < 1e-12) {
#ifndef SILENCE_ARGUS
                        std::cout << "The faces are parallel. Ignore it..." << std::endl;
#endif
                        continue;
                    }
                    Ri_vec = normalize(Ri_vec);
                    // Check the direction of Ri is correct or not.
                    if (dot(cross(E_vec, Ri_vec), Mi) < 0.0)
                        Ri_vec = -Ri_vec;

                    // Use edge_opp_vert(E, face_idx) to double-check the orientation.
                    Vert *opp_vert = edge_opp_vert(E, face_idx);
                    assert(opp_vert != NULL);
                    Vec3 inside_dir = normalize(opp_vert->node->x - E->n[0]->x);
                    // if (dot(inside_dir, Ri_vec) < 0.0)
                    //    Ri_vec = - Ri_vec;
                    // assert(dot(inside_dir, Ri_vec) >= 0.0);

                    const double E_dot_N = dot(E_vec, N);
                    if (std::abs(E_dot_N) < 1e-12) {
#ifndef SILENCE_ARGUS
                        std::cout << "The edge and the face are parallel. Ignore it..."
                                  << std::endl;
#endif
                        continue;
                    }

                    // Compute gradient vector G.
                    G += Ri_vec - 2.0 * dot(E_vec, Ri_vec) / E_dot_N * N;
                }
                // Store the gradient vector multiplied by the step size.
                it->G = G;

                if (iter == 0) {
                    Vec3 n0 = A->v[0]->node->x;
                    Vec3 n1 = A->v[1]->node->x;
                    Vec3 n2 = A->v[2]->node->x;
                    vd->addVisualLine3(n0, n1, Vec3(0, 0, 1), 'P');
                    vd->addVisualLine3(n1, n2, Vec3(0, 0, 1), 'P');
                    vd->addVisualLine3(n2, n0, Vec3(0, 0, 1), 'P');
                    n0 = E->n[0]->x;
                    n1 = E->n[1]->x;
                    vd->addVisualLine3(n0, n1, Vec3(0, 0, 1), 'P');
                    Vec3 mid = (n0 + n1) * .5;
                    vd->addVisualLine3(mid, mid + it->G * 1e-3, Vec3(0, 0, 1), 'P');
                }
            }

            //     if (print_logs)
            //         std::cout << "n_ixns: " << new_ixns.size()
            //             << " n_edge_ixns: " << edge_ixns.size()
            //             << " x_new: " << x_new << std::endl;

            //     // Aggregate the gradient vectors.
            node_to_G.clear();
            for (std::set<EdgeIxn>::const_iterator it = edge_ixns.begin();
                 it != edge_ixns.end(); ++it) {
                Edge *E = it->edge;  // Edge E
                Face *A = it->face;  // Face A
                const Vec3 &G = it->G;

                for (size_t axis = 0; axis < 2; ++axis)
                    node_to_G[E->n[axis]] += G;
                for (size_t axis = 0; axis < 3; ++axis)
                    node_to_G[A->v[axis]->node] -= G;
            }

            //     // Compute the norm of the aggregated gradient vector.
            double norm2_Gs = 0.0;
            for (std::map<Node *, Vec3>::const_iterator it = node_to_G.begin();
                 it != node_to_G.end(); ++it) {
                norm2_Gs += norm2(it->second);
            }
            const double norm_Gs = sqrt(norm2_Gs);

            // Normalize the displacement D using the norm so that
            // |D*(aggregated_gradient_vector)/(norm_Gs)| == D.
            const double D_ = D / norm_Gs;

            // Apply the displacements.
            for (std::map<Node *, Vec3>::const_iterator it = node_to_G.begin();
                 it != node_to_G.end(); ++it) {
                if (is_free(it->first)) {
                    it->first->x += it->second * D_;
                    if (iter > 0) {
                        continue;
                    }
                    vd->addVisualPoint3(it->first->x, Vec3(1, 0, 0), 'I');
                    vd->addVisualLine3(it->first->x, it->first->x + it->second * 1e-3, Vec3(1, 0, 0), 'I');
                }
            }

            // Update the world-space information.
            for (int m = 0; m < meshes.size(); m++) {
                compute_ws_data(*meshes[m]);
                update_accel_struct(*accs[m]);
            }
            for (int o = 0; o < obs_meshes.size(); o++) {
                compute_ws_data(*obs_meshes[o]);
                update_accel_struct(*obs_accs[o]);
            }
        }

        // // Wrap-up.
        for (int m = 0; m < meshes.size(); m++) {
            // Save collision-free configurations.
            update_x0(*meshes[m]);
        }
        for (int o = 0; o < obs_meshes.size(); o++) {
            // Save collision-free configurations.
            update_x0(*obs_meshes[o]);
        }
        destroy_accel_structs(accs);
        destroy_accel_structs(obs_accs);

#ifndef SILENCE_ARGUS
        std::cout << "Separation done." << std::endl;
#endif
    }

    void separate(vector<Mesh *> &meshes, const vector<Mesh *> &old_meshes,
                  const vector<Mesh *> &obs_meshes) {
        arcsim::meshes = &meshes;
        arcsim::old_meshes = &old_meshes;
        arcsim::obs_meshes = &obs_meshes;
        arcsim::xold = node_positions(meshes);
        vector<AccelStruct *> accs = create_accel_structs(meshes, false),
                obs_accs = create_accel_structs(obs_meshes, false);
        vector<Ixn> ixns;
        int iter;
        for (iter = 0; iter < max_iter; iter++) {
            if (!ixns.empty())
                update_active(accs, ixns);
            vector<Ixn> new_ixns = find_intersections(accs, obs_accs);
            if (new_ixns.empty())
                break;
            append(ixns, new_ixns);
            solve_ixns(ixns);
            for (int m = 0; m < meshes.size(); m++) {
                compute_ws_data(*meshes[m]);
                update_accel_struct(*accs[m]);
            }
        }
        if (iter == max_iter) {
            cerr << "Post-remeshing separation failed to converge!" << endl;
            debug_save_meshes(meshes, "meshes");
            debug_save_meshes(old_meshes, "oldmeshes");
            debug_save_meshes(obs_meshes, "obsmeshes");
            exit(1);
        }
        for (int m = 0; m < meshes.size(); m++) {
            compute_ws_data(*meshes[m]);
            update_x0(*meshes[m]);
        }
        destroy_accel_structs(accs);
        destroy_accel_structs(obs_accs);
    }

    Vec3 pos(const Face *face, const Bary &b) {
        return b[0] * face->v[0]->node->x
               + b[1] * face->v[1]->node->x
               + b[2] * face->v[2]->node->x;
    }

    Vec3 old_pos(const Face *face, const Bary &b) {

        if (!is_free(face))
            return pos(face, b);

        Vec2 u = b[0] * face->v[0]->u + b[1] * face->v[1]->u + b[2] * face->v[2]->u;

        int m;
        for (m = 0; m < arcsim::meshes->size(); m++) {
            if ((*arcsim::meshes)[m]->faces[face->index] == face) {
                break;
            }
        }

        Mesh *oldMesh = (*arcsim::old_meshes)[m];
        if (oldMesh == 0 || oldMesh == NULL) {
            std::cout << "Error: old mesh is bad \n";
            exit(0);
        }
        Face *old_face = get_enclosing_face(*(*arcsim::old_meshes)[m], u);
        Bary old_b = get_barycentric_coords(u, old_face);

        return pos(old_face, old_b);
    }

    void update_active(const vector<AccelStruct *> &accs, const vector<Ixn> &ixns) {
        // for (int a = 0; a < accs.size(); a++)
        //     mark_all_inactive(*accs[a]);
        // for (int i = 0; i < ixns.size(); i++)
        //     for (int f = 0; f < 2; f++) {
        //         const Face *face = f==0 ? ixns[i].f0 : ixns[i].f1;
        //         int m = find(face, *::meshes);
        //         if (m == -1)
        //             continue;
        //         mark_active(*accs[m], face);
        //     }
    }

    static int nthreads = 0;
    static vector<Ixn> *ixns = NULL;

    void find_face_intersection(const Face *face0, const Face *face1);

    int intersection_number(const vector<AccelStruct *> &accs,
                            const vector<AccelStruct *> &obs_accs) {
        return find_intersections(accs, obs_accs).size();
    }

    vector<Ixn> find_intersections(const vector<AccelStruct *> &accs,
                                   const vector<AccelStruct *> &obs_accs) {

        if (!arcsim::ixns) {
            arcsim::nthreads = omp_get_max_threads();
            arcsim::ixns = new vector<Ixn>[arcsim::nthreads];
        }
        for (int t = 0; t < arcsim::nthreads; t++)
            arcsim::ixns[t].clear();
        if (arcsim::magic.use_stack_overloop_face_check) {
            find_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_intersection);
        } else {
            for_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_intersection);
        }

        vector<Ixn> ixns;
        for (int t = 0; t < arcsim::nthreads; t++)
            append(ixns, arcsim::ixns[t]);
        return ixns;
    }

    vector<Ixn> find_intersections(const vector<Face *> &faces, const vector<AccelStruct *> &obs_accs) {

        if (!arcsim::ixns) {
            arcsim::nthreads = omp_get_max_threads();
            arcsim::ixns = new vector<Ixn>[arcsim::nthreads];
        }
        for (int t = 0; t < arcsim::nthreads; t++)
            arcsim::ixns[t].clear();
        for_overlapping_faces(faces, obs_accs, arcsim::thickness, find_face_intersection);
        vector<Ixn> ixns;
        for (int t = 0; t < arcsim::nthreads; t++)
            append(ixns, arcsim::ixns[t]);
        return ixns;
    }

    bool adjacent(const Face *face0, const Face *face1);

    bool intersection_midpoint(const Face *face0, const Face *face1,
                               Bary &b0, Bary &b1);

    bool farthest_points(const Face *face0, const Face *face1, const Vec3 &d,
                         Bary &b0, Bary &b1);

    void find_face_intersection(const Face *face0, const Face *face1) {

        if (adjacent(face0, face1))
            return;

        int t = omp_get_thread_num();

        Bary b0, b1;
        bool is_ixn = intersection_midpoint(face0, face1, b0, b1);
        if (!is_ixn)
            return;

        if (face0 == NULL || face0 == 0) {
            std::cout << "Error: Face0 is bad.\n";
            exit(0);
        } else if (face1 == NULL || face1 == 0) {
            std::cout << "Error: Face1 is bad.\n";
            exit(0);
        }

        // Vec3 oldPos0 = old_pos(face0,b0);
        // Vec3 oldPos1 = old_pos(face1,b1);
        // Vec3 posDiff = oldPos0 - oldPos1;
        //    Vec3 n = normalize(oldPos0 - oldPos1);
        //    farthest_points(face0, face1, n, b0, b1);
        Vec3 n;
        arcsim::ixns[t].push_back(Ixn(face0, b0, face1, b1, n));
    }

    bool adjacent(const Face *face0, const Face *face1) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                if (face0->v[i]->node == face1->v[j]->node)
                    return true;
        return false;
    }

    bool face_plane_intersection(const Face *face, const Face *plane,
                                 Bary &b0, Bary &b1);

    bool intersection_midpoint(const Face *face0, const Face *face1,
                               Bary &b0, Bary &b1) {
        if (norm2(cross(face0->n, face1->n)) < 1e-12)
            return false;
        Bary b00, b01, b10, b11;
        bool ix0 = face_plane_intersection(face0, face1, b00, b01),
                ix1 = face_plane_intersection(face1, face0, b10, b11);
        if (!ix0 || !ix1)
            return false;
        int axis = major_axis(cross(face0->n, face1->n));
        double a00 = pos(face0, b00)[axis], a01 = pos(face0, b01)[axis],
                a10 = pos(face1, b10)[axis], a11 = pos(face1, b11)[axis];
        double amin = std::max(std::min(a00, a01), std::min(a10, a11)),
                amax = std::min(std::max(a00, a01), std::max(a10, a11)),
                amid = (amin + amax) / 2;
        if (amin > amax)
            return false;
        b0 = (a01 == a00) ? b00 : b00 + (amid - a00) / (a01 - a00) * (b01 - b00);
        b1 = (a11 == a10) ? b10 : b10 + (amid - a10) / (a11 - a10) * (b11 - b10);
        return true;
    }

    bool face_plane_intersection(const Face *face, const Face *plane,
                                 Bary &b0, Bary &b1) {
        const Vec3 &x0 = plane->v[0]->node->x, &n = plane->n;
        double h[3];
        int sign_sum = 0;
        for (int v = 0; v < 3; v++) {
            h[v] = dot(face->v[v]->node->x - x0, n);
            sign_sum += sgn(h[v]);
        }
        if (sign_sum == -3 || sign_sum == 3)
            return false;
        int v0 = -1;
        for (int v = 0; v < 3; v++)
            if (sgn(h[v]) == -sign_sum)
                v0 = v;
        double t0 = h[v0] / (h[v0] - h[NEXT(v0)]), t1 = h[v0] / (h[v0] - h[PREV(v0)]);
        b0[v0] = 1 - t0;
        b0[NEXT(v0)] = t0;
        b0[PREV(v0)] = 0;
        b1[v0] = 1 - t1;
        b1[PREV(v0)] = t1;
        b1[NEXT(v0)] = 0;
        return true;
    }

    bool get_line_of_intersection(const Face *face0, const Face *face1, Line3 &line) {
        if (norm2(cross(face0->n, face1->n)) < 1e-12)
            return false;

        Bary b00, b01, b10, b11;
        bool ix0 = face_plane_intersection(face0, face1, b00, b01),
                ix1 = face_plane_intersection(face1, face0, b10, b11);
        if (!ix0 || !ix1)
            return false;

        const Vec3 &p00 = pos(face0, b00);
        const Vec3 &p01 = pos(face0, b01);
        const Vec3 &p10 = pos(face1, b10);
        const Vec3 &p11 = pos(face1, b11);

        int axis = major_axis(cross(face0->n, face1->n));
        double a00 = p00[axis], a01 = p01[axis],
                a10 = p10[axis], a11 = p11[axis];

        // If there is no overlapping between the two lines of intersection,
        if (std::max(a00, a01) < std::min(a10, a11) || std::max(a10, a11) < std::min(a00, a01))
            return false;

        // Sort to pick the overlapping region.
        std::vector<Vec3> pts;
        pts.push_back(p00);
        pts.push_back(p01);
        pts.push_back(p10);
        pts.push_back(p11);
        std::sort(pts.begin(), pts.end(), cmpOneAxis(axis));

        line = std::make_pair(pts[2], pts[1]);
        return true;
    }

    int major_axis(const Vec3 &v) {
        return (abs(v[0]) > abs(v[1]) && abs(v[0]) > abs(v[2])) ? 0
                                                                : (abs(v[1]) > abs(v[2])) ? 1 : 2;
    }

    double vf_clear_distance(const Face *face0, const Face *face1, const Vec3 &d,
                             double last_dist, Bary &b0, Bary &b1);

    double ee_clear_distance(const Face *face0, const Face *face1, const Vec3 &d,
                             double last_dist, Bary &b0, Bary &b1);

    bool farthest_points(const Face *face0, const Face *face1, const Vec3 &d,
                         Bary &b0, Bary &b1) {
        double last_dist = 0;
        last_dist = vf_clear_distance(face0, face1, d, last_dist, b0, b1);
        last_dist = vf_clear_distance(face1, face0, -d, last_dist, b1, b0);
        last_dist = ee_clear_distance(face0, face1, d, last_dist, b0, b1);
        return last_dist > 0;
    }

    double vf_clear_distance(const Face *face0, const Face *face1, const Vec3 &d,
                             double last_dist, Bary &b0, Bary &b1) {
        for (int v = 0; v < 3; v++) {
            const Vec3 &xv = face0->v[v]->node->x, &x0 = face1->v[0]->node->x,
                    &x1 = face1->v[1]->node->x, &x2 = face1->v[2]->node->x;
            const Vec3 &n = face1->n;
            double h = dot(xv - x0, n), dh = dot(d, n);
            if (h * dh >= 0)
                continue;
            double a0 = stp(x2 - x1, xv - x1, d),
                    a1 = stp(x0 - x2, xv - x2, d),
                    a2 = stp(x1 - x0, xv - x0, d);
            if (a0 <= 0 || a1 <= 0 || a2 <= 0)
                continue;
            double dist = -h / dh;
            if (dist > last_dist) {
                last_dist = dist;
                b0 = Bary(0);
                b0[v] = 1;
                b1[0] = a0 / (a0 + a1 + a2);
                b1[1] = a1 / (a0 + a1 + a2);
                b1[2] = a2 / (a0 + a1 + a2);
            }
        }
        return last_dist;
    }

    double ee_clear_distance(const Face *face0, const Face *face1, const Vec3 &d,
                             double last_dist, Bary &b0, Bary &b1) {
        for (int e0 = 0; e0 < 3; e0++) {
            for (int e1 = 0; e1 < 3; e1++) {
                const Vec3 &x00 = face0->v[e0]->node->x,
                        &x01 = face0->v[NEXT(e0)]->node->x,
                        &x10 = face1->v[e1]->node->x,
                        &x11 = face1->v[NEXT(e1)]->node->x;
                Vec3 n = cross(normalize(x01 - x00), normalize(x11 - x10));
                double h = dot(x00 - x10, n), dh = dot(d, n);
                if (h * dh >= 0)
                    continue;
                double a00 = stp(x01 - x10, x11 - x10, d),
                        a01 = stp(x11 - x10, x00 - x10, d),
                        a10 = stp(x01 - x00, x11 - x00, d),
                        a11 = stp(x10 - x00, x01 - x00, d);
                if (a00 * a01 <= 0 || a10 * a11 <= 0)
                    continue;
                double dist = -h / dh;
                if (dist > last_dist) {
                    last_dist = dist;
                    b0 = Bary(0);
                    b0[e0] = a00 / (a00 + a01);
                    b0[NEXT(e0)] = a01 / (a00 + a01);
                    b1 = Bary(0);
                    b1[e1] = a10 / (a10 + a11);
                    b1[NEXT(e1)] = a11 / (a10 + a11);
                }
            }
        }
        return last_dist;
    }

    struct SeparationOpt : public NLConOpt {
        const vector<Ixn> &ixns;
        vector<Node *> nodes;
        double inv_m;

        SeparationOpt(const vector<Ixn> &ixns) : ixns(ixns), inv_m(0) {
            for (int i = 0; i < ixns.size(); i++) {
                if (is_free(ixns[i].f0))
                    for (int v = 0; v < 3; v++)
                        include(ixns[i].f0->v[v]->node, nodes);
                if (is_free(ixns[i].f1))
                    for (int v = 0; v < 3; v++)
                        include(ixns[i].f1->v[v]->node, nodes);
            }
            nvar = nodes.size() * 3;
            ncon = ixns.size();
            for (int n = 0; n < nodes.size(); n++)
                inv_m += 1 / nodes[n]->a;
        }

        void initialize(double *x) const;

        double objective(const double *x) const;

        void obj_grad(const double *x, double *grad) const;

        double constraint(const double *x, int j, int &sign) const;

        void con_grad(const double *x, int j, double factor, double *grad) const;

        void finalize(const double *x) const;
    };

    void solve_ixns(const vector<Ixn> &ixns) {
        augmented_lagrangian_method(SeparationOpt(ixns));
    }

    void SeparationOpt::initialize(double *x) const {
        for (int n = 0; n < nodes.size(); n++)
            set_subvec(x, n, nodes[n]->x);
    }

    double SeparationOpt::objective(const double *x) const {
        double f = 0;
        for (int n = 0; n < nodes.size(); n++) {
            const Node *node = nodes[n];
            Vec3 dx = get_subvec(x, n) - arcsim::xold[get_index(node, *arcsim::meshes)];
            f += inv_m * node->a * dot(dx, dx) / 2;
        }
        return f;
    }

    void SeparationOpt::obj_grad(const double *x, double *grad) const {
        for (int n = 0; n < nodes.size(); n++) {
            const Node *node = nodes[n];
            Vec3 dx = get_subvec(x, n) - arcsim::xold[get_index(node, *arcsim::meshes)];
            set_subvec(grad, n, inv_m * node->a * dx);
        }
    }

    double SeparationOpt::constraint(const double *x, int j, int &sign) const {
        const Ixn &ixn = ixns[j];
        sign = 1;
        double c = -arcsim::thickness;
        for (int v = 0; v < 3; v++) {
            int i0 = find(ixn.f0->v[v]->node, nodes),
                    i1 = find(ixn.f1->v[v]->node, nodes);
            Vec3 x0 = (i0 != -1) ? get_subvec(x, i0) : ixn.f0->v[v]->node->x,
                    x1 = (i1 != -1) ? get_subvec(x, i1) : ixn.f1->v[v]->node->x;
            c += ixn.b0[v] * dot(ixn.n, x0);
            c -= ixn.b1[v] * dot(ixn.n, x1);
        }
        return c;
    }

    void SeparationOpt::con_grad(const double *x, int j, double factor,
                                 double *grad) const {
        const Ixn &ixn = ixns[j];
        for (int v = 0; v < 3; v++) {
            int i0 = find(ixn.f0->v[v]->node, nodes),
                    i1 = find(ixn.f1->v[v]->node, nodes);
            if (i0 != -1)
                add_subvec(grad, i0, factor * ixn.b0[v] * ixn.n);
            if (i1 != -1)
                add_subvec(grad, i1, -factor * ixn.b1[v] * ixn.n);
        }
    }

    void SeparationOpt::finalize(const double *x) const {
        for (int n = 0; n < nodes.size(); n++)
            nodes[n]->x = get_subvec(x, n);
    }
}