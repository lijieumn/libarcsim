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
#include "collision.hpp"
#include "dynamicremesh.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "remesh.hpp"
#include "tensormax.hpp"
#include "timer.hpp"
#include "util.hpp"
#include "io.hpp"
#include <algorithm>
#include <cstdlib>
#include <map>
using namespace std;

namespace arcsim {
    static const bool verbose = false;

    static Cloth::Remeshing *remeshing;
    static bool plasticity;

    static CuttingPlaneSet planeSet;
    static vector<Plane> planes;


    void setPlaneSet(CuttingPlaneSet &planeSet) {
        arcsim::planeSet = planeSet;
    }

    void setPlanes(vector<Plane> &planes) {
        arcsim::planes = planes;
    }

    void create_vert_sizing(Mesh &mesh);

    void destroy_vert_sizing(Mesh &mesh);

    Sizing &operator+=(Sizing &s1, const Sizing &s2) {
        s1.M += s2.M;
        return s1;
    }

    Sizing operator+(const Sizing &s1, const Sizing &s2) {
        Sizing s = s1;
        s += s2;
        return s;
    }

    Sizing &operator*=(Sizing &s, double a) {
        s.M *= a;
        return s;
    }

    Sizing operator*(const Sizing &s, double a) {
        Sizing s2 = s;
        s2 *= a;
        return s2;
    }

    Sizing operator*(double a, const Sizing &s) {
        return s * a;
    }

    Sizing &operator/=(Sizing &s, double a) {
        s.M /= a;
        return s;
    }

    Sizing operator/(const Sizing &s, double a) {
        Sizing s2 = s;
        s2 /= a;
        return s2;
    }

    double norm2(const Vec2 &u, const Sizing &s) {
        return dot(u, s.M * u);
    }

    double norm(const Vec2 &u, const Sizing &s) {
        return sqrt(std::max(norm2(u, s), 0.));
    }

    Mat2x2 mean(const Sizing &s) {
        return s.M;
    }

    template<int n>
    Mat<n, n> sqrt(const Mat<n, n> &A);

// The algorithm

    bool
    fix_up_mesh(vector<Face *> &active, Mesh &mesh, const vector<AccelStruct *> &obs_accs, vector<Edge *> *edges = 0);

    bool split_worst_edge(Mesh &mesh, const vector<AccelStruct *> &obs_accs);

    typedef bool (*CollapsableCallback)(const Node *node);

    bool normal_callback(const Node *node) {
        return true;
    }

    bool VF_callback(const Node *node) {
        return node->temp;
    }

    bool EE_callback(const Node *node) {
        return node->temp2;
    }

    void print_op(RemeshOp &op);

    bool improve_some_face(vector<Face *> &active, Cloth &cloth, CollapsableCallback callback, bool check_collision,
                           const vector<AccelStruct *> &obs_accs);


//std::vector<ArgusImpact> norefine_impacts(const vector<Impact>& impacts, Cloth &cloth, vector<AccelStruct*> obs_accs);

//std::vector<ArgusImpact> norefine_face_impacts(const vector<Impact>& impacts, Cloth &cloth, vector<AccelStruct*> obs_accs);

    std::vector<ArgusImpact>
    refine_impacts(const vector<Impact> &impacts, Cloth &cloth, vector<AccelStruct *> obs_accs);

//std::vector<ArgusImpact> refine_face_impacts(const vector<Impact>& impacts, Cloth &cloth, vector<AccelStruct*> obs_accs);

//std::vector<ArgusImpact> refine_edge_impacts(const vector<Impact>& impacts, Cloth &cloth, vector<AccelStruct*> obs_accs);

    void collapse_cloth_impact_verts(Cloth &cloth, vector<AccelStruct *> &obs_accs);

    void collapse_cloth_impact_edges(Cloth &cloth, vector<AccelStruct *> &obs_accs);

    void set_identity_sizing(Mesh &mesh);

    void static_remesh(Cloth &cloth, const vector<AccelStruct *> &obs_accs, bool collapse) {
        arcsim::remeshing = &cloth.remeshing;
        Mesh &mesh = cloth.mesh;
        for (int v = 0; v < mesh.verts.size(); v++) {
            Sizing *sizing = new Sizing;
            sizing->M = Mat2x2(1.f / sq(remeshing->size_min));
            mesh.verts[v]->sizing = sizing;
        }
        while (split_worst_edge(mesh, obs_accs));
        vector<Face *> active = mesh.faces;
        if (collapse) {
            while (improve_some_face(active, cloth, normal_callback, false, obs_accs));
        }
        for (int v = 0; v < mesh.verts.size(); v++)
            delete mesh.verts[v]->sizing;
        update_indices(mesh);
        compute_ms_data(mesh);
        compute_masses(cloth);
    }

    RemeshOp flip_edges(vector<Face *> &active, Mesh &mesh);


/*
std::vector<ArgusImpact> collision_fixed(Cloth &cloth, std::vector<Mesh*> obs_meshes, const vector<Plane> &planes,
					  const std::vector<Impact>& impacts, bool plasticity){

	::remeshing = &cloth.remeshing;
    ::plasticity = plasticity;
    Mesh &mesh = cloth.mesh;
    set_identity_sizing(mesh);
    // vector<Face*> active = mesh.faces;
    vector<AccelStruct*> obs_accs = create_accel_structs(obs_meshes, false);

	//std::vector<ArgusImpact> argus_vf_impacts = norefine_face_impacts(impacts,cloth,obs_accs);
	std::vector<ArgusImpact> argus_impacts = norefine_impacts(impacts,cloth,obs_accs);


	set_identity_sizing(mesh);

	update_indices(mesh);
    compute_ms_data(mesh);
    compute_masses(cloth);

	return argus_vf_impacts;

}*/

    std::vector<ArgusImpact> collision_refine(Cloth &cloth, std::vector<Mesh *> obs_meshes, const vector<Plane> &planes,
                                              const std::vector<Impact> &impacts, bool plasticity) {

        arcsim::remeshing = &cloth.remeshing;
        arcsim::plasticity = plasticity;
        Mesh &mesh = cloth.mesh;
        set_identity_sizing(mesh);
        // vector<Face*> active = mesh.faces;
        vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, false);

        std::vector<ArgusImpact> argusImpacts = refine_impacts(impacts, cloth, obs_accs);


        set_identity_sizing(mesh);
        destroy_accel_structs(obs_accs);
        destroy_vert_sizing(mesh);

        update_indices(mesh);
        compute_ms_data(mesh);
        compute_masses(cloth);

        return argusImpacts;

    }

    void collision_coarsen(Cloth &cloth, std::vector<Mesh *> obs_meshes, const vector<Plane> &planes, bool plasticity) {
        // ::remeshing = &cloth.remeshing;
        //    ::plasticity = plasticity;
        //    Mesh &mesh = cloth.mesh;
        //    create_vert_sizing(mesh, planes);
        //    vector<Face*> active = mesh.faces;

        // set_identity_sizing(mesh);

        //    vector<AccelStruct*> obs_accs = create_accel_structs(obs_meshes, false);
        //    // VERT-FACE COARSEN
        //    collapse_cloth_impact_verts(cloth, obs_accs);
        //    // EDGE-EDGE COARSEN
        // collapse_cloth_impact_edges(cloth, obs_accs);

        // destroy_accel_structs(obs_accs);
        // destroy_vert_sizing(mesh);
        // update_indices(mesh);
        //    compute_ms_data(mesh);
        //    compute_masses(cloth);

    }

// void test_eigen_decomposition(Mat2x2 &A) {
//     Eig<2> eig1 = eigen_decomposition(A);
//     Eig<2> eig2 = eigen_decomposition_new(A);
//     cout << "A:" << endl << A << endl;
//     cout << "l1:" << eig1.l << endl;
//     cout << "Q1:" << eig1.Q << endl;
//     cout << "----" << endl;
//     cout << "l2:" << eig2.l << endl;
//     cout << "Q2:" << eig2.Q << endl;
//     cout << "----" << endl;
//     cout << "diff l:" << norm(eig1.l - eig2.l) << endl;
//     cout << "diff Q:" << norm_F(eig1.Q - eig2.Q) << endl;
//     cout << endl;
// }

    void dynamic_remesh(Cloth &cloth, bool plasticity, const vector<AccelStruct *> &obs_accs, bool collapse) {
        arcsim::remeshing = &cloth.remeshing;
        arcsim::plasticity = plasticity;
        Mesh &mesh = cloth.mesh;
        VisualDebugger *vd = VisualDebugger::getInstance();

        vector<CuttingPlanes> allPlanes[3];
        allPlanes[0] = arcsim::planeSet.nodePlanes;
        allPlanes[1] = arcsim::planeSet.edgePlanes;
        allPlanes[2] = arcsim::planeSet.facePlanes;
        char keys[] = {'N', 'E', 'F'};
        for (int i = 0; i < 3; i++) {
            for (int s = 0; s < allPlanes[i].size(); s++) {
                CuttingPlanes &planes = allPlanes[i][s];
                for (int p = 0; p < planes.size(); p++) {
                    CuttingPlane &plane = planes[p];
                    Vec3 x = plane.pos;
                    vd->addVisualPoint3(x, Vec3(1, 0, 0), keys[i]);
                    for (int d = 0; d < plane.directions.size(); d++) {
                        Vec3 dir = plane.directions[d];
                        vd->addVisualLine3(x, x + 1e-1 * dir, Vec3(1, 0, 0), keys[i]);
                    }
                }
            }
        }
        create_vert_sizing(mesh);
        vector<Face *> active = mesh.faces;
        fix_up_mesh(active, mesh, obs_accs);
        while (split_worst_edge(mesh, obs_accs));
        active = mesh.faces;
        if (collapse) {
            while (improve_some_face(active, cloth, normal_callback, false, obs_accs));
        }
        for (int i = 0; i < cloth.mesh.verts.size(); i++) {
            Vert *vert = cloth.mesh.verts[i];
            Mat2x2 M = vert->sizing->M;
            Eig<2> eig = eigen_decomposition(M);
            Mat2x2 sqrtM = sqrt(M);
            Mat2x2 diff = sqrtM * sqrtM - M;
            // if (norm_F(diff) > 1e-10) {
            //     cout << "M:" << M << endl;
            //     cout << "sqrtM:" << sqrtM << endl;
            //     cout << "diff:" << diff << endl;
            // }
            if (eig.l[0] < 0 || eig.l[1] < 0) {
                continue;
            }
            vd->addVisualLine2(vert->u, vert->u + eig.Q.col(0) * std::sqrt(1 / eig.l[0]) * .1, Vec3(0, 0, 1), 'm');
            vd->addVisualLine2(vert->u, vert->u - eig.Q.col(0) * std::sqrt(1 / eig.l[0]) * .1, Vec3(0, 0, 1), 'm');
            vd->addVisualLine2(vert->u, vert->u + eig.Q.col(1) * std::sqrt(1 / eig.l[1]) * .1, Vec3(0, 0, 1), 'm');
            vd->addVisualLine2(vert->u, vert->u - eig.Q.col(1) * std::sqrt(1 / eig.l[1]) * .1, Vec3(0, 0, 1), 'm');
        }
        destroy_vert_sizing(mesh);
        update_indices(mesh);
        compute_ms_data(mesh);
        compute_masses(cloth);
    }


    void collapse_edges(Cloth &cloth, const std::vector<Plane> &planes, const vector<AccelStruct *> &obs_accs) {
        // create_vert_sizing(cloth.mesh, planes);
        // vector<Face*> active = cloth.mesh.faces;
        // while (improve_some_face(active, cloth, normal_callback, true, obs_accs));
        // destroy_vert_sizing(cloth.mesh);
        // update_indices(cloth.mesh);
        // compute_ms_data(cloth.mesh);
        // compute_masses(cloth);
    }

// Sizing

    void set_identity_sizing(Mesh &mesh) {

        for (int v = 0; v < mesh.verts.size(); v++) {
            Sizing *sizing = new Sizing;
            sizing->M = Mat2x2(1);
            mesh.verts[v]->sizing = sizing;
        }

    }


    double angle(const Vec3 &n1, const Vec3 &n2) {
        return acos(clamp(dot(n1, n2), -1., 1.));
    }

    template<int n>
    Mat<n, n> sqrt(const Mat<n, n> &A) {
        Eig<n> eig = eigen_decomposition(A);
        for (int i = 0; i < n; i++)
            eig.l[i] = eig.l[i] >= 0 ? std::sqrt(eig.l[i]) : -std::sqrt(-eig.l[i]);
        return eig.Q * diag(eig.l) * eig.Q.t();
    }

    template<int n>
    Mat<n, n> pos(const Mat<n, n> &A) {
        Eig<n> eig = eigen_decomposition(A);
        for (int i = 0; i < n; i++)
            eig.l[i] = std::max(eig.l[i], 0.);
        return eig.Q * diag(eig.l) * eig.Q.t();
    }

    Mat2x2 perp(const Mat2x2 &A) {
        return Mat2x2(Vec2(A(1, 1), -A(1, 0)),
                      Vec2(-A(0, 1), A(0, 0)));
    }

    Mat2x2 compression_metric(const Mat2x2 &e, const Mat2x2 &S2, double c) {
        Mat2x2 D = e.t() * e - 4 * sq(c) * perp(S2) * arcsim::magic.rib_stiffening;
        return pos(-e + sqrt(D)) / (2 * sq(c));
    }

    Mat2x2 obstacle_metric(const Face *face) {
        // tensor_max (const std::vector<Mat2x2> &Ms)
        if (arcsim::magic.use_proximity_metric) {
            vector<Mat2x2> Ms;
            vector<pair<Vec3, Vec3> > ps;
            for (int v = 0; v < 3; v++) {
                CuttingPlanes planes = planeSet.nodePlanes[face->v[v]->node->index];
                for (int p = 0; p < planes.size(); p++) {
                    CuttingPlane plane = planes[p];
                    assert(plane.directions.size() == 1);
                    ps.push_back(make_pair(plane.pos, plane.directions[0]));
                }
            }
            for (int e = 0; e < 3; e++) {
                CuttingPlanes planes = planeSet.edgePlanes[face->adje[e]->index];
                for (int p = 0; p < planes.size(); p++) {
                    CuttingPlane plane = planes[p];
                    assert(plane.directions.size() != 0);
                    assert(plane.directions.size() == 1 || plane.directions.size() == 2);
                    for (int d = 0; d < plane.directions.size(); d++) {
                        ps.push_back(make_pair(plane.pos, plane.directions[d]));
                    }
                }
            }
            CuttingPlanes planes = planeSet.facePlanes[face->index];
            for (int p = 0; p < planes.size(); p++) {
                CuttingPlane plane = planes[p];
                for (int d = 0; d < plane.directions.size(); d++) {
                    ps.push_back(make_pair(plane.pos, plane.directions[d]));
                }
            }
            for (int p = 0; p < ps.size(); p++) {
                double h[3];
                double max_h = 0;
                for (int v = 0; v < 3; v++) {
                    h[v] = dot(face->v[v]->node->x - ps[p].first, ps[p].second);
                    if (fabs(h[v]) > max_h) {
                        max_h = fabs(h[v]);
                    }
                }
                Vec2 dh = derivative(h[0], h[1], h[2], face);
                Ms.push_back(outer(dh, dh) / sq(max_h));
            }

            return tensor_max(Ms);
        } else {
            Mat2x2 o = Mat2x2(0);
            for (int v = 0; v < 3; v++) {
                Plane p = planes[face->v[v]->node->index];
                if (norm2(p.second) == 0)
                    continue;
                double h[3];
                for (int v1 = 0; v1 < 3; v1++)
                    h[v1] = dot(face->v[v1]->node->x - p.first, p.second);
                Vec2 dh = derivative(h[0], h[1], h[2], face);
                o += outer(dh, dh) / sq(h[v]);
            }
            return o / 3.;
        }
    }

    Sizing compute_face_sizing(const Face *face) {
        Sizing s;
        Mat2x2 Sp = curvature<PS>(face);
        Mat2x2 Sw1 = curvature<WS>(face);
        Mat3x2 Sw2 = derivative(face->v[0]->node->n, face->v[1]->node->n,
                                face->v[2]->node->n, face);
        Mat2x2 Mcurvp = !arcsim::plasticity ? Mat2x2(0)
                                      : (Sp.t() * Sp) / sq(remeshing->refine_angle);
        Mat2x2 Mcurvw1 = (Sw1.t() * Sw1) / sq(remeshing->refine_angle);
        Mat2x2 Mcurvw2 = (Sw2.t() * Sw2) / sq(remeshing->refine_angle);
        Mat3x2 V = derivative(face->v[0]->node->v, face->v[1]->node->v,
                              face->v[2]->node->v, face);
        Mat2x2 Mvel = (V.t() * V) / sq(remeshing->refine_velocity);
        Mat3x2 F = derivative(face->v[0]->node->x, face->v[1]->node->x,
                              face->v[2]->node->x, face);
        // Mat2x2 Mcomp = compression_metric(F.t()*F)
        //                / sq(remeshing->refine_compression);
        Mat2x2 Mcomp = compression_metric(F.t() * F - Mat2x2(1), Sw2.t() * Sw2,
                                          remeshing->refine_compression);
        // Mat2x2 Mobs = (planes.empty()) ? Mat2x2(0) : obstacle_metric(face, planeSet);
        Mat2x2 Mobs = obstacle_metric(face);
        // Mat2x2 Mcloth = (cloth_planes.empty()) ? Mat2x2(0) : obstacle_metric(face, cloth_planes);
        Mat2x2 Mcloth = Mat2x2(0);

        // impacts related sizing
        // Mat2x2 Mimp = face->hasImpact() ? Mat2x2(1)*3.0/sq(arcsim::magic.merge_radius): Mat2x2(0);
        //Mat2x2 Mimp = Mat2x2(0);

        vector<Mat2x2> Ms(7);
        Ms[0] = Mcurvp;
        Ms[1] = Mcurvw1;
        Ms[2] = Mcurvw2;
        Ms[3] = Mvel;
        Ms[4] = Mcomp;
        Ms[5] = Mobs;
        Ms[6] = Mcloth;

        s.M = arcsim::magic.combine_tensors ? tensor_max(Ms)
                                      : Ms[0] + Ms[1] + Ms[2] + Ms[3] + Ms[4] + Ms[5] + Ms[6];
        Eig<2> eig = eigen_decomposition(s.M);
        for (int i = 0; i < 2; i++)
            eig.l[i] = clamp(eig.l[i],
                             1.f / sq(remeshing->size_max),
                             1.f / sq(remeshing->size_min));
        double lmax = std::max(eig.l[0], eig.l[1]);
        double lmin = lmax * sq(remeshing->aspect_min);
        for (int i = 0; i < 2; i++)
            if (eig.l[i] < lmin)
                eig.l[i] = lmin;
        s.M = eig.Q * diag(eig.l) * eig.Q.t();
        return s;
    }

    Sizing compute_vert_sizing(const Vert *vert,
                               const map<Face *, Sizing> &face_sizing) {
        Sizing sizing;
        for (int f = 0; f < vert->adjf.size(); f++) {
            const Face *face = vert->adjf[f];
            sizing += face->a / 3. * face_sizing.find((Face *) face)->second;
        }
        sizing /= vert->a;
        return sizing;
    }

// Cache

    void create_vert_sizing(Mesh &mesh) {
        map<Face *, Sizing> face_sizing;
        for (int f = 0; f < mesh.faces.size(); f++)
            face_sizing[mesh.faces[f]] = compute_face_sizing(mesh.faces[f]);
        for (int v = 0; v < mesh.verts.size(); v++)
            mesh.verts[v]->sizing =
                    new Sizing(compute_vert_sizing(mesh.verts[v], face_sizing));
    }

    void destroy_vert_sizing(Mesh &mesh) {
        for (int v = 0; v < mesh.verts.size(); v++)
            delete mesh.verts[v]->sizing;
    }

    double edge_metric(const Vert *vert0, const Vert *vert1) {
        if (!vert0 || !vert1)
            return 0;
        Vec2 du = vert0->u - vert1->u;
        return std::sqrt((norm2(du, *vert0->sizing) + norm2(du, *vert1->sizing)) / 2.);
    }

    double edge_metric(const Edge *edge) {
        return std::max(edge_metric(edge_vert(edge, 0, 0), edge_vert(edge, 0, 1)),
                   edge_metric(edge_vert(edge, 1, 0), edge_vert(edge, 1, 1)));
        // return (edge->adjf[0] && edge->adjf[1]) ? m/2 : m;
    }

// Helpers

    template<typename T>
    void include_all(const vector<T> &u, vector<T> &v) { for (int i = 0; i < u.size(); i++) include(u[i], v); }

    template<typename T>
    void exclude_all(const vector<T> &u, vector<T> &v) { for (int i = 0; i < u.size(); i++) exclude(u[i], v); }

    template<typename T>
    void set_null_all(const vector<T> &u, vector<T> &v) { for (int i = 0; i < u.size(); i++) exclude(u[i], v); }

    void update_active(const RemeshOp &op, vector<Node *> &active) {
        exclude_all(op.removed_nodes, active);
        include_all(op.added_nodes, active);
    }

    void update_active(const RemeshOp &op, vector<Face *> &active) {
        exclude_all(op.removed_faces, active);
        include_all(op.added_faces, active);
    }

    void update_active(const vector<RemeshOp> &ops, vector<Face *> &active) {
        for (int i = 0; i < ops.size(); i++)
            update_active(ops[i], active);
    }

// Fixing-upping

    Vert *most_valent_vert(const vector<Face *> &faces);

// Vert *farthest_neighbor (const Vert *vert);

    bool fix_up_mesh(vector<Face *> &active, Mesh &mesh, const vector<AccelStruct *> &obs_accs, vector<Edge *> *edges) {
        RemeshOp flip_ops = argus_based_flip_edges(active, mesh, obs_accs);
        update_active(flip_ops, active);
        if (edges)
            set_null_all(flip_ops.removed_edges, *edges);
        flip_ops.done();
        return !flip_ops.empty();
    }

    RemeshOp flip_some_edges(vector<Face *> &active, Mesh &mesh);

    RemeshOp flip_edges(vector<Face *> &active, Mesh &mesh) {
        RemeshOp ops;
        for (int i = 0; i < 3 * mesh.verts.size(); i++) {// don't loop without bound
            RemeshOp op = flip_some_edges(active, mesh);
            if (op.empty())
                break;
            ops = compose(ops, op);
        }
        return ops;
    }


    RemeshOp argus_based_flip_edges(vector<Face *> &active, Mesh &mesh, const vector<AccelStruct *> &obs_accs);

    vector<Edge *> find_edges_to_flip(const vector<Face *> &active);

    vector<Edge *> independent_edges(const vector<Edge *> &edges);

    bool inverted(const Face *face) { return area(face) < 1e-12; }

    bool degenerate(const Face *face) {
        return aspect(face) < remeshing->aspect_min / 4;
    }

    bool any_inverted(const vector<Face *> faces) {
        for (int i = 0; i < faces.size(); i++) if (inverted(faces[i])) return true;
        return false;
    }

    bool any_degenerate(const vector<Face *> faces) {
        for (int i = 0; i < faces.size(); i++) if (degenerate(faces[i])) return true;
        return false;
    }

    RemeshOp flip_some_edges(vector<Face *> &active, Mesh &mesh) {
        RemeshOp ops;
        static int n_edges_prev = 0;
        vector<Edge *> edges = independent_edges(find_edges_to_flip(active));
        if (edges.size() == n_edges_prev) // probably infinite loop
            return ops;
        n_edges_prev = edges.size();
        for (int e = 0; e < edges.size(); e++) {
            Edge *edge = edges[e];
            RemeshOp op = flip_edge(edge);
            op.apply(mesh);
            if (any_inverted(op.added_faces)) {
                op.inverse().apply(mesh);
                op.inverse().done();
                continue;
            }
            update_active(op, active);
            ops = compose(ops, op);
        }
        return ops;
    }

    bool should_flip(const Edge *edge);

    vector<Edge *> find_edges_to_flip(const vector<Face *> &active) {
        vector<Edge *> edges;
        map<Edge *, int> eMap;
        int idx = 0;
        for (int f = 0; f < active.size(); f++) {
            for (int e = 0; e < 3; e++) {
                Edge *edge = active[f]->adje[e];
                eMap[edge] = idx++;
            }
        }
        for (auto &search : eMap) {
            edges.push_back(search.first);
        }
        // for (int f = 0; f < active.size(); f++) {
        //     include(active[f]->adje[0], edges);
        //     include(active[f]->adje[1], edges);
        //     include(active[f]->adje[2], edges);
        // }
        vector<Edge *> fedges;
        for (int e = 0; e < edges.size(); e++) {
            Edge *edge = edges[e];
            if (is_seam_or_boundary(edge) || edge->label != 0
                || !should_flip(edge))
                continue;
            fedges.push_back(edge);
        }
        return fedges;
    }

    bool independent(const Edge *edge, const vector<Edge *> &edges) {
        for (int i = 0; i < edges.size(); i++) {
            Edge *edge1 = edges[i];
            if (edge->n[0] == edge1->n[0] || edge->n[0] == edge1->n[1]
                || edge->n[1] == edge1->n[0] || edge->n[1] == edge1->n[1])
                return false;
        }
        return true;
    }

    vector<Edge *> independent_edges(const vector<Edge *> &edges) {
        vector<Edge *> iedges;
        for (int e = 0; e < edges.size(); e++)
            if (independent(edges[e], iedges))
                iedges.push_back(edges[e]);
        return iedges;
    }

    vector<Edge *> independent_edges_new(const vector<Edge *> &edges) {
        vector<Edge *> iedges;
        bool *visited = new bool[edges.size()];
        for (int e = 0; e < edges.size(); e++) {
            visited[e] = false;
        }
        map<Edge *, int> eMap;
        for (int e = 0; e < edges.size(); e++) {
            eMap[edges[e]] = e;
        }
        for (int e = 0; e < edges.size(); e++) {
            if (visited[e]) {
                continue;
            }
            iedges.push_back(edges[e]);
            for (int n = 0; n < 2; n++) {
                Node *node = edges[e]->n[n];
                for (int a = 0; a < node->adje.size(); a++) {
                    Edge *ae = node->adje[a];
                    auto search = eMap.find(ae);
                    if (search != eMap.end()) { //found
                        visited[search->second] = true;;
                    }
                }
            }
        }
        delete[] visited;
        return iedges;
    }

    double cross(const Vec2 &u, const Vec2 &v) { return u[0] * v[1] - u[1] * v[0]; }

// from Bossen and Heckbert 1996
    bool should_flip(const Edge *edge) {
        const Vert *vert0 = edge_vert(edge, 0, 0), *vert1 = edge_vert(edge, 0, 1),
                *vert2 = edge_opp_vert(edge, 0), *vert3 = edge_opp_vert(edge, 1);
        Vec2 x = vert0->u, z = vert1->u, w = vert2->u, y = vert3->u;
        Mat2x2 M = (mean(*vert0->sizing) + mean(*vert1->sizing)
                    + mean(*vert2->sizing) + mean(*vert3->sizing)) / 4.;
        return wedge(z - y, x - y) * dot(x - w, M * (z - w)) + dot(z - y, M * (x - y)) * wedge(x - w, z - w)
               < -arcsim::magic.edge_flip_threshold * (wedge(z - y, x - y) + wedge(x - w, z - w));
    }


// Splitting

    vector<Edge *> find_bad_edges(const Mesh &mesh);

    Sizing mean_vert_sizing(const Vert *vert0, const Vert *vert1);

    Vert *adjacent_vert(const Node *node, const Vert *vert);

    bool split_worst_edge(Mesh &mesh, const vector<AccelStruct *> &obs_accs) {
        vector<Edge *> edges = find_bad_edges(mesh);
        for (int e = 0; e < edges.size(); e++) {
            Edge *edge = edges[e];
            if (!edge) continue;
            Node *node0 = edge->n[0], *node1 = edge->n[1];
            RemeshOp op = split_edge(edge);
            op.apply(mesh);
            for (int v = 0; v < op.added_verts.size(); v++) {
                Vert *vertnew = op.added_verts[v];
                Vert *v0 = adjacent_vert(node0, vertnew),
                        *v1 = adjacent_vert(node1, vertnew);
                vertnew->sizing = new Sizing(mean_vert_sizing(v0, v1));
            }
            set_null_all(op.removed_edges, edges);
            op.done();
            if (verbose)
                cout << "Split " << node0 << " and " << node1 << endl;
            vector<Face *> active = op.added_faces;
            fix_up_mesh(active, mesh, obs_accs, &edges);
        }
        return !edges.empty();
    }

// don't use edge pointer as secondary sort key, otherwise not reproducible
    struct Deterministic_sort {
        inline bool operator()(const std::pair<double, Edge *> &left, const std::pair<double, Edge *> &right) {
            return left.first < right.first;
        }
    } deterministic_sort;

    vector<Edge *> find_bad_edges(const Mesh &mesh) {
        vector<pair<double, Edge *> > edgems;
        for (int e = 0; e < mesh.edges.size(); e++) {
            Edge *edge = mesh.edges[e];
            double m = edge_metric(edge);
            if (m > 1)
                edgems.push_back(make_pair(m, edge));
        }
        sort(edgems.begin(), edgems.end(), deterministic_sort);
        vector<Edge *> edges(edgems.size());
        for (int e = 0; e < edgems.size(); e++)
            edges[e] = edgems[edgems.size() - e - 1].second;
        return edges;
    }

    Sizing mean_vert_sizing(const Vert *vert0, const Vert *vert1) {
        Sizing sizing = *vert0->sizing;
        sizing += *vert1->sizing;
        return sizing / 2.;
    }

    Vert *adjacent_vert(const Node *node, const Vert *vert) {
        const Edge *edge = get_edge(node, vert->node);
        for (int i = 0; i < 2; i++)
            for (int s = 0; s < 2; s++)
                if (edge_vert(edge, s, i) == vert)
                    return edge_vert(edge, s, 1 - i);
        return NULL;
    }


// Splitting Faces (for impact handling)


    bool find_in_mesh2(const Node *node, const Mesh &mesh) {
        for (int n = 0; n < mesh.nodes.size(); n++) {
            if (node == mesh.nodes[n]) {
                return true;
            }
        }
        return false;
    }


// Produces ArgusImpacts from impacts which require no refinement. These are VF collisions where the vertex
// belongs to the cloth and the face belongs to an external object
/*std::vector<ArgusImpact> norefine_face_impacts(const vector<Impact>& impacts, Cloth &cloth, vector<AccelStruct*> obs_accs){

	std::vector<ArgusImpact> argusImpacts;

	// Looping through all the raw impacts and splitting enclosing faces
	// We also generate new ArgusImpacts from the raw impacts
	for(int i = 0; i < impacts.size(); i++){

		// Only considering impacts which are between a vertex and a face
		if( impacts[i].type == Impact::VF ){

			Impact impact = impacts[i];

			// Analyzing the four nodes in the standard impact signature.
			bool firstNodeInMesh = find_in_mesh2(impact.nodes[0],cloth.mesh);
			bool secondNodeInMesh = find_in_mesh2(impact.nodes[1],cloth.mesh);
			bool thirdNodeInMesh = find_in_mesh2(impact.nodes[2],cloth.mesh);
			bool fourthNodeInMesh = find_in_mesh2(impact.nodes[3],cloth.mesh);

			// If vertex in the cloth and face in external object. No refinement needed, simply convert to ArgusImpact as-is
			if( firstNodeInMesh && !secondNodeInMesh && !thirdNodeInMesh && !fourthNodeInMesh ){
				ArgusImpact aImpact;
				aImpact.nodeA = impact.nodes[0];
				aImpact.nodeB = 0;
				aImpact.posB = -impact.w[1]*impact.nodes[1]->x + -impact.w[2]*impact.nodes[2]->x + -impact.w[3]*impact.nodes[3]->x;
				aImpact.velB = -impact.w[1]*impact.nodes[1]->v + -impact.w[2]*impact.nodes[2]->v + -impact.w[3]*impact.nodes[3]->v;
				aImpact.normal = impact.n;
				aImpact.obsNormal = impact.obsNormal;
				aImpact.type = "VF_STD";
				aImpact.debug = impact.debug;
				argusImpacts.push_back( aImpact );

			}
		}
	}

	return argusImpacts;

}
*
*/

    void move_node_vf(Impact &impact, int i) {
        Node *to_move = impact.nodes[i];
        Node *n1 = impact.nodes[1];
        Node *n2 = impact.nodes[2];
        Node *n3 = impact.nodes[3];

        double b1 = -impact.w[1];
        double b2 = -impact.w[2];
        double b3 = -impact.w[3];

        to_move->y = b1 * n1->y + b2 * n2->y + b3 * n3->y;
        to_move->x = b1 * n1->x + b2 * n2->x + b3 * n3->x;
        to_move->x0 = b1 * n1->x0 + b2 * n2->x0 + b3 * n3->x0;
        to_move->v = b1 * n1->v + b2 * n2->v + b3 * n3->v;
        to_move->acceleration = b1 * n1->acceleration + b2 * n2->acceleration + b3 * n3->acceleration;
        // to_move->label = n1->label == n2->label ? (n2->label == n3->label ? n1->label : 0) : 0;

        Node *obsNode = impact.nodes[0];
        Face *obsFace = obsNode->verts[0]->adjf[0];
        impact.inverted = false;
        impact.nodes[0] = to_move;
        impact.nodes[1] = obsNode;
        for (int n = 0; n < 3; n++) {
            if (obsFace->v[n]->node == obsNode) {
                impact.nodes[2] = obsFace->v[(n + 1) % 3]->node;
                impact.nodes[3] = obsFace->v[(n + 2) % 3]->node;
                break;
            }
        }
        impact.w[0] = 1.0;
        impact.w[1] = -1.0;
        impact.w[2] = 0;
        impact.w[3] = 0;
        impact.n = -impact.n;

        // Un-inverting the impact
        impact.inverted = false;
    }

    void move_node_ee(Impact &impact, int i) {
        Node *to_move = impact.nodes[i];
        Node *n1 = impact.nodes[0];
        Node *n2 = impact.nodes[1];

        double b1 = impact.w[0];
        double b2 = impact.w[1];

        to_move->y = b1 * n1->y + b2 * n2->y;
        to_move->x = b1 * n1->x + b2 * n2->x;
        to_move->x0 = b1 * n1->x0 + b2 * n2->x0;
        to_move->v = b1 * n1->v + b2 * n2->v;
        to_move->acceleration = b1 * n1->acceleration + b2 * n2->acceleration;

        Node *obsNode1 = impact.nodes[2];
        Node *obsNode2 = impact.nodes[3];

        Edge *obsEdge = getCommonEdge(obsNode1, obsNode2);
        Face *obsFace = obsEdge->adjf[0];
        impact.nodes[0] = to_move;
        impact.nodes[1] = obsNode1;
        impact.nodes[2] = obsNode2;
        for (int n = 0; n < 3; n++) {
            if (obsFace->v[n]->node != obsNode1 && obsFace->v[n]->node != obsNode2) {
                impact.nodes[3] = obsFace->v[n]->node;
                break;
            }
        }
        impact.w[0] = 1.0;
        impact.w[1] = impact.w[2];
        impact.w[2] = impact.w[3];
        impact.w[3] = 0;
        impact.type = Impact::VF;

        impact.inverted = false;
    }

    void print_impacts(const Impact &impact) {
        cout << "impact:---------------------------" << endl;
        cout << "type:" << impact.type << endl;
        std::cout << "self? " << impact.self << std::endl;
        cout << (impact.inverted ? "inverted" : "not inverted") << endl;
        for (int ii = 0; ii < 4; ii++) {
            if (!impact.nodes[ii] || impact.nodes[ii] == NULL || impact.nodes[ii] == 0) {
                continue;
            }
            cout << "node" << ii << ":" << impact.nodes[ii] << "\t" << impact.w[ii] << "\t" << impact.nodes[ii]->x
                 << "\t" << impact.nodes[ii]->verts[0]->u << endl;
        }
        cout << "matPosA:" << impact.matPosA << endl;
        cout << "matPosB:" << impact.matPosB << endl;
        cout << "normal:" << impact.n << endl;
        cout << "time:" << impact.t << endl;
    }

    void print_argus_impacts(ArgusImpact &aImpact) {
        std::cout << "argusImpact:---------------------\n";
        std::cout << "type: " << aImpact.type << std::endl;
        std::cout << "nodeA: " << aImpact.nodeA << std::endl;
        std::cout << "posA: " << aImpact.nodeA->x << std::endl;
        if (aImpact.nodeB != 0 && aImpact.nodeB != NULL) {
            std::cout << "nodeB: " << aImpact.nodeB << std::endl;
            std::cout << "posB: " << aImpact.nodeB->x << std::endl;
        } else {
            std::cout << "posB: " << aImpact.posB << std::endl;
            std::cout << "velB: " << aImpact.velB << std::endl;
        }
        cout << "normal:" << aImpact.normal << endl;

    }

    RemeshOp refine_edge_impact(Cloth &cloth, Node *nodeA, Node *nodeB, double bA, double bB, Vec2 impMatPos) {

        // Find the common edge (if any) between nodeA and nodeB in the impact signature
        Edge *commonEdge = getCommonEdge(nodeA, nodeB);
        Face *encFace = get_enclosing_face(cloth.mesh, impMatPos);

        // If there is a common edge, then the impact data is still valid, so we can just split the edge
        if (commonEdge != 0 && commonEdge != NULL) {
            if (commonEdge->n[0] == nodeA) {
                return split_edge(commonEdge, bA);
            } else {
                return split_edge(commonEdge, bB);
            }
        }

            // Else find the enclosing face of the impact point
        else if (encFace != 0 && encFace != NULL) {

            Vec3 b = get_barycentric_coords(impMatPos, encFace);

            // If the impact is along edge e12 of face, refine that edge
            if (fabs(b[0]) < 1.e-3 && b[1] > 1.e-3 && b[2] > 1.e-3) {
                commonEdge = getCommonEdge(encFace->v[1]->node, encFace->v[2]->node);
                if (commonEdge->n[0] == encFace->v[1]->node) {
                    return split_edge(commonEdge, b[1]);
                } else {
                    return split_edge(commonEdge, b[2]);
                }
            }

                // If the impact is along edge e02 of face, refine that edge
            else if (fabs(b[1]) < 1.e-3 && b[0] > 1.e-3 && b[2] > 1.e-3) {
                commonEdge = getCommonEdge(encFace->v[0]->node, encFace->v[2]->node);
                if (commonEdge->n[0] == encFace->v[0]->node) {
                    return split_edge(commonEdge, b[0]);
                } else {
                    return split_edge(commonEdge, b[2]);
                }

            }

                // If the impact is along edge e01 of face, refine that edge
            else if (fabs(b[2]) < 1.e-3 && b[0] > 1.e-3 && b[1] > 1.e-3) {
                commonEdge = getCommonEdge(encFace->v[0]->node, encFace->v[1]->node);
                if (commonEdge->n[0] == encFace->v[0]->node) {
                    return split_edge(commonEdge, b[0]);
                } else {
                    return split_edge(commonEdge, b[1]);
                }
            }

                // Else there simply is no edge to refine, so do nothing
            else {
                return RemeshOp();
            }

        } else {
            std::cout
                    << "Error: Tried to refine edge impact but there was no edge nor enclosing face at impact position.\n";
            exit(0);
        }
    }


/*


	// Bary pos of impact (along first edge) in material space
	Vec2 u = bA*nodeA->verts[0]->u + bB*nodeB->verts[0]->u;

	std::cout << "u is: " << u << std::endl;

	// Finding the edge these nodes share in common (if any)
	Edge* commonEdge = getCommonEdge(nodeA,nodeB);

	std::cout << "commonEdge is found. is it 0 or null? " << (commonEdge == 0 || commonEdge == NULL) << std::endl;

	if( commonEdge == 0 || commonEdge == NULL ){
        // cout << "commonEdge not found ---------------!!!!!!!!!!!" << endl;
		// return RemeshOp();
		return RemeshOp();

	}	else {

		return split_edge(commonEdge,0.2);
	}

	std::cout << "about to get enclosing face \n";

	// Finding the enclosing face that the impact pos is within (if any)
	Face* encFace = get_enclosing_face(cloth.mesh,u);

	std::cout << "got enclosing face \n";

	// If a common edge exists, we can use the barycentric coordinate data from the impact signature
	double bary = -1;
	if( commonEdge != 0 && commonEdge != NULL ){
		if( commonEdge->n[0] == nodeA ){
			bary = bA;
		} else {
			bary = bB;
		}

		std::cout << "there was a common edge, so set bary \n";

	}


	// Else, if impact pos is at least inside of some enclosing face
	else if( encFace != 0 && encFace != NULL ){

		std::cout << "there was not a common edge, so getting bary coords from enc face \n";

		Vec3 b = get_barycentric_coords(u,encFace);

		double b0;
		double b1;
		Node* node0;
		Node* node1;

		if( b[0] < b[1] && b[0] < b[2] ){
			node0 = encFace->v[1]->node;
			node1 = encFace->v[2]->node;
			b0 = b[1];
			b1 = b[2];
		} else if( b[1] < b[0] && b[1] < b[2] ){
			node0 = encFace->v[0]->node;
			node1 = encFace->v[2]->node;
			b0 = b[0];
			b1 = b[2];
		} else if( b[2] < b[0] && b[2] < b[1] ){
			node0 = encFace->v[0]->node;
			node1 = encFace->v[1]->node;
			b0 = b[0];
			b1 = b[1];
		} else {
			return RemeshOp();
			exit(0);
		}

		commonEdge = getCommonEdge(node0,node1);

		if( commonEdge->n[0] == node0 ){
			bary = b0;
		} else {
			bary = b1;
		}

	}


	std::cout << "declaring edge split op \n";

	// A splitting (and ArgusImpact generation) only occurs if the barys aren't too small or large
	RemeshOp edgeSplitOp;

	// if( bary > 0.1 && 1.0-bary > 0.1){ //no need to check bary. It's already handled by impacts merging
	int numEdgesBefore = cloth.mesh.edges.size();

	std::cout << "is common edge 0 or null? \n";
	if( commonEdge == 0 || commonEdge == NULL ){
		std::cout << "yes \n";
	} else {
		std::cout << "no \n";
	}

	std::cout << "about to grab common edge -> n[0] \n";

	Node* node0 = commonEdge->n[0];

	std::cout << "about to grab common edge -> n[1] \n";

	Node* node1 = commonEdge->n[1];
	//edgeSplitOp = split_edge(commonEdge);

	std::cout << "setting the edge split op \n";

	edgeSplitOp = split_edge(commonEdge,bary);
	// }

	std::cout << "returning edge split op \n";

	std::cout << std::endl << std::endl;
	return edgeSplitOp;
*/





    bool create_argus_vf_self(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        Vec2 u = -impact.w[1] * impact.nodes[1]->verts[0]->u +
                 -impact.w[2] * impact.nodes[2]->verts[0]->u +
                 -impact.w[3] * impact.nodes[3]->verts[0]->u;

        Face *encFace = get_enclosing_face(cloth.mesh, u);
        Vec3 b = get_barycentric_coords(u, encFace);

        if (fabs(b[0] - 1.f) < 0.1) {
            aImpact.nodeA = impact.nodes[0];
            aImpact.nodeB = encFace->v[0]->node;
            aImpact.normal = impact.n;
            aImpact.debug = impact.debug;
        } else if (fabs(b[1] - 1.f) < 0.1) {
            aImpact.nodeA = impact.nodes[0];
            aImpact.nodeB = encFace->v[1]->node;
            aImpact.normal = impact.n;
            aImpact.debug = impact.debug;
        } else if (fabs(b[2] - 1.f) < 0.1) {
            aImpact.nodeA = impact.nodes[0];
            aImpact.nodeB = encFace->v[2]->node;
            aImpact.normal = impact.n;
            aImpact.debug = impact.debug;
        } else {

            RemeshOp faceSplitOp = split_face(encFace, b[0], b[1], b[2]);
            faceSplitOp.apply(cloth.mesh);

            aImpact.nodeA = impact.nodes[0];
            aImpact.nodeB = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
            aImpact.posB = -impact.w[1] * impact.nodes[1]->x + -impact.w[2] * impact.nodes[2]->x +
                           -impact.w[3] * impact.nodes[3]->x;
            aImpact.velB = -impact.w[1] * impact.nodes[1]->v + -impact.w[2] * impact.nodes[2]->v +
                           -impact.w[3] * impact.nodes[3]->v;
            aImpact.normal = impact.n;
            // aImpact.obsNormal = impact.obsNormal;
            aImpact.debug = impact.debug;

            // Updating sizing data
            set_identity_sizing(cloth.mesh);

            // Updating indices and other data
            update_indices(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);

            Vert *newVert = aImpact.nodeA->verts[0];
            // std::vector<Face*> adjf = newVert -> adjf;
            vector<Face *> activeFaces = cloth.mesh.faces;

            argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);
            // RemeshOp flip_ops = flip_edges(adjf, cloth.mesh);
            // compute_ms_data(cloth.mesh);

            // if (has_proximity(flip_ops.added_faces, obs_accs)) {
            //       if (has_intersection(flip_ops.added_faces, obs_accs)) {
            // 	flip_ops.inverse().apply(cloth.mesh);
            // 	flip_ops.inverse().done();
            // } else {
            // 	flip_ops.done();
            // }

        }

        aImpact.self = true;
        aImpact.inverted = false;

        return true;

    }

    bool create_argus_vf_standard(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        aImpact.nodeA = impact.nodes[0];
        aImpact.nodeB = 0;
        aImpact.posB = -impact.w[1] * impact.nodes[1]->x + -impact.w[2] * impact.nodes[2]->x +
                       -impact.w[3] * impact.nodes[3]->x;
        aImpact.velB = -impact.w[1] * impact.nodes[1]->v + -impact.w[2] * impact.nodes[2]->v +
                       -impact.w[3] * impact.nodes[3]->v;
        aImpact.normal = impact.n;
        aImpact.debug = impact.debug;

        return true;
    }

    bool create_argus_vf_inverted(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {


        //Vec2 u = -impact.w[1]*impact.nodes[1]->verts[0]->u +
        //		 -impact.w[2]*impact.nodes[2]->verts[0]->u +
        //		 -impact.w[3]*impact.nodes[3]->verts[0]->u;

        Vec2 u = impact.matPosB;

        Face *encFace = get_enclosing_face(cloth.mesh, u);
        Vec3 b = get_barycentric_coords(u, encFace);

        // Barycentrically the splitting point shouldn't be too close to a vert
        // in the enclosing face. If it is, then don't add the argus impact
        // TODO: remove this restriction and just rely on flipping
        // if( fabs(b[0]) < 0.1 || fabs(b[1]) < 0.1 || fabs(b[2]) < 0.1 ){
        // 	return false;
        // }

        RemeshOp faceSplitOp = split_face(encFace, b[0], b[1], b[2]);

        faceSplitOp.apply(cloth.mesh);
        aImpact.nodeA = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
        aImpact.nodeB = 0;
        aImpact.posB = impact.nodes[0]->x;
        aImpact.velB = impact.nodes[0]->v;
        aImpact.normal = -impact.n;
        // aImpact.obsNormal = impact.obsNormal;  // should this be negative?
        aImpact.debug = impact.debug;
        aImpact.inverted = true;

        // Updating sizing data
        set_identity_sizing(cloth.mesh);

        // Updating indices and other data
        update_indices(cloth.mesh);
        compute_ms_data(cloth.mesh);
        compute_masses(cloth);

        Vert *newVert = aImpact.nodeA->verts[0];
        // std::vector<Face*> adjf = newVert -> adjf;

        vector<Face *> activeFaces = cloth.mesh.faces;

        argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);
        // RemeshOp flip_ops = flip_edges(adjf, cloth.mesh);
        // compute_ms_data(cloth.mesh);
        // std::vector<Face*> addedFaces = flip_ops.added_faces;


        // if (has_intersection(flip_ops.added_faces, obs_accs)) {
        // 	flip_ops.inverse().apply(cloth.mesh);
        // 	flip_ops.inverse().done();
        // } else {
        // 	flip_ops.done();
        // }

        return true;
    }

    bool create_argus_ee_self(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        {
            // Splitting the first edge
            RemeshOp splitOp0 = refine_edge_impact(cloth, impact.nodes[0], impact.nodes[1], impact.w[0], impact.w[1],
                                                   impact.matPosA);
            splitOp0.apply(cloth.mesh);
            Node *newNode = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
            newNode->temp2 = true;
            Vert *newVert = newNode->verts[0];

            aImpact.nodeA = newNode;
            aImpact.normal = impact.n;
            aImpact.self = true;

            // Updating sizing data
            set_identity_sizing(cloth.mesh);

            // Updating indices and other data
            update_indices(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);

            // Flipping edges to maintain Delaunay triangulation
            // std::vector<Face*> adjf = newVert -> adjf;
            vector<Face *> activeFaces = cloth.mesh.faces;
            argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);

            // RemeshOp flip_ops = flip_edges(adjf, cloth.mesh);
            // compute_ms_data(cloth.mesh);
            // // if (has_proximity(flip_ops.added_faces, obs_accs)) {
            //       if (has_intersection(flip_ops.added_faces, obs_accs)) {
            // 	flip_ops.inverse().apply(cloth.mesh);
            // 	flip_ops.inverse().done();
            // } else {
            //           flip_ops.done();
            // }

            // // Updating indices and other data
            // update_indices(cloth.mesh);
            // compute_ms_data(cloth.mesh);
            // compute_masses(cloth);

        }

        {
            // Splitting the second edge
            RemeshOp splitOp1 = refine_edge_impact(cloth, impact.nodes[2], impact.nodes[3], -impact.w[2], -impact.w[3],
                                                   impact.matPosB);
            splitOp1.apply(cloth.mesh);
            Node *newNode = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
            newNode->temp2 = true;
            Vert *newVert = newNode->verts[0];

            aImpact.nodeB = newNode;

            // Updating sizing data
            set_identity_sizing(cloth.mesh);

            // Updating indices and other data
            update_indices(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);

            // Flipping edges to maintain Delaunay triangulation
            // std::vector<Face*> adjf = newVert -> adjf;
            vector<Face *> activeFaces = cloth.mesh.faces;
            argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);
            // RemeshOp flip_ops = flip_edges(adjf, cloth.mesh);
            // compute_ms_data(cloth.mesh);
            // // if (has_proximity(flip_ops.added_faces, obs_accs)) {
            //       if (has_intersection(flip_ops.added_faces, obs_accs)) {
            // 	flip_ops.inverse().apply(cloth.mesh);
            // 	flip_ops.inverse().done();
            // } else {
            //           flip_ops.done();
            // }

            // // Updating indices and other data
            // update_indices(cloth.mesh);
            // compute_ms_data(cloth.mesh);
            // compute_masses(cloth);
        }

        return true;

    }

    bool create_argus_ee_standard(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        // Try to refine the edge at the impact site
        RemeshOp splitOp0 = refine_edge_impact(cloth, impact.nodes[0], impact.nodes[1], impact.w[0], impact.w[1],
                                               impact.matPosA);

        // If that failed, try to refine the enclosing face
        if (splitOp0.empty()) {

            Face *encFace = get_enclosing_face(cloth.mesh, impact.matPosA);

            if (encFace != 0 && encFace != NULL) {
                Vec3 b = get_barycentric_coords(impact.matPosA, encFace);
                splitOp0 = split_face(encFace, b[0], b[1], b[2]);

            }

        }

        // If some form of refinement succeeded, apply the operation and create the argus impact
        if (!splitOp0.empty()) {

            splitOp0.apply(cloth.mesh);

            Node *newNode = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
            newNode->temp2 = true;

            Vert *newVert = newNode->verts[0];

            // Updating sizing data
            set_identity_sizing(cloth.mesh);

            // Updating indices and other data
            update_indices(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);


            aImpact.nodeA = newNode;
            aImpact.nodeB = 0;
            aImpact.posB = -impact.w[2] * impact.nodes[2]->x + -impact.w[3] * impact.nodes[3]->x;
            aImpact.velB = -impact.w[2] * impact.nodes[2]->v + -impact.w[3] * impact.nodes[3]->v;
            aImpact.normal = impact.n;
            aImpact.debug = impact.debug;

            // Flipping edges to maintain delaunay triangulation
            vector<Face *> activeFaces = cloth.mesh.faces;
            argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);


            // Updating indices and other data
            // update_indices(cloth.mesh);
            // compute_ms_data(cloth.mesh);
            // compute_masses(cloth);


            return true;

        }

        return false;
    }

    bool create_argus_ee_inverted(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        aImpact.inverted = true;
        return false;
    }

    bool create_argus_ve_self(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {


        aImpact.self = true;
        return false;

    }

    bool create_argus_ve_standard(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        aImpact.nodeA = impact.nodes[0];
        aImpact.nodeB = 0;
        aImpact.posB = impact.w[1] * impact.nodes[1]->x + impact.w[2] * impact.nodes[2]->x;
        aImpact.velB = impact.w[1] * impact.nodes[1]->v + impact.w[2] * impact.nodes[2]->v;
        aImpact.normal = impact.n;
        // aImpact.obsNormal = impact.obsNormal;
        aImpact.debug = impact.debug;

        return true;
    }

    bool create_argus_ve_inverted(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        RemeshOp edgeSplitOp = refine_edge_impact(cloth, impact.nodes[1], impact.nodes[2], impact.w[1], impact.w[2],
                                                  impact.matPosB);
        edgeSplitOp.apply(cloth.mesh);

        aImpact.nodeA = cloth.mesh.nodes[cloth.mesh.nodes.size() - 1];
        aImpact.nodeB = 0;
        aImpact.posB = impact.nodes[0]->x;
        aImpact.velB = impact.nodes[0]->v;
        aImpact.normal = -impact.n;
        // aImpact.obsNormal = impact.obsNormal; // should this be negative?
        aImpact.debug = impact.debug;
        aImpact.inverted = true;

        return true;
    }

    bool create_argus_vv_self(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {


        aImpact.self = true;
        return false;

    }

    bool create_argus_vv_standard(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        aImpact.nodeA = impact.nodes[0];
        aImpact.nodeB = 0;
        aImpact.posB = impact.nodes[1]->x;
        aImpact.velB = impact.nodes[1]->v;
        aImpact.normal = impact.n;
        // aImpact.obsNormal = impact.obsNormal;
        aImpact.debug = impact.debug;

        return true;
    }

    bool create_argus_vv_inverted(ArgusImpact &aImpact, Impact &impact, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        aImpact.nodeA = impact.nodes[1];
        aImpact.nodeB = 0;
        aImpact.posB = impact.nodes[0]->x;
        aImpact.velB = impact.nodes[0]->v;
        aImpact.normal = -impact.n;
        // aImpact.obsNormal = impact.obsNormal;
        aImpact.debug = impact.debug;
        aImpact.inverted = true;
        return true;
    }

    void save_mesh(Mesh &mesh) {
        update_indices(mesh);
        compute_ms_data(mesh);
        save_obj(mesh, stringf("%s%06d.obj", "temp/", mesh.sub_step++));
    }

// General refinement, produces Argus impacts
    std::vector<ArgusImpact>
    refine_impacts(const vector<Impact> &impacts, Cloth &cloth, vector<AccelStruct *> obs_accs) {

        std::vector<ArgusImpact> argusImpacts;

        for (int i = 0; i < impacts.size(); i++) {

            ArgusImpact argusImpact;
            Impact rawImpact = impacts[i];
            bool created = false;

            switch (rawImpact.type) {

                // VERTEX-FACE PROXIMITY
                case Impact::VF: {


                    // SELF VERTEX-FACE
                    if (rawImpact.self) {
                        created = create_argus_vf_self(argusImpact, rawImpact, cloth, obs_accs);
                    }
                        // INVERTED VERTEX-FACE
                    else if (rawImpact.inverted) {
                        created = create_argus_vf_inverted(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // STANDARD VERTEX-FACE
                    else {
                        created = create_argus_vf_standard(argusImpact, rawImpact, cloth, obs_accs);
                    }

                    argusImpact.type = ArgusImpact::VF;

                    break;
                }

                    // EDGE-EDGE PROXIMITY
                case Impact::EE: {

                    // SELF EDGE-EDGE
                    if (rawImpact.self) {
                        created = create_argus_ee_self(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // INVERTED EDGE-EDGE
                    else if (rawImpact.inverted) {
                        created = create_argus_ee_inverted(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // STANDARD EDGE-EDGE
                    else {
                        created = create_argus_ee_standard(argusImpact, rawImpact, cloth, obs_accs);
                    }

                    argusImpact.type = ArgusImpact::EE;

                    break;
                }

                    // VERTEX-EDGE PROXIMITY
                case Impact::VE: {

                    // SELF VERTEX-EDGE
                    if (rawImpact.self) {
                        created = create_argus_ve_self(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // INVERTED VERTEX-EDGE
                    else if (rawImpact.inverted) {
                        created = create_argus_ve_inverted(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // STANDARD VERTEX-EDGE
                    else {
                        created = create_argus_ve_standard(argusImpact, rawImpact, cloth, obs_accs);
                    }

                    argusImpact.type = ArgusImpact::VE;

                    break;

                }

                    // VERTEX-VERTEX PROXIMITY
                case Impact::VV: {

                    // SELF VERTEX-VERTEX
                    if (rawImpact.self) {
                        created = create_argus_vv_self(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // INVERTED VERTEX-VERTEX
                    else if (rawImpact.inverted) {
                        created = create_argus_vv_inverted(argusImpact, rawImpact, cloth, obs_accs);
                    }

                        // STANDARD VERTEX-VERTEX
                    else {
                        created = create_argus_vv_standard(argusImpact, rawImpact, cloth, obs_accs);
                    }

                    argusImpact.type = ArgusImpact::VV;

                    break;

                }

                default: {
                    std::cout << "Error: Invalid impact type in refine_impacts.\n";
                    exit(0);
                }

            } // end of switch


            //	if( created && Impact

            if (created) {
                argusImpact.nodeA->freeAge = 0;
                if (argusImpact.nodeB != 0 && argusImpact.nodeB != NULL) {
                    argusImpact.nodeB->freeAge = 0;
                }
                argusImpacts.push_back(argusImpact);
                for (int j = 0; j < rawImpact.adjacentImpacts.size(); j++) {
                    ArgusImpact adjImpact;
                    adjImpact.nodeB = 0;
                    adjImpact.posB = rawImpact.adjacentImpacts[j].obsX;
                    adjImpact.velB = rawImpact.adjacentImpacts[j].obsV;
                    if (rawImpact.adjacentImpacts[j].inverted) {
                        adjImpact.normal = -rawImpact.adjacentImpacts[j].n;
                    } else {
                        adjImpact.normal = rawImpact.adjacentImpacts[j].n;
                    }
                    adjImpact.nodeA = argusImpact.nodeA;
                    argusImpacts.push_back(adjImpact);
                }
            }

        }


        // Updating indices and other data
        update_indices(cloth.mesh);
        compute_ms_data(cloth.mesh);
        compute_masses(cloth);

        return argusImpacts;

    }



// Refining edges (for impact handling)







// Collapsing


    Vec3 preCollapseMomentum(Node *nodeToCollapse, std::vector<Node *> &oneRing) {

        // Computing oneRing momentum before collapse
        oneRing.clear();
        Vec3 pBef = nodeToCollapse->m * nodeToCollapse->v;
        for (int e = 0; e < nodeToCollapse->adje.size(); e++) {
            Edge *edge1 = nodeToCollapse->adje[e];
            Node *node2 = (edge1->n[0] != nodeToCollapse) ? edge1->n[0] : edge1->n[1];
            oneRing.push_back(node2);
            pBef += node2->m * node2->v;
        }

        return pBef;

    }


    Vec3 postCollapseMomentum(std::vector<Node *> &oneRing) {

        Vec3 pAft = Vec3(0, 0, 0);
        for (int i = 0; i < oneRing.size(); i++) {
            pAft += oneRing[i]->m * oneRing[i]->v;
        }

        return pAft;

    }


    void applyEdgeCollapseOp(RemeshOp &op, Cloth &cloth, Node *nodeToCollapse) {
        if (arcsim::magic.conserve_momentum) {

            std::vector<Node *> oneRing;
            Vec3 pBef = preCollapseMomentum(nodeToCollapse, oneRing);
            Vec2 uCollapse = nodeToCollapse->verts[0]->u;
            op.apply(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);
            Vec3 pAft = postCollapseMomentum(oneRing);
            Vec3 dp = pAft - pBef;

            //std::cout << "pBef: " << pBef << std::endl;
            //std::cout << "pAft: " << pAft << std::endl;
            //std::cout << "dp: " << dp << std::endl;

            // Buffering velocities in case this operation is reverted later
            for (int i = 0; i < oneRing.size(); i++) {
                oneRing[i]->vBuffer = oneRing[i]->v;
            }

            Face *face = get_enclosing_face(cloth.mesh, uCollapse);
            Vec3 b = get_barycentric_coords(uCollapse, face);


            double baryMass = b[0] * face->v[0]->node->m + b[1] * face->v[1]->node->m + b[2] * face->v[2]->node->m;
            Vec3 dv0 = b[0] * dp / baryMass;
            Vec3 dv1 = b[1] * dp / baryMass;
            Vec3 dv2 = b[2] * dp / baryMass;

            face->v[0]->node->v -= dv0;
            face->v[1]->node->v -= dv1;
            face->v[2]->node->v -= dv2;

            Vec3 pFix = postCollapseMomentum(oneRing);

            //std::cout << "pFix: " << pFix << std::endl;
            //std::cout << "pFix - pBef: " << pFix - pBef << std::endl;



        } else {
            op.apply(cloth.mesh);
        }
    }

    void applyInverseEdgeCollapseOp(RemeshOp &op, Cloth &cloth) {

        Node *collapsedNode = op.removed_nodes.size() == 1 ? op.removed_nodes[0] : 0;
        if (collapsedNode == 0) {
            std::cout << "Error: Cannot invert an edge collapse op which didn't collapse anything.\n";
            exit(0);
        }

        op.inverse().apply(cloth.mesh);

        if (arcsim::magic.conserve_momentum) {
            for (int e = 0; e < collapsedNode->adje.size(); e++) {
                Edge *edge1 = collapsedNode->adje[e];
                Node *node2 = (edge1->n[0] != collapsedNode) ? edge1->n[0] : edge1->n[1];
                node2->v = node2->vBuffer;
            }
        }

        //std::cout << "applied inverse edge collapse op ... \n";
        //std::cout << "the momentum has been reverted to ... \n";
        std::vector<Node *> oneRing;
        Vec3 pBef = preCollapseMomentum(collapsedNode, oneRing);
        //std::cout << "pRevert: " << pBef << std::endl;

    }



// Uniform dv update idea to edge flip momentum conservation
/*void applyFlipEdgeOp(RemeshOp& op, Cloth& cloth){

	if( ::magic.conserve_momentum ){

		std::cout << "applying flip edge op \n";

		Edge* edgeToFlip = op.removed_edges[0];
		Edge* newEdge = op.added_edges[0];

		// 1. Identify the four verts
		Face* f0 = op.removed_faces[0];
		Face* f1 = op.removed_faces[1];
		std::vector<Vert*> verts;
		include(f0->v[0],verts);
		include(f0->v[1],verts);
		include(f0->v[2],verts);
		include(f1->v[0],verts);
		include(f1->v[1],verts);
		include(f1->v[2],verts);

		// 2. Compute momentum of the verts
		Vec3 pBef = Vec3(0);
		for(int i = 0; i < verts.size(); i++){
			pBef += verts[i]->node->m * verts[i]->node->v;
		}

		std::cout << "pBef: " << pBef << std::endl;

		// 3. Flip
		op.apply(cloth.mesh);
		compute_ms_data(cloth.mesh);
		compute_masses(cloth);

		// 4. Compute new momentum there
		Vec3 pAft = Vec3(0);
		for(int i = 0; i < verts.size(); i++){
			pAft += verts[i]->node->m * verts[i]->node->v;
		}

		std::cout << "pAft: " << pAft << std::endl;

		// 5. Compute dp
		Vec3 dp = pAft - pBef;

		std::cout << "dp: " << dp << std::endl;

		// 6. Compute dv
		Vec3 dv = dp / (verts[0]->node->m + verts[1]->node->m + verts[2]->node->m + verts[3]->node->m );

		std::cout << "dv: " << dv << std::endl;

		// 7. Buffer v
		verts[0]->node->vBuffer = verts[0]->node->v;
		verts[1]->node->vBuffer = verts[1]->node->v;
		verts[2]->node->vBuffer = verts[2]->node->v;
		verts[3]->node->vBuffer = verts[3]->node->v;

		// 8. Apply dv
		verts[0]->node->v -= dv;
		verts[1]->node->v -= dv;
		verts[2]->node->v -= dv;
		verts[3]->node->v -= dv;

		Vec3 pFix = Vec3(0);
		for(int i = 0; i < verts.size(); i++){
			pFix += verts[i]->node->m * verts[i]->node->v;
		}

		std::cout << "pFix: " << pFix << std::endl;
		std::cout << "net momentum change: " << pFix - pBef << std::endl;

	} else {
		op.apply(cloth.mesh);
	}

}*/




// METHOD OF PAPER, AREA WEIGHTING
/*
void applyFlipEdgeOp(RemeshOp& op, Cloth& cloth){

	if( ::magic.conserve_momentum ){

		std::cout << "applying flip edge op \n";

		// 0. Defining the initial and final edges
		Edge* ei = op.removed_edges[0];
		Edge* ef = op.added_edges[0];

		// 1. Defining the source and target faces
		Face* fs[2];
		Face* ft[2];
	    fs[0] = op.removed_faces[0];
		fs[1] = op.removed_faces[1];
		ft[0] = op.added_faces[0];
		ft[1] = op.added_faces[1];


		// 2. Computing momentum of each source face
		Vec3 ps[2];
		Vec3 pt[2];
		for(int f = 0; f < 2; f++){
			Vec3 velNodeA = fs[f]->v[0]->node->v;
			Vec3 velNodeB = fs[f]->v[1]->node->v;
			Vec3 velNodeC = fs[f]->v[2]->node->v;
			ps[f] = (fs[f]->m / 3.0) * (velNodeA + velNodeB + velNodeC);
		}

		// 3. Getting mat space positions of the inital edge verts
		Vec2 us0 = edge_vert(ei,0,0)->u;
		Vec2 us1 = edge_vert(ei,0,1)->u;

		std::cout << "flipping ... \n";

		// 4. Flip
		op.apply(cloth.mesh);
		compute_ms_data(cloth.mesh);
		compute_masses(cloth);

		std::cout << "flipped \n";

		// 5. Computing momentum of the target faces
		for(int f = 0; f < 2; f++){
			Vec3 velNodeA = ft[f]->v[0]->node->v;
			Vec3 velNodeB = ft[f]->v[1]->node->v;
			Vec3 velNodeC = ft[f]->v[2]->node->v;
			pt[f] = (ft[f]->m / 3.0) * (velNodeA + velNodeB + velNodeC);
		}

		std::cout << "computed momentum \n";

		// 6. Getting mat space positions of the final edge verts
		Vec2 ut0 = edge_vert(ef,0,0)->u;
		Vec2 ut1 = edge_vert(ef,0,1)->u;

		std::cout << "got mat space positions \n";

		std::cout << "ut0: " << us0 << std::endl;
		std::cout << "ut1: " << us1 << std::endl;
		std::cout << "ut2: " << ut0 << std::endl;
		std::cout << "ut3: " << ut1 << std::endl;

		// 7. Computing the bary weights of intersection
		Vec3 x0 = Vec3(ut0[0],ut0[1],1.0);
		Vec3 x1 = Vec3(ut1[0],ut1[1],1.0);
		Vec3 y0 = Vec3(us0[0],us0[1],0);
		Vec3 y1 = Vec3(us1[0],us1[1],0);
		Vec3 n = cross(normalize(x1-x0), normalize(y1-y0));
		n = normalize(n);
		double h = dot(x0-y0, n);
		double a0 = stp(y1-x1, y0-x1, n), a1 = stp(y0-x0, y1-x0, n),
			   b0 = stp(x0-y1, x1-y1, n), b1 = stp(x1-y0, x0-y0, n);
		double w[4];
		w[0] = a0/(a0 + a1);
		w[1] = a1/(a0 + a1);
		w[2] = -b0/(b0 + b1);
		w[3] = -b1/(b0 + b1);


		std::cout << "x0: " << x0 << std::endl;
		std::cout << "x1: " << x1 << std::endl;
		std::cout << "y0: " << y0 << std::endl;
		std::cout << "y1: " << y1 << std::endl;

		std::cout << "w0: " << w[0] << std::endl;
		std::cout << "w1: " << w[1] << std::endl;
		std::cout << "w2: " << w[2] << std::endl;
		std::cout << "w3: " << w[3] << std::endl;



		// 8. Computing the desired momentum for each face
		Vec3 ptPrime;
		if( is_inside(ut0,fs[0]) ){
			ptPrime = w[1]*ps[0] + w[0]*ps[1];
		} else if( is_inside(ut0,fs[1]) ){
			ptPrime = w[1]*ps[1] + w[0]*ps[0];
		} else {
			std::cout << "error not inside anything!\n";
			exit(0);
		}

		std::cout << "ps[0]: " << ps[0] << std::endl;
		std::cout << "ps[1]: " << ps[1] << std::endl;
		std::cout << "pt[0]: " << pt[0] << std::endl;
		std::cout << "pt[1]: " << pt[1] << std::endl;
		std::cout << "ptprime: " << ptPrime << std::endl;

		// 9. Buffering each node velocity before applying adjustment to conserve momentum
		std::vector<Vert*> verts;
		include(ft[0]->v[0],verts);
		include(ft[0]->v[1],verts);
		include(ft[0]->v[2],verts);
		include(ft[1]->v[0],verts);
		include(ft[1]->v[1],verts);
		include(ft[1]->v[2],verts);
		verts[0]->node->vBuffer = verts[0]->node->v;
		verts[1]->node->vBuffer = verts[1]->node->v;
		verts[2]->node->vBuffer = verts[2]->node->v;
		verts[3]->node->vBuffer = verts[3]->node->v;

		// 10. Applying a velocity adjustment to each vert
		for(int f = 0; f < 2; f++){
			Vec3 dp = ptPrime - pt[f];
			for(int i = 0; i < 3; i++){
				Node* node = ft[f]->v[i]->node;
				node->v += dp/(node->m * 3.0);
			}
		}

	} else {
		op.apply(cloth.mesh);
	}

}
* */




// MY IDEA
    void applyFlipEdgeOp(RemeshOp &op, Cloth &cloth) {

        if (arcsim::magic.conserve_momentum) {

            // 0. Defining the initial and final edges
            Edge *ei = op.removed_edges[0];
            Edge *ef = op.added_edges[0];

            // 1. Defining the source and target faces
            Face *fs[2];
            Face *ft[2];
            fs[0] = op.removed_faces[0];
            fs[1] = op.removed_faces[1];
            ft[0] = op.added_faces[0];
            ft[1] = op.added_faces[1];

            // 2. Getting the verts before flip
            std::vector<Vert *> verts;
            include(fs[0]->v[0], verts);
            include(fs[0]->v[1], verts);
            include(fs[0]->v[2], verts);
            include(fs[1]->v[0], verts);
            include(fs[1]->v[1], verts);
            include(fs[1]->v[2], verts);

            // 3. Computing momentum before flip
            Vec3 pBef = Vec3(0, 0, 0);
            for (int i = 0; i < verts.size(); i++) {
                pBef += verts[i]->node->v * verts[i]->node->m;
            }

            // 4. Getting mat space positions of the inital edge verts
            Vec2 us0 = edge_vert(ei, 0, 0)->u;
            Vec2 us1 = edge_vert(ei, 0, 1)->u;

            // 5. Apply the flip
            op.apply(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);

            // 6. Getting the verts after flip
            verts.clear();
            include(ft[0]->v[0], verts);
            include(ft[0]->v[1], verts);
            include(ft[0]->v[2], verts);
            include(ft[1]->v[0], verts);
            include(ft[1]->v[1], verts);
            include(ft[1]->v[2], verts);

            // 7. Computing momentum after flip
            Vec3 pAft = Vec3(0, 0, 0);
            for (int i = 0; i < verts.size(); i++) {
                pAft += verts[i]->node->v * verts[i]->node->m;
            }

            // 8. Computing momentum gain/loss
            Vec3 dp = pAft - pBef;

            // 9. Getting mat space positions of the final edge verts
            Vec2 ut0 = edge_vert(ef, 0, 0)->u;
            Vec2 ut1 = edge_vert(ef, 0, 1)->u;

            // 10. Computing the intersection point of the edges
            Vec3 x0 = Vec3(ut0[0], ut0[1], 1.0);
            Vec3 x1 = Vec3(ut1[0], ut1[1], 1.0);
            Vec3 y0 = Vec3(us0[0], us0[1], 0);
            Vec3 y1 = Vec3(us1[0], us1[1], 0);
            Vec3 n = cross(normalize(x1 - x0), normalize(y1 - y0));
            n = normalize(n);
            double h = dot(x0 - y0, n);
            double a0 = stp(y1 - x1, y0 - x1, n), a1 = stp(y0 - x0, y1 - x0, n),
                    b0 = stp(x0 - y1, x1 - y1, n), b1 = stp(x1 - y0, x0 - y0, n);
            double w[4];
            w[0] = a0 / (a0 + a1);
            w[1] = a1 / (a0 + a1);
            w[2] = -b0 / (b0 + b1);
            w[3] = -b1 / (b0 + b1);
            Vec2 intersectPt = w[0] * ut0 + w[1] * ut1;

            // 11. Compute barycentric coords and weighted mass
            Face *face = get_enclosing_face(cloth.mesh, intersectPt);
            Vec3 b = get_barycentric_coords(intersectPt, face);
            double baryMass = b[0] * face->v[0]->node->m + b[1] * face->v[1]->node->m + b[2] * face->v[2]->node->m;

            // 12. Compute node velocity corrections
            Vec3 dv0 = b[0] * dp / baryMass;
            Vec3 dv1 = b[1] * dp / baryMass;
            Vec3 dv2 = b[2] * dp / baryMass;

            // 13. Buffer velocities before updating
            for (int i = 0; i < verts.size(); i++) {
                verts[i]->node->vBuffer = verts[i]->node->v;
            }

            // 14. Update velocities
            face->v[0]->node->v -= dv0;
            face->v[1]->node->v -= dv1;
            face->v[2]->node->v -= dv2;

            // 15. Momentum after fix
            Vec3 pFix = Vec3(0, 0, 0);
            for (int i = 0; i < verts.size(); i++) {
                pFix += verts[i]->node->m * verts[i]->node->v;
            }


        } else {
            op.apply(cloth.mesh);
        }

    }


    void applyInverseFlipEdgeOp(RemeshOp &op, Cloth &cloth) {

        op.inverse().apply(cloth.mesh);

        if (arcsim::magic.conserve_momentum) {

            Face *f0 = op.added_faces[0];
            Face *f1 = op.added_faces[1];
            std::vector<Vert *> verts;
            include(f0->v[0], verts);
            include(f0->v[1], verts);
            include(f0->v[2], verts);
            include(f1->v[0], verts);
            include(f1->v[1], verts);
            include(f1->v[2], verts);

            verts[0]->node->v = verts[0]->node->vBuffer;
            verts[1]->node->v = verts[1]->node->vBuffer;
            verts[2]->node->v = verts[2]->node->vBuffer;
            verts[3]->node->v = verts[3]->node->vBuffer;

        }


    }


    void revertMomentum(std::vector<Node *> &oneRing) {
        for (int i = 0; i < oneRing.size(); i++) {
            oneRing[i]->v = oneRing[i]->vBuffer;
        }
    }


    vector<int> sort_edges_by_length(const Face *face);

    RemeshOp try_edge_collapse(Edge *edge, int which, Mesh &mesh);

    bool improve_some_face(vector<Face *> &active, Cloth &cloth, CollapsableCallback callback, bool check_collision,
                           const vector<AccelStruct *> &obs_accs) {

        for (int f = 0; f < active.size(); f++) {
            Face *face = active[f];
            for (int e = 0; e < 3; e++) {
                Edge *edge = face->adje[e];
                RemeshOp op;

                Node *nodeToCollapse = 0;

                if (op.empty() && callback(edge->n[0])) {
                    nodeToCollapse = edge->n[0];
                    op = try_edge_collapse(edge, 0, cloth.mesh);
                }
                if (op.empty() && callback(edge->n[1])) {
                    nodeToCollapse = edge->n[1];
                    op = try_edge_collapse(edge, 1, cloth.mesh);
                }

                if (op.empty()) continue;


                applyEdgeCollapseOp(op, cloth, nodeToCollapse);
                // compute_ms_data(cloth.mesh);
                // update_indices(cloth.mesh);

                if (false) {
                    // if(has_intersection(op.added_faces, &cloth.mesh, obs_accs)) {
                    op.inverse().apply(cloth.mesh);
                    op.inverse().done();
                    update_indices(cloth.mesh);
                } else {
                    op.done();
                    update_active(op, active);
                    vector<Face *> fix_active = op.added_faces;

                    RemeshOp flip_ops = argus_based_flip_edges(fix_active, cloth.mesh, obs_accs);
                    update_active(flip_ops, active);
                    flip_ops.done();
                    return true;
                }

            }
            remove(f--, active);
        }

        return false;
    }


    bool has_labeled_edges(const Node *node);

    bool can_collapse(const Edge *edge, int which);

    bool any_nearly_invalid(const vector<Edge *> edges) {
        for (int i = 0; i < edges.size(); i++)
            if (edge_metric(edges[i]) > 0.9) return true;
        return false;
    }

///////////////////////
/*
if( ::magic.conserve_momentum ){

		RemeshOp op;
		Node *node0 = edge->n[i], *node1 = edge->n[1-i];

		// Computing oneRing momentum before collapse
		Vec3 pBef = node0->m * node0->v;
		std::vector<Node*> oneRing;
		for(int e = 0; e < node0->adje.size(); e++){
			Edge *edge1 = node0->adje[e];
			Node *node2 = (edge1->n[0]!=node0) ? edge1->n[0] : edge1->n[1];
			oneRing.push_back(node2);
			pBef += node2->m * node2->v;
		}

		std::cout << "pBef: " << pBef << std::endl;

		op.removed_nodes.push_back(node0);
		for (int e = 0; e < node0->adje.size(); e++) {
			Edge *edge1 = node0->adje[e];
			op.removed_edges.push_back(edge1);
			Node *node2 = (edge1->n[0]!=node0) ? edge1->n[0] : edge1->n[1];
			if (node2 != node1 && !get_edge(node1, node2))
				op.added_edges.push_back(new Edge(node1, node2, edge1->theta_ideal,
												  edge1->label));
		}
		for (int s = 0; s < 2; s++) {
			Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1-i);
			if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i)))
				continue;
			op.removed_verts.push_back(vert0);
			for (int f = 0; f < vert0->adjf.size(); f++) {
				Face *face = vert0->adjf[f];
				op.removed_faces.push_back(face);
				if (!is_in(vert1, face->v)) {
					Vert *verts[3] = {face->v[0], face->v[1], face->v[2]};
					replace(vert0, vert1, verts);
					op.added_faces.push_back(new Face(verts[0], verts[1], verts[2],
													  face->label));
				}
			}
		}

		// local mass update
		for(int i = 0; i < oneRing.size(); i++){
			oneRing[i]->m = 0;
		}

		for(int i = 0; i < op.added_faces.size(); i++){
			Face* face = op.added_faces[i];
			compute_ms_data(face);
			face->a
			face->m =
			face->v[0]->node->m = (1.0/3.0)*face->m;

		Vec3 pAft = Vec3(0,0,0);
		for(int i = 0; i < oneRing.size(); i++){
			pAft += oneRing[i]->m * oneRing[i]->v;
		}



		std::cout << "pAft: " << pAft << std::endl;
		std::cout << "dp: " << pAft - pBef << std::endl;
		Vec3 dp = pAft - pBef;

		for(int i = 0; i < op.added_faces.size(); i++){
			Face* face = op.added_faces[i];
			if( is_inside(node0->verts[0]->u,face) ){
				Vec3 b = get_barycentric_coords(node0->verts[0]->u,face);
				face->v[0]->node->v -= b[0]*dp/face->v[0]->node->m;
				face->v[1]->node->v -= b[1]*dp/face->v[1]->node->m;
				face->v[2]->node->v -= b[2]*dp/face->v[2]->node->m;
			}
		}

	    Vec3 pFix = Vec3(0,0,0);
		for(int i = 0; i < oneRing.size(); i++){
			pFix += oneRing[i]->m * oneRing[i]->v;
		}

		std::cout << "dp after fix: " << pFix - pBef;




		return op;

		///////////////////////
*/

    RemeshOp try_edge_collapse(Edge *edge, int which, Mesh &mesh) {
        Node *node0 = edge->n[which], *node1 = edge->n[1 - which];

        if (node0->preserve
            || (is_seam_or_boundary(node0) && !is_seam_or_boundary(edge))
            || (has_labeled_edges(node0) && !edge->label) || node0->temp || node0->temp2 || node0->temp3) {
            return RemeshOp();
        }

        if (!can_collapse(edge, which)) {
            return RemeshOp();
        }

        return collapse_edge(edge, which);


    }

    bool has_labeled_edges(const Node *node) {
        for (int e = 0; e < node->adje.size(); e++)
            if (node->adje[e]->label)
                return true;
        return false;
    }

    bool can_collapse(const Edge *edge, int i) {
        for (int s = 0; s < 2; s++) {
            const Vert *vert0 = edge_vert(edge, s, i), *vert1 = edge_vert(edge, s, 1 - i);
            if (!vert0 || (s == 1 && vert0 == edge_vert(edge, 0, i)))
                continue;
            for (int f = 0; f < vert0->adjf.size(); f++) {
                const Face *face = vert0->adjf[f];
                if (is_in(vert1, face->v))
                    continue;
                const Vert *vs[3] = {face->v[0], face->v[1], face->v[2]};
                replace(vert0, vert1, vs);
                double a = wedge(vs[1]->u - vs[0]->u, vs[2]->u - vs[0]->u) / 2;
                double asp = aspect(vs[0]->u, vs[1]->u, vs[2]->u);
                if (a < 1e-6 || asp < remeshing->aspect_min)
                    return false;
                for (int e = 0; e < 3; e++)
                    if (vs[e] != vert1 && edge_metric(vs[NEXT(e)], vs[PREV(e)]) > 0.9)
                        return false;
            }
        }
        return true;
    }

    void print_op(RemeshOp &op) {
        cout << "added_verts:" << endl;
        for (int i = 0; i < op.added_verts.size(); i++) {
            cout << "\t" << op.added_verts[i] << endl;
        }
        cout << "removed_verts:" << endl;
        for (int i = 0; i < op.removed_verts.size(); i++) {
            cout << "\t" << op.removed_verts[i] << endl;
        }
        cout << "added_nodes:" << endl;
        for (int i = 0; i < op.added_nodes.size(); i++) {
            cout << "\t" << op.added_nodes[i] << endl;
        }
        cout << "removed_nodes:" << endl;
        for (int i = 0; i < op.removed_nodes.size(); i++) {
            cout << "\t" << op.removed_nodes[i] << endl;
        }
        cout << "added_faces:" << endl;
        for (int i = 0; i < op.added_faces.size(); i++) {
            cout << "\t" << op.added_faces[i] << endl;
        }
        cout << "removed_faces:" << endl;
        for (int i = 0; i < op.removed_faces.size(); i++) {
            cout << "\t" << op.removed_faces[i] << endl;
        }
        cout << "added_edges:" << endl;
        for (int i = 0; i < op.added_edges.size(); i++) {
            cout << "\t" << op.added_edges[i] << endl;
        }
        cout << "removed_edges:" << endl;
        for (int i = 0; i < op.removed_edges.size(); i++) {
            cout << "\t" << op.removed_edges[i] << endl;
        }
    }


    std::vector<Node *> get_vf_impact_nodes(Mesh &mesh) {

        std::vector<Node *> nodes;
        for (int i = 0; i < mesh.nodes.size(); i++) {
            mesh.nodes[i]->keep = false;
            if (mesh.nodes[i]->temp) {
                nodes.push_back(mesh.nodes[i]);
            }
        }
        for (int i = 0; i < mesh.edges.size(); i++) {
            mesh.edges[i]->collapsed = false;
        }
        return nodes;

    }


    std::vector<Node *> get_ee_impact_nodes(Mesh &mesh) {

        std::vector<Node *> nodes;
        for (int i = 0; i < mesh.nodes.size(); i++) {
            mesh.nodes[i]->keep = false;
            if (mesh.nodes[i]->temp2) {
                nodes.push_back(mesh.nodes[i]);
            }
        }
        for (int i = 0; i < mesh.edges.size(); i++) {
            mesh.edges[i]->collapsed = false;
        }
        return nodes;

    }

    RemeshOp argus_based_flip_edges(vector<Face *> &active, Mesh &mesh, const vector<AccelStruct *> &obs_accs) {

        bool flippedSomething;
        RemeshOp ops;
        unsigned int count = 0;
        int n_edges_prev = 0;
        do {

            //break;
            flippedSomething = false;
            std::vector<Edge *> indepEdges = independent_edges(find_edges_to_flip(active));
            if (indepEdges.size() == n_edges_prev) {
                break;
            }
            n_edges_prev = indepEdges.size();
            for (int k = 0; k < indepEdges.size(); k++) {

                RemeshOp flipEdgeOp = flip_edge(indepEdges[k]);
                if (flipEdgeOp.empty()) {
                    continue;
                }

                // compute_ws_data(mesh);
                Face *old_f1 = flipEdgeOp.removed_edges[0]->adjf[0];
                Face *old_f2 = flipEdgeOp.removed_edges[0]->adjf[1];
                double curv = fabs(dot(old_f1->n, old_f2->n));
                Edge *oldEdge = flipEdgeOp.removed_edges[0];
                Edge *newEdge = flipEdgeOp.added_edges[0];
                // double old_distance = distanceToObs(oldEdge, obs_accs);

                // applyFlipEdgeOp(flipEdgeOp,cloth); // with momentum conservation
                flipEdgeOp.apply(mesh);  // without momentum conservation

                // update_indices(mesh);
                // compute_ms_data(mesh);

                // double new_distance = distanceToObs(newEdge, obs_accs);

                if (any_inverted(flipEdgeOp.added_faces)) {
                    // if ( !has_proximity(flipEdgeOp.added_faces, obs_accs) ) {
                    // if ( has_intersection(flipEdgeOp.added_faces, &mesh, obs_accs) || any_inverted(flipEdgeOp.added_faces)) {
                    flipEdgeOp.inverse().apply(mesh);  // without momentum conservation
                    // applyInverseFlipEdgeOp(flipEdgeOp,cloth);  // with momentum conservation

                    flipEdgeOp.inverse().done();
                    // update_indices(mesh);
                    // RemeshOp op = split_edge(oldEdge);
                    // // if (oldEdge->collapsed) {
                    // //     op.added_nodes[0]->keep = true;
                    // // }
                    // op.apply(mesh);
                    // update_active(op, active);
                    // op.done();

                    // Node *node0 = oldEdge->n[0], *node1 = oldEdge->n[1];
                    // for (int v = 0; v < op.added_verts.size(); v++) {
                    //     Vert *vertnew = op.added_verts[v];
                    //     Vert *v0 = adjacent_vert(node0, vertnew),
                    //          *v1 = adjacent_vert(node1, vertnew);
                    //     vertnew->sizing = new Sizing(mean_vert_sizing(v0, v1));
                    // }
                    // op.added_nodes[0]->temp2 = true;
                    // std::cout << "didn't allow a flip. \n";
                    //exit(0);
                } else {
                    // if (new_distance - old_distance >= 0 || new_distance - .5*0.25*::magic.projection_thickness >= 0) {
                    // if (true) {
                    // save_mesh(mesh);
                    update_active(flipEdgeOp, active);
                    ops = compose(ops, flipEdgeOp);
                    // flipEdgeOp.done();
                    flippedSomething = true;
                    // std::cout << "allowed a flip.\n";
                }
                // flippedSomething = true;
            }

            count++;
        } while (active.size() > 0 && flippedSomething && count < mesh.verts.size());
        return ops;
    }

    bool collapse_some_verts(vector<Node *> &active, Cloth &cloth, vector<AccelStruct *> &obs_accs,
                             CollapsableCallback callback) {

        for (int a = 0; a < active.size(); a++) {
            Node *vfNode = active[a];
            if (arcsim::magic.consider_contact_force && norm(vfNode->r) > 1e-10) {
                continue;
            }

            if (vfNode->freeAge < arcsim::magic.node_lifetime) {
                continue;
            }

            for (int j = 0; j < vfNode->adje.size(); j++) {

                Edge *edge = vfNode->adje[j];

                RemeshOp edgeCollOp;


                Node *nodeToCollapse = 0;
                if (edge->n[0] == vfNode && callback(edge->n[0])) {
                    edgeCollOp = try_edge_collapse(edge, 0, cloth.mesh);
                    nodeToCollapse = edge->n[0];
                } else if (edge->n[1] == vfNode && callback(edge->n[1])) {
                    edgeCollOp = try_edge_collapse(edge, 1, cloth.mesh);
                    nodeToCollapse = edge->n[1];
                }
                Vec3 pos1 = edge->n[0]->x;
                Vec3 pos2 = edge->n[1]->x;

                if (edgeCollOp.empty()) {
                    continue;
                }
                //edgeCollOp.apply(mesh);
                applyEdgeCollapseOp(edgeCollOp, cloth, nodeToCollapse);

                compute_ms_data(cloth.mesh);

                // for (int i = 0; i < edgeCollOp.added_faces.size(); i++) {
                //     Face* af = edgeCollOp.added_faces[i];
                // }
                // If edge collapse failed, or the operation would result in an intersection, do nothing

                if (false) {
                    // if(has_intersection(edgeCollOp.added_faces, obs_accs) && !arcsim::magic.consider_contact_force){
                    //edgeCollOp.inverse().apply(mesh);
                    applyInverseEdgeCollapseOp(edgeCollOp, cloth);
                    edgeCollOp.inverse().done();
                    continue;

                } else {

                    // edgeCollOp.apply(mesh);
                    // cout << "added_nodes:" << edgeCollOp.added_nodes.size() << endl;
                    // cout << "removed_nodes:" << edgeCollOp.removed_nodes.size() << endl;
                    // cout << "before:" << active.size() << endl;
                    update_active(edgeCollOp, active);
                    // cout << "after:" << active.size() << endl;
                    // cout << endl;
                    edgeCollOp.done();
                    for (int i = 0; i < edgeCollOp.added_edges.size(); i++) {
                        edgeCollOp.added_edges[i]->collapsed = true;
                    }

                    set_identity_sizing(cloth.mesh);

                    // std::vector<Face*> active = edgeCollOp.added_faces;
                    vector<Face *> activeFaces = cloth.mesh.faces;

                    argus_based_flip_edges(activeFaces, cloth.mesh, obs_accs);

                    // for (int i = 0; i < mesh.nodes.size(); i++) {
                    //     if (mesh.nodes[i]->keep) {
                    //         exclude(mesh.nodes[i], active);
                    //     }
                    // }

                    return true;
                }
            }
            remove(a--, active);
        }
        return false;
    }


// Collapsing Faces (for impact handling)
    void collapse_cloth_impact_verts(Cloth &cloth, vector<AccelStruct *> &obs_accs) {

        std::vector<Node *> vfImpactNodes = get_vf_impact_nodes(cloth.mesh);
        while (collapse_some_verts(vfImpactNodes, cloth, obs_accs, VF_callback));
    }


// Collapsing Edges (for impact handling)
    void collapse_cloth_impact_edges(Cloth &cloth, vector<AccelStruct *> &obs_accs) {

        std::vector<Node *> eeImpactNodes = get_ee_impact_nodes(cloth.mesh);
        while (collapse_some_verts(eeImpactNodes, cloth, obs_accs, EE_callback));

    }

    void insert_nodes(vector<Cluster *> &clusters, Cloth &cloth, const vector<Mesh *> obs_meshes) {
        vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, false);
        VisualDebugger *vd = VisualDebugger::getInstance();
        for (int c = 0; c < clusters.size(); c++) {
            Cluster *cluster = clusters[c];
            if (cluster->node) {
                continue;
            }
            RemeshOp op;
            Face *encFace = get_enclosing_face(cloth.mesh, cluster->u);
            Vec3 b = get_barycentric_coords(cluster->u, encFace);
            // assert(norm(cluster->u - encFace->v[0]->u) > ::magic.merge_radius);
            // assert(norm(cluster->u - encFace->v[1]->u) > ::magic.merge_radius);
            // assert(norm(cluster->u - encFace->v[2]->u) > ::magic.merge_radius);
            int index = b[0] > b[1] ? (b[1] > b[2] ? 2 : 1) : (b[0] > b[2] ? 2 : 0);
            int i = (index + 1) % 3;
            int j = (index + 2) % 3;
            Vert *v1 = encFace->v[i];
            Vert *v2 = encFace->v[j];
            Vec2 v1_p = cluster->u - v1->u;
            Vec2 v1_v2 = v2->u - v1->u;
            double proj_len = dot(v1_p, v1_v2) / norm(v1_v2);
            Vec2 proj = proj_len * normalize(v1_v2);
            double d = norm(v1_p - proj);
            if (d <= arcsim::magic.merge_radius) {    // View as edge splitting
                Edge *edge = getCommonEdge(v1->node, v2->node);
                assert(edge != NULL);
                if (proj_len <= arcsim::magic.merge_radius) {
                    cluster->node = v1->node;
                    continue;
                }
                if (norm(v1_v2) - proj_len <= arcsim::magic.merge_radius) {
                    cluster->node = v2->node;
                    continue;
                }
                double bi = 1 - proj_len / norm(v1_v2);
                double w = edge->n[0] == v1->node ? bi : 1 - bi;
                op = split_edge(edge, w);
                op.added_nodes[0]->temp2 = true;
            } else {
                op = split_face(encFace, b[0], b[1], b[2]);
                op.added_nodes[0]->temp = true;
            }
            assert(!op.empty());
            op.apply(cloth.mesh);
            op.done();
            // Updating sizing data
            set_identity_sizing(cloth.mesh);
            // Updating indices and other data
            update_indices(cloth.mesh);
            compute_ms_data(cloth.mesh);
            compute_masses(cloth);
            assert(op.added_nodes.size() == 1);
            cluster->node = op.added_nodes[0];

            argus_based_flip_edges(cloth.mesh.faces, cloth.mesh, obs_accs);
        }
    }

    void insert_nodes(ImpactPoint *point, Cloth &cloth) {
        RemeshOp op;
        if (point->verts.size() == 2) { // edge impacts
            const Vert *v1 = point->verts[0];
            const Vert *v2 = point->verts[1];
            Edge *cEdge = getCommonEdge(v1->node, v2->node);
            if (cEdge) {
                Vec2 proj;
                material_ve_projection(point->u, v1->u, v2->u, proj);
                Vec2 end = cEdge->n[1] == v1->node ? v1->u : v2->u;
                double w = norm(proj - end) / norm(v1->u - v2->u);
                op = split_edge(cEdge, w);
                op.added_nodes[0]->temp2 = true;
            }
        }
        if (point->verts.size() == 3 || op.empty()) {
            Face *encFace = get_enclosing_face(cloth.mesh, point->u);
            Vec3 b = get_barycentric_coords(point->u, encFace);
            op = split_face(encFace, b[0], b[1], b[2]);
            op.added_nodes[0]->temp = true;
            op.added_nodes[0]->preserve = true;
        }
        assert(!op.empty());
        op.apply(cloth.mesh);
        op.done();
        // Updating sizing data
        set_identity_sizing(cloth.mesh);
        // Updating indices and other data
        update_indices(cloth.mesh);
        compute_ms_data(cloth.mesh);
        compute_masses(cloth);
        point->node = op.added_nodes[0];
        // argus_based_flip_edges(op.added_faces, cloth, obs_accs);
    }

    vector<ArgusImpact> convert_to_argus(const vector<NodalImpact *> nImpacts) {
        vector<ArgusImpact> aImpacts;
        for (int n = 0; n < nImpacts.size(); n++) {
            ArgusImpact aImpact;
            NodalImpact *nImpact = nImpacts[n];


            aImpact.nodeA = nImpact->p1->node;
            aImpact.nodeB = nImpact->p2->type == ImpactPoint::CLOTH ? nImpact->p2->node : NULL;
            aImpact.posB = nImpact->p2->x;  // posB and velB won't be used if nodeB != NULL
            aImpact.nodeA->freeAge = 0;
            if (aImpact.nodeB != NULL) {
                aImpact.nodeB->freeAge = 0;
            }

            aImpact.velB = nImpact->p2->v;
            aImpact.normal = nImpact->normal;
            aImpacts.push_back(aImpact);
        }
        return aImpacts;
    }

}