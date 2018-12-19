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

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "mesh.hpp"
#include "util.hpp"

namespace arcsim {
    double signed_vf_distance(const Vec3 &x,
                              const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
                              Vec3 *n, double *w);

    double signed_ee_distance(const Vec3 &x0, const Vec3 &x1,
                              const Vec3 &y0, const Vec3 &y1,
                              Vec3 *n, double *w);

    double unsigned_vf_distance(const Vec3 &x,
                                const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
                                Vec3 *n, double *w);

    double unsigned_ee_distance(const Vec3 &x0, const Vec3 &x1,
                                const Vec3 &y0, const Vec3 &y1,
                                Vec3 *n, double *w);

    bool set_unsigned_ve_distance(const Vec3 &x, const Vec3 &y0, const Vec3 &y1,
                                  double *_d, Vec3 *_n,
                                  double *_wx, double *_wy0, double *_wy1);


    Vec3 get_barycentric_coords(const Vec2 &point, const Face *face);

    bool is_inside(const Vec2 &point, const Face *f);

    Face *get_enclosing_face(const Mesh &mesh, const Vec2 &u,
                             Face *starting_face_hint = NULL);

    Edge *get_nearest_edge(const Mesh &mesh, const Vec2 &);

// Get the projection of a node to an edge, return the distance.
    double material_ve_projection(const Vec2 &x, const Vec2 &y0, const Vec2 &y1, Vec2 &proj);

    enum Space {
        PS, WS
    }; // plastic space, world space

    template<Space s>
    const Vec3 &pos(const Node *node);

    template<Space s>
    Vec3 &pos(Node *node);

    template<Space s>
    Vec3 nor(const Face *face);

    template<Space s>
    double dihedral_angle(const Edge *edge);

    template<Space s>
    Mat2x2 curvature(const Face *face);

    inline double area(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2) {
        return 0.5 * wedge(u1 - u0, u2 - u0);
    }

    inline double area(const Face *face) {
        return area(face->v[0]->u, face->v[1]->u, face->v[2]->u);
    }

    inline double perimeter(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2) {
        return norm(u0 - u1) + norm(u1 - u2) + norm(u2 - u0);
    }

    inline double aspect(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2) {
        double a = area(u0, u1, u2);
        double p = perimeter(u0, u1, u2);
        return 12 * sqrt(3) * a / sq(p);
    }

    inline double aspect(const Face *face) {
        return aspect(face->v[0]->u, face->v[1]->u, face->v[2]->u);
    }

    double unwrap_angle(double theta, double theta_ref);

    Face *getCommonFace(const Vert *v1, const Vert *v2, const Vert *v3);

    Edge *getCommonEdge(const Node *nodeA, const Node *nodeB);

    Edge *getCommonEdge(const Face *face0, const Face *face1);

    double getSharpAngle(Edge *edge);

    double getSharpAngle(Node *node);
}
#endif
