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

#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "remesh.hpp"
#include "cloth.hpp"
#include "constraint.hpp"
#include "collisionutil.hpp"
#include "simulation.hpp"
#include <vector>
using namespace std;


namespace arcsim {
//	extern Simulation sim;

// This impact struct is used to generate initial impact data. In the Argus project we analyze this
// and produce a corresponding ArgusImpact 
	struct Impact {
		bool inverted;
		bool self;
		bool debug;
		enum Type {
			VF, EE, VE, VV
		} type;
		double t;
		Node *nodes[4];
		const Vert *verts[4];
		double w[4];
		double dvN; // for proximity impacts, computed change in velocity in normal direction from repulsion impulse
		Vec3 n;
		Vec3 obsX;
		Vec3 obsV;
		vector<Impact> adjacentImpacts;
		Vec2 matPosA; // material space coordinates of contact location on first primitive
		Vec2 matPosB; // material space coordinates of contact location on second primitive
		// Vec2 u[4];
		// Vec3 obsNormal;  // normal of the obstacle  (this is just here for debugging purposes)
		Impact() {
			inverted = false;
			self = false;
			debug = false;
		}

		Impact(Type type, const Node *n0, const Node *n1, const Node *n2,
			   const Node *n3) : type(type) {
			nodes[0] = (Node *) n0;
			nodes[1] = (Node *) n1;
			nodes[2] = (Node *) n2;
			nodes[3] = (Node *) n3;
			inverted = false;
			self = false;
			debug = false;
			// obsNormal = Vec3(0);

		}
	};

// ArgusImpact data structure 
	struct ArgusImpact {

		enum Type {
			VF, EE, VE, VV
		} type;
		Node *nodeA;  // node at colliding point A
		Node *nodeB;  // node at colliding point B (if there's a node. else pointer defaults to 0)
		Vec3 posB;    // world-space position of colliding point B
		Vec3 velB;    // velocity of colliding point B
		Vec3 normal;  // impact normal
		// Vec3 obsNormal; // normal of the obstacle surface (this is just here for debugging purposes, and assumes no self-contact)
		bool inverted;
		bool self;
		bool debug;

		ArgusImpact() {
			nodeA = 0;
			nodeB = 0;
			posB = Vec3(0);
			velB = Vec3(0);
			normal = Vec3(0);
			// obsNormal = Vec3(0);
			inverted = false;
			self = false;
			debug = false;
		}
	};


	vector<Impact> find_continuous_impacts(const vector<AccelStruct *> &acc,
										   const vector<AccelStruct *> &obs_accs);

	vector<Impact> find_proximity_impacts(const vector<AccelStruct *> &acc,
										  const vector<AccelStruct *> &obs_accs);

	vector<Impact> find_proximity_impacts(const vector<Face *> &added_faces,
										  const vector<AccelStruct *> &obs_accs);

	vector<Impact> prepare_impacts(std::vector<Mesh *> &meshes,
								   const std::vector<Mesh *> &obs_meshes);

	bool has_proximity(const std::vector<Face *> &added_faces, const std::vector<AccelStruct *> &obs_accs);

	double distanceToObs(Edge *edge, const std::vector<AccelStruct *> &obs_accs, double *w = NULL,
						 Vec3 *n = NULL, Vec3 *obsX = NULL);

/*
void collision_response_bh (std::vector<Mesh*> &meshes,
                         const std::vector<Constraint*> &cons,
                         const std::vector<Mesh*> &obs_meshes,
                         double dt);
*/
}
#endif
