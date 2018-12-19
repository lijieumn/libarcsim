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

#include "collision.hpp"

#include "geometry.hpp"
#include "magic.hpp"
#include "optimization.hpp"
#include "simulation.hpp"
#include "timer.hpp"
#include "proximity.hpp"
#include <algorithm>
#include <fstream>
#include <omp.h>
using namespace std;


namespace arcsim {

	static const int max_iter = 100;

	static double thickness = arcsim::magic.detection_thickness;
	static double real_thickness = arcsim::magic.projection_thickness;

//static const double &thickness = ::magic.repulsion_thickness; // use this for bridson harmon
	static double ratio = .5;
	static double &dt = magic.dt;

	static double obs_mass;
	static bool deform_obstacles;

	static vector<Vec3> xold;
	static vector<Vec3> xold_obs;

	vector<Impact> find_continuous_impacts(const vector<AccelStruct *> &acc,
										   const vector<AccelStruct *> &obs_accs);

	vector<Impact> find_proximity_impacts(const vector<AccelStruct *> &acc,
										  const vector<AccelStruct *> &obs_accs);

	vector<Impact> find_proximity_impacts(const vector<Face *> &added_faces,
										  const vector<AccelStruct *> &obs_accs);

	vector<Impact> find_impacts(const vector<AccelStruct *> &acc,
								const vector<AccelStruct *> &obs_accs);

	vector<Impact> independent_impacts(const vector<Impact> &impacts);


	ostream &operator<<(ostream &out, const Impact &imp);

	static double *w = NULL;
	static Vec3 n;
	static Vec3 obsX;
	static double closest = 0xffff;

	void for_overlapping_edges(Edge *edge, BVHNode *node, double thickness);

	double distanceToObs(Edge *edge, const std::vector<AccelStruct *> &obs_accs,
						 double *w, Vec3 *n, Vec3 *obsX) {

		if (!edge) {
			return 0xffff;
		}
		closest = 0xffff;
		for (int o = 0; o < obs_accs.size(); o++) {
			if (obs_accs[o]->root) {
				for_overlapping_edges(edge, obs_accs[o]->root, arcsim::thickness);
			}
		}
		if (w && n && obsX) {
			w[0] = arcsim::w[0];
			w[1] = arcsim::w[1];
			*n = arcsim::n;
			*obsX = arcsim::obsX;
		}

		return closest;
	}

	void for_overlapping_edges(Edge *edge, BVHNode *node, double thickness) {
		if (!node->_active) {
			return;
		}
		BOX box;
		box += edge->n[0]->x;
		box += edge->n[0]->x;
		if (!overlap(box, node->_box, thickness)) {
			return;
		}
		if (node->isLeaf()) {
			// callback(face, node->getFace());
			for (int e = 0; e < 3; e++) {
				Edge *obsEdge = node->getFace()->adje[e];
				Node *n1 = edge->n[0];
				Node *n2 = edge->n[1];
				Node *n3 = obsEdge->n[0];
				Node *n4 = obsEdge->n[1];

				Vec3 n;
				double *w = new double[4];
				Vec3 x0 = edge->n[0]->x;
				Vec3 x1 = edge->n[1]->x;
				Vec3 x2 = obsEdge->n[0]->x;
				Vec3 x3 = obsEdge->n[1]->x;
				double d = signed_ee_distance(x0, x1, x2, x3, &n, w);
				d = fabs(d);
				bool inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);
				if (inside && (d < closest || closest < 0)) {
					closest = d;
                    arcsim::obsX = -x2 * w[2] - x3 * w[3];
                    arcsim::n = normalize(w[0] * x0 + w[1] * x1 + w[2] * x2 + w[3] * x3);
                    arcsim::w = w;
				}
			}
		} else {
			for_overlapping_edges(edge, node->getLeftChild(), thickness);
			for_overlapping_edges(edge, node->getRightChild(), thickness);
		}
	}

	bool find_node_in_mesh(const Node *node, const Mesh &mesh) {
		for (int n = 0; n < mesh.nodes.size(); n++) {
			if (node == mesh.nodes[n]) {
				return true;
			}
		}
		return false;
	}


	bool find_node_in_meshes(const Node *node, const std::vector<Mesh *> &meshes) {
		for (int m = 0; m < meshes.size(); m++) {
			if (find_node_in_mesh(node, *meshes[m]))
                return true;
		}
        return false;
	}

	bool has_proximity(const std::vector<Face *> &added_faces, const std::vector<AccelStruct *> &obs_accs) {
		vector<Impact> impacts = find_proximity_impacts(added_faces, obs_accs);
		return !impacts.empty();
	}

	vector<Impact> prepare_impacts(std::vector<Mesh *> &meshes, const std::vector<Mesh *> &obs_meshes) {

		vector<AccelStruct *> accs = create_accel_structs(meshes, true);
		vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, true);

		// Find impacts from collisions between cloths and obstacles
		vector<Impact> impacts = find_impacts(accs, obs_accs);
		// vector<Impact> impacts = find_proximity_impacts(accs, obs_accs);

		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);

		// Looping through all the raw impacts
		for (int i = int(impacts.size()) - 1; i >= 0; i--) {

			Impact &impact = impacts[i];
			// bool first = find_node_in_meshes(impact.nodes[0],meshes);
			// bool second = find_node_in_meshes(impact.nodes[1],meshes);
			// bool third = find_node_in_meshes(impact.nodes[2],meshes);
			// bool fourth = find_node_in_meshes(impact.nodes[3],meshes);
			bool first = impact.nodes[0] ? impact.nodes[0]->inMesh : false;
			bool second = impact.nodes[1] ? impact.nodes[1]->inMesh : false;
			bool third = impact.nodes[2] ? impact.nodes[2]->inMesh : false;
			bool fourth = impact.nodes[3] ? impact.nodes[3]->inMesh : false;
			bool valid = false;

			// Switch by impact type (VF, EE, VE, VV)
			switch (impact.type) {

				// Vertex-Face Impact
				case Impact::VF: {

					bool standardVF = first && !second && !third && !fourth && arcsim::magic.vfImpacts;
					bool invertedVF = !first && second && third && fourth && arcsim::magic.fvImpacts;
					bool selfVF = first && second && third && fourth && arcsim::magic.vfImpacts && arcsim::magic.self_contact;
					valid = standardVF || invertedVF || selfVF;
					impact.inverted = invertedVF;
					impact.self = selfVF;
					break;
				}

					// Edge-Edge Impact
				case Impact::EE: {

					bool standardEE = first && second && !third && !fourth && arcsim::magic.eeImpacts;
					bool selfEE = first && second && third && fourth && arcsim::magic.eeImpacts && arcsim::magic.self_contact;
					valid = standardEE || selfEE;
					impact.self = selfEE;
					break;
				}

					// Vertex-Edge Impact
				case Impact::VE: {

					bool standardVE = first && !second && !third && arcsim::magic.veImpacts;
					bool invertedVE = !first && second && third && arcsim::magic.evImpacts;
					bool selfVE = first && second && third && arcsim::magic.veImpacts && arcsim::magic.self_contact;
					valid = standardVE || invertedVE || selfVE;
					impact.inverted = invertedVE;
					impact.self = selfVE;
					break;
				}

					// Vertex-Vertex Impact
				case Impact::VV: {

					bool standardVV = first && !second && arcsim::magic.vvImpacts;
					bool selfVV = first && second && arcsim::magic.vvImpacts && arcsim::magic.self_contact;
					valid = standardVV || selfVV;
					impact.self = selfVV;
					break;
				}

			}

			// If the impact didn't have a valid or enabled signature, erase it
			if (!valid) {
				impacts.erase(impacts.begin() + i);
			}

		}


		return impacts;

	}










// Impacts

	static int nthreads = 0;
	static vector<Impact> *impacts = NULL;

	void find_face_continuous_impacts(const Face *face0, const Face *face1);

	void find_face_proximity_impacts(const Face *face0, const Face *face1);

	void find_face_impacts(const Face *face0, const Face *face1);

	vector<Impact> find_continuous_impacts(const vector<AccelStruct *> &accs,
										   const vector<AccelStruct *> &obs_accs) {
		if (!impacts) {
			arcsim::nthreads = omp_get_max_threads();
			arcsim::impacts = new vector<Impact>[arcsim::nthreads];
		}
		for (int t = 0; t < arcsim::nthreads; t++)
			arcsim::impacts[t].clear();
		for_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_continuous_impacts);
		vector<Impact> impacts;
		for (int t = 0; t < arcsim::nthreads; t++)
			append(impacts, arcsim::impacts[t]);

		return impacts;
	}

	static int frame = 0;

	vector<Impact> find_impacts(const vector<AccelStruct *> &accs, const vector<AccelStruct *> &obs_accs) {
		if (!impacts) {
			arcsim::nthreads = omp_get_max_threads();
			arcsim::impacts = new vector<Impact>[arcsim::nthreads];
		}
		for (int t = 0; t < arcsim::nthreads; t++)
			arcsim::impacts[t].clear();
		arcsim::thickness = (frame < arcsim::magic.init_frames) ? 4 * arcsim::magic.detection_thickness : arcsim::magic.detection_thickness;
		if (arcsim::magic.use_stack_overloop_face_check) {
			find_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_impacts);
		} else {
			for_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_impacts);
		}
		//for_overlapping_faces(accs, accs, arcsim::thickness, find_face_proximity_impacts);
		vector<Impact> impacts;
		for (int t = 0; t < arcsim::nthreads; t++)
			append(impacts, arcsim::impacts[t]);
		arcsim::frame++;

		return impacts;

	}

	vector<Impact> find_proximity_impacts(const vector<AccelStruct *> &accs,
										  const vector<AccelStruct *> &obs_accs) {
		if (!impacts) {
			arcsim::nthreads = omp_get_max_threads();
			arcsim::impacts = new vector<Impact>[arcsim::nthreads];
		}
		for (int t = 0; t < arcsim::nthreads; t++)
			arcsim::impacts[t].clear();
		if (arcsim::magic.use_stack_overloop_face_check) {
			find_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_proximity_impacts);
		} else {
			for_overlapping_faces(accs, obs_accs, arcsim::thickness, find_face_proximity_impacts);
		}
		//for_overlapping_faces(accs, accs, arcsim::thickness, find_face_proximity_impacts);
		vector<Impact> impacts;
		for (int t = 0; t < arcsim::nthreads; t++)
			append(impacts, arcsim::impacts[t]);

		return impacts;
	}

	vector<Impact> find_proximity_impacts(const vector<Face *> &added_faces,
										  const vector<AccelStruct *> &obs_accs) {
		if (!impacts) {
			arcsim::nthreads = omp_get_max_threads();
			arcsim::impacts = new vector<Impact>[arcsim::nthreads];
		}
		for (int t = 0; t < arcsim::nthreads; t++)
			arcsim::impacts[t].clear();
		for_overlapping_faces(added_faces, obs_accs, arcsim::thickness, find_face_proximity_impacts);
		vector<Impact> impacts;
		for (int t = 0; t < arcsim::nthreads; t++)
			append(impacts, arcsim::impacts[t]);

		return impacts;

	}

	Edge *getEdge(Node *nodeA, Node *nodeB) {
		for (int i = 0; i < nodeA->adje.size(); i++) {
			for (int j = 0; j < nodeB->adje.size(); j++) {
				if (nodeA->adje[i] == nodeB->adje[j]) {
					return nodeA->adje[i];
				}
			}
		}

		return NULL;

	}


	Vec3 projectToPlane(Vec3 point, Face *face, bool forward = false) {

		Vec3 faceCenter;
		if (forward) {
			faceCenter = (1.0 / 3.0) * (face->v[0]->node->x + face->v[0]->node->v * dt
										+ face->v[1]->node->x + face->v[1]->node->v * dt
										+ face->v[2]->node->x + face->v[1]->node->v * dt);
		} else {
			faceCenter = (1.0 / 3.0) * (face->v[0]->node->x
										+ face->v[1]->node->x
										+ face->v[2]->node->x);
		}

		Vec3 faceNormal = face->n;
		Vec3 disp = point - faceCenter;
		double dist = dot(disp, faceNormal);

		return point - dist * faceNormal;

	}

	bool inWedge(Vec3 pos, Edge *edge, bool forward = false) {

		Vec3 x = pos;
		bool in = true;
		for (int s = 0; s < 2; s++) {
			const Face *face = edge->adjf[s];
			if (!face)
				continue;

			Vec3 pointInPlane = projectToPlane(pos, (Face *) face, forward);
			Vec3 x0, x1, x2;
			if (forward) {
				x0 = face->v[0]->node->x + face->v[0]->node->v * dt;
				x1 = face->v[1]->node->x + face->v[1]->node->v * dt;
				x2 = face->v[2]->node->x + face->v[1]->node->v * dt;
			} else {
				x0 = face->v[0]->node->x;
				x1 = face->v[1]->node->x;
				x2 = face->v[2]->node->x;
			}
			Vec3 n;
			double w[4];
			double *wptr = w;
			signed_vf_distance(pointInPlane, x0, x1, x2, &n, wptr);

			int i = 0;
			for (; i < 3; i++) {
				if (face->v[i]->node != edge->n[0] && face->v[i]->node != edge->n[1])
					break;
			}

			double baryOppositeSide = -wptr[i + 1];

			in &= (baryOppositeSide < 0);

		}

		return in;

	}


	bool inCone(Vec3 pos, Node *node, bool forward = false) {

		bool in = true;

		for (int i = 0; i < node->adje.size(); i++) {

			Edge *edge = node->adje[i];
			Node *otherNode = edge->n[0] == node ? edge->n[1] : edge->n[0];

			Vec3 a, b;
			if (forward) {
				a = otherNode->x + otherNode->v * dt;
				b = node->x + node->v * dt;
			} else {
				a = otherNode->x;
				b = node->x;
			}
			Vec3 e = b - a;
			Vec3 r = pos - a;

			Vec3 projPt = a + dot(r, e) / dot(e, e) * e;

			double distPA = norm(a - projPt);
			double distPB = norm(b - projPt);

			Vec3 unitVecPA = normalize(a - projPt);
			Vec3 unitVecPB = normalize(b - projPt);

			bool bCloser = distPB < distPA;
			bool sameDir = fabs(dot(unitVecPA, unitVecPB) - 1.0) < 0.001;

			in &= (bCloser && sameDir);

		}

		return in;

	}


	Vec2 getMatPos(const Edge *edge, double b) {

		double b0 = b;
		double b1 = 1.0 - b;

		Vert *vtest0 = edge_vert(edge, 0, 0);
		Vert *vtest1 = edge_vert(edge, 0, 1);
		Vert *vtest2 = edge_vert(edge, 1, 0);
		Vert *vtest3 = edge_vert(edge, 1, 1);

		for (int s = 0; s < 2; s++) {

			if (!edge->adjf[s]) {
				continue;
			}

			Vert *v0 = edge_vert(edge, s, 0);
			Vert *v1 = edge_vert(edge, s, 1);

			return b0 * v0->u + b1 * v1->u;

		}

        assert(false);
        return Vec2(0,0);

	}

	bool continuous_collision_test(Impact::Type type, Node *nodes[4], Impact &impact);

	bool vf_continuous_collision_test_new(const Vert *vert, const Face *face, Impact &impact);

	bool ee_continuous_collision_test_new(const Edge *edge0, const Edge *edge1, Impact &impact);

	Vert *getVert(const Edge *edge, int i);

// Continuous Collision tests
	bool vf_continuous_collision_test(const Vert *vert, const Face *face, Impact &impact);

	bool ee_continuous_collision_test(const Edge *edge0, const Edge *edge1, Impact &impact);

// Proximity Collision tests
	bool vf_proximity_collision_test(const Vert *vert, const Face *face, Impact &impact, bool forward = false);

	bool ee_proximity_collision_test(const Edge *edge0, const Edge *edge1, Impact &impact, bool forward = false);

	bool ve_proximity_collision_test(const Vert *vert, const Edge *edge, Impact &impact, bool forward = false);

	bool vv_proximity_collision_test(const Vert *vert0, const Vert *vert1, Impact &impact, bool forward = false);


// Identifying impacts between two faces
	void find_face_proximity_impacts(const Face *face0, const Face *face1) {

		int t = omp_get_thread_num();
		Impact impact;

		// FV test : face1 vert with face0
		if (arcsim::magic.fvImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int n = 0; n < face1->rNode.size(); n++) {
					Node *node = face1->rNode[n];
					if (vf_proximity_collision_test(node->verts[0], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					}
				}
			} else {
				for (int v = 0; v < 3; v++)
					if (vf_proximity_collision_test(face1->v[v], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					}
			}
		}

		// VF test : face0 vert with face1
		if (arcsim::magic.vfImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int n = 0; n < face0->rNode.size(); n++) {
					Node *node = face0->rNode[n];
					if (vf_proximity_collision_test(node->verts[0], face1, impact)) {
						arcsim::impacts[t].push_back(impact);
					}
				}
			} else {
				for (int v = 0; v < 3; v++)
					if (vf_proximity_collision_test(face0->v[v], face1, impact))
						arcsim::impacts[t].push_back(impact);
			}
		}


		// EE test : face0 edge with face1 edge
		if (arcsim::magic.eeImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int e0 = 0; e0 < face0->rEdge.size(); e0++) {
					for (int e1 = 0; e1 < face1->rEdge.size(); e1++) {
						Edge *edge0 = face0->rEdge[e0];
						Edge *edge1 = face1->rEdge[e1];
						if (ee_proximity_collision_test(edge0, edge1, impact))
							arcsim::impacts[t].push_back(impact);
					}
				}
			} else {
				for (int e0 = 0; e0 < 3; e0++)
					for (int e1 = 0; e1 < 3; e1++)
						if (ee_proximity_collision_test(face0->adje[e0], face1->adje[e1], impact))
							arcsim::impacts[t].push_back(impact);
			}
		}


		if (arcsim::magic.sim_type != "bridsonharmon") {

			// VE test : face0 vert with face1 edge
			if (arcsim::magic.veImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n = 0; n < face0->rNode.size(); n++) {
						for (int e = 0; e < face1->rEdge.size(); e++) {
							Node *node = face0->rNode[n];
							Edge *edge = face1->rEdge[e];
							if (ve_proximity_collision_test(node->verts[0], edge, impact))
								arcsim::impacts[t].push_back(impact);
						}
					}
				} else {
					for (int v = 0; v < 3; v++)
						for (int e = 0; e < 3; e++)
							if (ve_proximity_collision_test(face0->v[v], face1->adje[e], impact))
								arcsim::impacts[t].push_back(impact);
				}
			}

			// EV test : face1 vert with face0 edge
			if (arcsim::magic.evImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n = 0; n < face1->rNode.size(); n++) {
						for (int e = 0; e < face0->rEdge.size(); e++) {
							Node *node = face1->rNode[n];
							Edge *edge = face0->rEdge[e];
							if (ve_proximity_collision_test(node->verts[0], edge, impact))
								arcsim::impacts[t].push_back(impact);
						}
					}
				} else {
					for (int v = 0; v < 3; v++)
						for (int e = 0; e < 3; e++)
							if (ve_proximity_collision_test(face1->v[v], face0->adje[e], impact)) {
								impact.inverted = true;
								arcsim::impacts[t].push_back(impact);
							}
				}
			}

			// VV test : face0 vert with face1 vert
			if (arcsim::magic.vvImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n0 = 0; n0 < face0->rNode.size(); n0++) {
						for (int n1 = 0; n1 < face1->rNode.size(); n1++) {
							Node *node0 = face0->rNode[n0];
							Node *node1 = face1->rNode[n1];
							if (vv_proximity_collision_test(node0->verts[0], node1->verts[0], impact))
								arcsim::impacts[t].push_back(impact);
						}
					}
				} else {
					for (int v0 = 0; v0 < 3; v0++)
						for (int v1 = 0; v1 < 3; v1++)
							if (vv_proximity_collision_test(face0->v[v0], face1->v[v1], impact))
								arcsim::impacts[t].push_back(impact);
				}
			}
		}
	}

	void find_face_continuous_impacts(const Face *face0, const Face *face1) {

		int t = omp_get_thread_num();
		Impact impact;

		// VF test : face0 vert with face1
		if (arcsim::magic.vfImpacts) {
			for (int v = 0; v < 3; v++)
				if (vf_continuous_collision_test(face0->v[v], face1, impact))
					arcsim::impacts[t].push_back(impact);
		}

		// FV test : face1 vert with face0
		if (arcsim::magic.fvImpacts) {
			for (int v = 0; v < 3; v++)
				if (vf_continuous_collision_test(face1->v[v], face0, impact)) {
					impact.inverted = true;
					arcsim::impacts[t].push_back(impact);
				}
		}

		// EE test : face0 edge with face1 edge
		if (arcsim::magic.eeImpacts) {
			for (int e0 = 0; e0 < 3; e0++)
				for (int e1 = 0; e1 < 3; e1++)
					if (ee_continuous_collision_test(face0->adje[e0], face1->adje[e1], impact))
						arcsim::impacts[t].push_back(impact);
		}

	}

	void find_face_impacts(const Face *face0, const Face *face1) {
		int t = omp_get_thread_num();
		Impact impact;

		// FV test : face1 vert with face0
		if (arcsim::magic.fvImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int n = 0; n < face1->rNode.size(); n++) {
					Node *node = face1->rNode[n];
					if (vf_proximity_collision_test(node->verts[0], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (vf_continuous_collision_test_new(node->verts[0], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (arcsim::magic.forward_proximity &&
							   vf_proximity_collision_test(node->verts[0], face0, impact, true)) {
						arcsim::impacts[t].push_back(impact);
					}
				}
			} else {
				for (int v = 0; v < 3; v++)
					if (vf_proximity_collision_test(face1->v[v], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (vf_continuous_collision_test_new(face1->v[v], face0, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (arcsim::magic.forward_proximity &&
							   vf_proximity_collision_test(face1->v[v], face0, impact, true)) {
						arcsim::impacts[t].push_back(impact);
					}
			}
		}

		// VF test : face0 vert with face1
		if (arcsim::magic.vfImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int n = 0; n < face0->rNode.size(); n++) {
					Node *node = face0->rNode[n];
					if (vf_proximity_collision_test(node->verts[0], face1, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (vf_continuous_collision_test_new(node->verts[0], face1, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (arcsim::magic.forward_proximity &&
							   vf_proximity_collision_test(node->verts[0], face1, impact, true)) {
						arcsim::impacts[t].push_back(impact);
					}
				}
			} else {
				for (int v = 0; v < 3; v++)
					if (vf_proximity_collision_test(face0->v[v], face1, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (vf_continuous_collision_test_new(face0->v[v], face1, impact)) {
						arcsim::impacts[t].push_back(impact);
					} else if (arcsim::magic.forward_proximity &&
							   vf_proximity_collision_test(face0->v[v], face1, impact, true)) {
						arcsim::impacts[t].push_back(impact);
					}
			}
		}


		// EE test : face0 edge with face1 edge
		if (arcsim::magic.eeImpacts) {
			if (arcsim::magic.use_representative_triangles) {
				for (int e0 = 0; e0 < face0->rEdge.size(); e0++) {
					for (int e1 = 0; e1 < face1->rEdge.size(); e1++) {
						Edge *edge0 = face0->rEdge[e0];
						Edge *edge1 = face1->rEdge[e1];
						if (ee_proximity_collision_test(edge0, edge1, impact)) {
							arcsim::impacts[t].push_back(impact);
						} else if (ee_continuous_collision_test_new(edge0, edge1, impact)) {
							arcsim::impacts[t].push_back(impact);
						} else if (arcsim::magic.forward_proximity &&
								   ee_proximity_collision_test(edge0, edge1, impact, true)) {
							arcsim::impacts[t].push_back(impact);
						}
					}
				}
			} else {
				for (int e0 = 0; e0 < 3; e0++)
					for (int e1 = 0; e1 < 3; e1++)
						if (ee_proximity_collision_test(face0->adje[e0], face1->adje[e1], impact)) {
							arcsim::impacts[t].push_back(impact);
						} else if (ee_continuous_collision_test_new(face0->adje[e0], face1->adje[e1], impact)) {
							arcsim::impacts[t].push_back(impact);
						} else if (arcsim::magic.forward_proximity &&
								   ee_proximity_collision_test(face0->adje[e0], face1->adje[e1], impact, true)) {
							arcsim::impacts[t].push_back(impact);
						}
			}
		}


		if (arcsim::magic.sim_type != "bridsonharmon") {

			// VE test : face0 vert with face1 edge
			if (arcsim::magic.veImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n = 0; n < face0->rNode.size(); n++) {
						for (int e = 0; e < face1->rEdge.size(); e++) {
							Node *node = face0->rNode[n];
							Edge *edge = face1->rEdge[e];
							if (ve_proximity_collision_test(node->verts[0], edge, impact)) {
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   ve_proximity_collision_test(node->verts[0], edge, impact, true)) {
								arcsim::impacts[t].push_back(impact);
							}
						}
					}
				} else {
					for (int v = 0; v < 3; v++)
						for (int e = 0; e < 3; e++)
							if (ve_proximity_collision_test(face0->v[v], face1->adje[e], impact)) {
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   ve_proximity_collision_test(face0->v[v], face1->adje[e], impact, true)) {
								arcsim::impacts[t].push_back(impact);
							}
				}
			}

			// EV test : face1 vert with face0 edge
			if (arcsim::magic.evImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n = 0; n < face1->rNode.size(); n++) {
						for (int e = 0; e < face0->rEdge.size(); e++) {
							Node *node = face1->rNode[n];
							Edge *edge = face0->rEdge[e];
							if (ve_proximity_collision_test(node->verts[0], edge, impact)) {
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   ve_proximity_collision_test(node->verts[0], edge, impact, true)) {
								arcsim::impacts[t].push_back(impact);
							}
						}
					}
				} else {
					for (int v = 0; v < 3; v++)
						for (int e = 0; e < 3; e++)
							if (ve_proximity_collision_test(face1->v[v], face0->adje[e], impact)) {
								impact.inverted = true;
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   ve_proximity_collision_test(face1->v[v], face0->adje[e], impact, true)) {
								// impact.inverted = true;
								arcsim::impacts[t].push_back(impact);
							}
				}
			}

			// VV test : face0 vert with face1 vert
			if (arcsim::magic.vvImpacts) {
				if (arcsim::magic.use_representative_triangles) {
					for (int n0 = 0; n0 < face0->rNode.size(); n0++) {
						for (int n1 = 0; n1 < face1->rNode.size(); n1++) {
							Node *node0 = face0->rNode[n0];
							Node *node1 = face1->rNode[n1];
							if (vv_proximity_collision_test(node0->verts[0], node1->verts[0], impact)) {
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   vv_proximity_collision_test(node0->verts[0], node1->verts[0], impact, true)) {
								arcsim::impacts[t].push_back(impact);
							}
						}
					}
				} else {
					for (int v0 = 0; v0 < 3; v0++)
						for (int v1 = 0; v1 < 3; v1++)
							if (vv_proximity_collision_test(face0->v[v0], face1->v[v1], impact)) {
								arcsim::impacts[t].push_back(impact);
							} else if (arcsim::magic.forward_proximity &&
									   vv_proximity_collision_test(face0->v[v0], face1->v[v1], impact, true)) {
								arcsim::impacts[t].push_back(impact);
							}
				}
			}
		}
	}


	int solve_cubic(double a3, double a2, double a1, double a0, double t[3]);

	Vec3 pos(const Node *node, double t);

	Vec3 pos_forward(const Node *node, double t);

	bool continuous_collision_test(Impact::Type type, Node *nodes[4], Impact &impact) {

		impact.type = type;
		for (int n = 0; n < 4; n++) {
			impact.nodes[n] = nodes[n];
		}
		const Vec3 &x0 = nodes[0]->x, v0 = nodes[0]->v * dt;
		Vec3 x1 = nodes[1]->x - x0, x2 = nodes[2]->x - x0, x3 = nodes[3]->x - x0;
		Vec3 v1 = (nodes[1]->v * dt) - v0, v2 = (nodes[2]->v * dt) - v0,
				v3 = (nodes[3]->v * dt) - v0;
		double a0 = stp(x1, x2, x3),
				a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3),
				a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3),
				a3 = stp(v1, v2, v3);
		double t[4];
		int nsol = solve_cubic(a3, a2, a1, a0, t);
		t[nsol] = 1; // also check at end of timestep

		for (int i = 0; i < nsol; i++) {
			if (t[i] < 0 || t[i] > 1)
				continue;
			impact.t = t[i];
			Vec3 x0 = pos_forward(nodes[0], t[i] * dt), x1 = pos_forward(nodes[1], t[i] * dt),
					x2 = pos_forward(nodes[2], t[i] * dt), x3 = pos_forward(nodes[3], t[i] * dt);
			Vec3 &n = impact.n;
			double *w = impact.w;
			double d;
			bool inside;
			if (type == Impact::VF) {
			    w[0] = -1, w[1] = 1, w[2] = 1, w[3] = 1;
				d = signed_vf_distance(x0, x1, x2, x3, &n, w);
				inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);
			} else {// Impact::EE
			    w[0] = -1, w[1] = -1, w[2] = 1, w[3] = 1;
				d = signed_ee_distance(x0, x1, x2, x3, &n, w);
				inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);
			}
			if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > 0)
				n = -n;
			if (abs(d) < 1e-6 && inside)
				return true;
		}
		return false;
	}

	bool vf_continuous_collision_test_new(const Vert *vert, const Face *face, Impact &impact) {
		if (vert->node == face->v[0]->node || vert->node == face->v[1]->node || vert->node == face->v[2]->node)
			return false;

		Node *nodes[4];
		nodes[0] = vert->node;
		nodes[1] = face->v[0]->node;
		nodes[2] = face->v[1]->node;
		nodes[3] = face->v[2]->node;
		BOX n_box, f_box;
		n_box += nodes[0]->x;
		n_box += nodes[0]->x + nodes[0]->v * dt;
		for (int n = 0; n < 3; n++) {
			Node *node = nodes[n];
			f_box += node->x;
			f_box += node->x + node->v * dt;
		}
		if (!overlap(n_box, f_box, arcsim::thickness)) {
			return false;
		}

		bool collision = continuous_collision_test(Impact::VF, nodes, impact);
		impact.matPosA = vert->u;
		impact.matPosB = -impact.w[1] * face->v[0]->u + -impact.w[2] * face->v[1]->u + -impact.w[3] * face->v[2]->u;
		impact.verts[0] = vert;
		impact.verts[1] = face->v[0];
		impact.verts[2] = face->v[1];
		impact.verts[3] = face->v[2];
		return collision;
	}

	bool ee_continuous_collision_test_new(const Edge *edge0, const Edge *edge1, Impact &impact) {
		if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1]
			|| edge0->n[1] == edge1->n[0] || edge0->n[1] == edge1->n[1])
			return false;

		Node *nodes[4];
		nodes[0] = edge0->n[0];
		nodes[1] = edge0->n[1];
		nodes[2] = edge1->n[0];
		nodes[3] = edge1->n[1];
		BOX e_box0, e_box1;
		for (int n = 0; n < 2; n++) {
			Node *node = nodes[n];
			e_box0 += node->x;
			e_box0 += node->x + node->v * dt;
		}
		for (int n = 2; n < 4; n++) {
			Node *node = nodes[n];
			e_box1 += node->x;
			e_box1 += node->x + node->v * dt;
		}
		if (!overlap(e_box0, e_box1, arcsim::thickness)) {
			return false;
		}

		bool collision = continuous_collision_test(Impact::EE, nodes, impact);
		impact.matPosA = getMatPos(edge0, impact.w[0]);
		impact.matPosB = getMatPos(edge1, -impact.w[2]);
		impact.verts[0] = getVert(edge0, 0);
		impact.verts[1] = getVert(edge0, 1);
		impact.verts[2] = getVert(edge1, 0);
		impact.verts[3] = getVert(edge1, 1);
		return collision;
	}

// VERTEX-FACE CONTINUOUS COLLISION TEST
	bool vf_continuous_collision_test(const Vert *vert, const Face *face, Impact &impact) {

		if (vert->node == face->v[0]->node || vert->node == face->v[1]->node || vert->node == face->v[2]->node)
			return false;

		if (!overlap(node_box(vert->node, true), face_box(face, true), arcsim::thickness))
			return false;

		impact.type = Impact::VF;
		impact.nodes[0] = vert->node;
		impact.nodes[1] = face->v[0]->node;
		impact.nodes[2] = face->v[1]->node;
		impact.nodes[3] = face->v[2]->node;

		const Vec3 &x0 = impact.nodes[0]->x0;
		Vec3 x1 = impact.nodes[1]->x0 - x0;
		Vec3 x2 = impact.nodes[2]->x0 - x0;
		Vec3 x3 = impact.nodes[3]->x0 - x0;

		Vec3 v0 = impact.nodes[0]->x - x0;
		Vec3 v1 = (impact.nodes[1]->x - impact.nodes[1]->x0) - v0;
		Vec3 v2 = (impact.nodes[2]->x - impact.nodes[2]->x0) - v0;
		Vec3 v3 = (impact.nodes[3]->x - impact.nodes[3]->x0) - v0;

		double a0 = stp(x1, x2, x3);
		double a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3);
		double a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3);
		double a3 = stp(v1, v2, v3);
		double t[4];

		// Solving the cubic equation
		int nsol = solve_cubic(a3, a2, a1, a0, t);

		t[nsol] = 1; // also check at end of timestep

		// Looping through each solution from cubic solve
		for (int i = 0; i < nsol; i++) {
			if (t[i] < 0 || t[i] > 1)
				continue;
			impact.t = t[i];

			Vec3 x0 = pos(impact.nodes[0], t[i]);
			Vec3 x1 = pos(impact.nodes[1], t[i]);
			Vec3 x2 = pos(impact.nodes[2], t[i]);
			Vec3 x3 = pos(impact.nodes[3], t[i]);

			Vec3 &n = impact.n;
			double *w = impact.w;
			w[0] = -1, w[1] = 1, w[2] = 1, w[3] = 1;
			double d = signed_vf_distance(x0, x1, x2, x3, &n, w);
			bool inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);

			if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > 0)
				n = -n;
			if (abs(d) < 1e-6 && inside) {
				impact.matPosA = vert->u;
				impact.matPosB =
						-impact.w[1] * face->v[0]->u + -impact.w[2] * face->v[1]->u + -impact.w[3] * face->v[2]->u;
				impact.verts[0] = vert;
				impact.verts[1] = face->v[0];
				impact.verts[2] = face->v[1];
				impact.verts[3] = face->v[2];
				return true;
			}
		}
		return false;

	}


// EDGE-EDGE CONTINUOUS COLLISION TEST
	bool ee_continuous_collision_test(const Edge *edge0, const Edge *edge1, Impact &impact) {

		if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1]
			|| edge0->n[1] == edge1->n[0] || edge0->n[1] == edge1->n[1])
			return false;

		if (!overlap(edge_box(edge0, true), edge_box(edge1, true), arcsim::thickness))
			return false;

		impact.type = Impact::VF;
		impact.nodes[0] = edge0->n[0];
		impact.nodes[1] = edge0->n[1];
		impact.nodes[2] = edge1->n[0];
		impact.nodes[3] = edge1->n[1];

		const Vec3 &x0 = impact.nodes[0]->x0;
		Vec3 x1 = impact.nodes[1]->x0 - x0;
		Vec3 x2 = impact.nodes[2]->x0 - x0;
		Vec3 x3 = impact.nodes[3]->x0 - x0;

		Vec3 v0 = impact.nodes[0]->x - x0;
		Vec3 v1 = (impact.nodes[1]->x - impact.nodes[1]->x0) - v0;
		Vec3 v2 = (impact.nodes[2]->x - impact.nodes[2]->x0) - v0;
		Vec3 v3 = (impact.nodes[3]->x - impact.nodes[3]->x0) - v0;

		double a0 = stp(x1, x2, x3);
		double a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3);
		double a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3);
		double a3 = stp(v1, v2, v3);
		double t[4];

		// Solving the cubic equation
		int nsol = solve_cubic(a3, a2, a1, a0, t);

		t[nsol] = 1; // also check at end of timestep

		// Looping through each solution from cubic solve
		for (int i = 0; i < nsol; i++) {
			if (t[i] < 0 || t[i] > 1)
				continue;
			impact.t = t[i];

			Vec3 x0 = pos(impact.nodes[0], t[i]);
			Vec3 x1 = pos(impact.nodes[1], t[i]);
			Vec3 x2 = pos(impact.nodes[2], t[i]);
			Vec3 x3 = pos(impact.nodes[3], t[i]);

			Vec3 &n = impact.n;
			double *w = impact.w;
			w[0] = -1, w[1] = -1, w[2] = 1, w[3] = 1;
			double d = signed_ee_distance(x0, x1, x2, x3, &n, w);
			bool inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);

			if (dot(n, w[1] * v1 + w[2] * v2 + w[3] * v3) > 0)
				n = -n;
			if (abs(d) < 1e-6 && inside) {
				impact.matPosA = getMatPos(edge0, impact.w[0]);
				impact.matPosB = getMatPos(edge1, -impact.w[2]);
				impact.verts[0] = getVert(edge0, 0);
				impact.verts[1] = getVert(edge0, 1);
				impact.verts[2] = getVert(edge1, 0);
				impact.verts[3] = getVert(edge1, 1);
				return true;
			}
		}
		return false;
	}

// VERTEX-FACE PROXIMITY TEST
	bool vf_proximity_collision_test(const Vert *vert, const Face *face, Impact &impact, bool forward) {
		const Node *node = vert->node;
		if (node == face->v[0]->node || node == face->v[1]->node || node == face->v[2]->node)
			return false;
		double desired_thickness = forward ? arcsim::real_thickness : arcsim::thickness;
		BOX n_box, f_box;
		if (forward) {
			n_box += node->x + node->v * dt;
			for (int n = 0; n < 3; n++) {
				Node *f_node = face->v[n]->node;
				f_box += f_node->x + f_node->v * dt;
			}
		} else {
			n_box += node->x;
			for (int n = 0; n < 3; n++) {
				Node *f_node = face->v[n]->node;
				f_box += f_node->x;
			}
		}
		if (!overlap(n_box, f_box, desired_thickness))
			return false;

		Vec3 &n = impact.n;
		double *w = impact.w;
		Vec3 x0, x1, x2, x3;

		if (arcsim::magic.sim_type == "bridsonharmon") {
			x0 = node->x0;
			x1 = face->v[0]->node->x0;
			x2 = face->v[1]->node->x0;
			x3 = face->v[2]->node->x0;
		} else {
			if (forward) {
				x0 = node->x + node->v * dt;
				x1 = face->v[0]->node->x + face->v[0]->node->v * dt;
				x2 = face->v[1]->node->x + face->v[1]->node->v * dt;
				x3 = face->v[2]->node->x + face->v[2]->node->v * dt;
			} else {
				x0 = node->x;
				x1 = face->v[0]->node->x;
				x2 = face->v[1]->node->x;
				x3 = face->v[2]->node->x;
			}
		}

		double d = signed_vf_distance(x0, x1, x2, x3, &n, w);
		bool inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);


		if (fabs(d) < desired_thickness && inside) {

			if (dot(n, x2 - x0) > 0) {
				n = -n;
			}


			impact.type = Impact::VF;
			impact.matPosA = vert->u;
			impact.matPosB = -w[1] * face->v[0]->u + -w[2] * face->v[1]->u + -w[3] * face->v[2]->u;
			impact.nodes[0] = (Node *) node;
			impact.nodes[1] = (Node *) face->v[0]->node;
			impact.nodes[2] = (Node *) face->v[1]->node;
			impact.nodes[3] = (Node *) face->v[2]->node;
			impact.verts[0] = vert;
			impact.verts[1] = face->v[0];
			impact.verts[2] = face->v[1];
			impact.verts[3] = face->v[2];

			return true;
		}

		return false;

	}


	Vec2 getVert0Pos(const Edge *edge) {

		for (int s = 0; s < 2; s++) {

			if (!edge->adjf[s]) {
				continue;
			}

			Vert *v0 = edge_vert(edge, s, 0);
			Vert *v1 = edge_vert(edge, s, 1);

			return edge_vert(edge, s, 0)->u;

		}

        assert(false);
        return Vec2(0,0);
	}

	Vec2 getVert1Pos(const Edge *edge) {

		for (int s = 0; s < 2; s++) {

			if (!edge->adjf[s]) {
				continue;
			}

			Vert *v0 = edge_vert(edge, s, 0);
			Vert *v1 = edge_vert(edge, s, 1);

			return edge_vert(edge, s, 1)->u;

		}

        assert(false);
        return Vec2(0,0);

	}

	Vert *getVert(const Edge *edge, int i) {
		for (int s = 0; s < 2; s++) {
			if (!edge->adjf[s]) {
				continue;
			}
			return edge_vert(edge, s, i);
		}
        assert(false);
        return NULL;
	}

// Check if the two edges belong to adjacent faces
	bool is_adjacent(const Edge *edge0, const Edge *edge1) {
		for (int n = 0; n < 2; n++) {
			for (int f = 0; f < 2; f++) {
				Node *node = edge0->n[n];
				Face *face = edge1->adjf[f];
				if (face && (node == face->v[0]->node || node == face->v[1]->node || node == face->v[2]->node)) {
					return true;
				}
				node = edge1->n[n];
				face = edge0->adjf[f];
				if (face && (node == face->v[0]->node || node == face->v[1]->node || node == face->v[2]->node)) {
					return true;
				}
			}
		}
		return false;
	}

// EDGE-EDGE PROXIMITY TEST
	bool ee_proximity_collision_test(const Edge *edge0, const Edge *edge1, Impact &impact, bool forward) {
		if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1]
			|| edge0->n[1] == edge1->n[0] || edge0->n[1] == edge1->n[1])
			return false;
		double desired_thickness = forward ? arcsim::real_thickness : arcsim::thickness;
		BOX e_box0, e_box1;
		if (forward) {
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge0->n[n];
				e_box0 += e_node->x + e_node->v * dt;
			}
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge1->n[n];
				e_box1 += e_node->x + e_node->v * dt;
			}
		} else {
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge0->n[n];
				e_box0 += e_node->x;
			}
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge1->n[n];
				e_box1 += e_node->x;
			}
		}
		if (!overlap(e_box0, e_box1, desired_thickness))
			return false;

		Vec3 &n = impact.n;
		double *w = impact.w;
		Vec3 x0, x1, x2, x3;
		if (arcsim::magic.sim_type == "bridsonharmon") {
			x0 = edge0->n[0]->x0;
			x1 = edge0->n[1]->x0;
			x2 = edge1->n[0]->x0;
			x3 = edge1->n[1]->x0;
		} else {
			if (forward) {
				x0 = edge0->n[0]->x + edge0->n[0]->v * dt;
				x1 = edge0->n[1]->x + edge0->n[1]->v * dt;
				x2 = edge1->n[0]->x + edge1->n[0]->v * dt;
				x3 = edge1->n[1]->x + edge1->n[1]->v * dt;
			} else {
				x0 = edge0->n[0]->x;
				x1 = edge0->n[1]->x;
				x2 = edge1->n[0]->x;
				x3 = edge1->n[1]->x;
			}
		}
		double d = signed_ee_distance(x0, x1, x2, x3, &n, w);
		bool inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);

		bool is_adj = is_adjacent(edge0, edge1);
		is_adj = false;
		bool e0_in_e1_wedge = is_adj || inWedge(w[0] * x0 + w[1] * x1, (Edge *) edge1, forward);
		bool e1_in_e0_wedge = is_adj || inWedge(-w[2] * x2 + -w[3] * x3, (Edge *) edge0, forward);

		if (fabs(d) < desired_thickness && (e0_in_e1_wedge || e1_in_e0_wedge) && inside) {

			Vec3 x10 = edge0->n[1]->x0 - edge0->n[0]->x0;
			Vec3 y10 = edge1->n[1]->x0 - edge1->n[0]->x0;
			Vec3 nComputed = cross(normalize(x10), normalize(y10));

			n = w[0] * x0 + w[1] * x1 + w[2] * x2 + w[3] * x3;
			impact.n = normalize(n);

			impact.type = Impact::EE;
			impact.nodes[0] = (Node *) edge0->n[0];
			impact.nodes[1] = (Node *) edge0->n[1];
			impact.nodes[2] = (Node *) edge1->n[0];
			impact.nodes[3] = (Node *) edge1->n[1];
			impact.matPosA = getMatPos(edge0, w[0]);
			impact.matPosB = getMatPos(edge1, -w[2]);
			impact.verts[0] = getVert(edge0, 0);
			impact.verts[1] = getVert(edge0, 1);
			impact.verts[2] = getVert(edge1, 0);
			impact.verts[3] = getVert(edge1, 1);


			return true;
		}

		return false;
	}


// VERTEX-EDGE PROXIMITY TEST
	bool ve_proximity_collision_test(const Vert *vert, const Edge *edge, Impact &impact, bool forward) {
		const Node *node = vert->node;
		if (node == edge->n[0] || node == edge->n[1])
			return false;
		double desired_thickness = forward ? arcsim::real_thickness : arcsim::thickness;
		BOX n_box, e_box;
		if (forward) {
			n_box += node->x + node->v * dt;
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge->n[n];
				e_box += e_node->x + e_node->v * dt;
			}
		} else {
			n_box += node->x;
			for (int n = 0; n < 2; n++) {
				Node *e_node = edge->n[n];
				e_box += e_node->x;
			}
		}
		if (!overlap(n_box, e_box, desired_thickness))
			return false;

		Vec3 x, y0, y1;
		if (forward) {
			x = node->x;
			y0 = edge->n[0]->x + edge->n[0]->v * dt;
			y1 = edge->n[1]->x + edge->n[1]->v * dt;
		} else {
			x = node->x;
			y0 = edge->n[0]->x;
			y1 = edge->n[1]->x;
		}
		double d = 1.e8;
		Vec3 &n = impact.n;
		double wx;
		double wy0;
		double wy1;

		set_unsigned_ve_distance(x, y0, y1, &d, &n, &wx, &wy0, &wy1);

		bool nodeInWedge = inWedge(x, (Edge *) edge, forward);
		bool edgeInCone = wy0 < 1.0 && wy1 < 1.0 ? inCone(wy0 * y0 + wy1 * y1, (Node *) node, forward) : false;

		if (fabs(d) < desired_thickness && nodeInWedge && edgeInCone) {

			impact.type = Impact::VE;
			impact.nodes[0] = (Node *) node;
			impact.nodes[1] = (Node *) edge->n[0];
			impact.nodes[2] = (Node *) edge->n[1];
			impact.nodes[3] = 0;
			impact.matPosA = vert->u;
			impact.matPosB = getMatPos(edge, wy0);
			impact.verts[0] = vert;
			impact.verts[1] = getVert(edge, 0);
			impact.verts[2] = getVert(edge, 1);
			impact.verts[3] = NULL;

			impact.w[0] = wx;
			impact.w[1] = wy0;
			impact.w[2] = wy1;
			impact.w[3] = 0;

			return true;

		}

		return false;
	}

// VERTEX-VERTEX PROXIMITY TEST
	bool vv_proximity_collision_test(const Vert *vert0, const Vert *vert1, Impact &impact, bool forward) {
		if (arcsim::magic.sim_type == "bridsonharmon") {
			return false;
		}
		double desired_thickness = forward ? arcsim::real_thickness : arcsim::thickness;
		const Node *node0 = vert0->node;
		const Node *node1 = vert1->node;
		if (node0 == node1)
			return false;
		BOX n_box0, n_box1;
		if (forward) {
			n_box0 += node0->x + node0->v * dt;
			n_box1 += node1->x + node1->v * dt;
		} else {
			n_box0 += node0->x;
			n_box1 += node1->x;
		}
		if (!overlap(n_box0, n_box1, desired_thickness))
			return false;

		double d = norm(node1->x - node0->x);
		bool node0InsideCone = inCone(node0->x, (Node *) node1, forward);
		bool node1InsideCone = inCone(node1->x, (Node *) node0, forward);

		if (node0InsideCone && node1InsideCone && d < desired_thickness) {

			impact.type = Impact::VV;
			impact.nodes[0] = (Node *) node0;
			impact.nodes[1] = (Node *) node1;
			impact.nodes[2] = NULL;
			impact.nodes[3] = NULL;
			impact.n = normalize(node0->x - node1->x);  // todo: might need to flip sign on n?
			impact.w[0] = 1;
			impact.w[1] = 1;
			impact.w[2] = 0;
			impact.w[3] = 0;
			impact.matPosA = vert0->u;
			impact.matPosB = vert1->u;
			impact.verts[0] = vert0;
			impact.verts[1] = vert1;
			impact.verts[2] = NULL;
			impact.verts[3] = NULL;

			return true;

		}

		return false;
	}


	Vec3 pos(const Node *node, double t) {
		return node->x0 + t * (node->x - node->x0);
	}

	Vec3 pos_forward(const Node *node, double t) {
		return node->x + t * node->v;
	}

// Solving cubic equations

	double newtons_method(double a, double b, double c, double d, double x0,
						  int init_dir);

// solves a x^3 + b x^2 + c x + d == 0
	int solve_cubic(double a, double b, double c, double d, double x[3]) {
		double xc[2];
		int ncrit = solve_quadratic(3 * a, 2 * b, c, xc);
		if (ncrit == 0) {
			x[0] = newtons_method(a, b, c, d, xc[0], 0);
			return 1;
		} else if (ncrit == 1) {// cubic is actually quadratic
			return solve_quadratic(b, c, d, x);
		} else {
			double yc[2] = {d + xc[0] * (c + xc[0] * (b + xc[0] * a)),
							d + xc[1] * (c + xc[1] * (b + xc[1] * a))};
			int i = 0;
			if (yc[0] * a >= 0)
				x[i++] = newtons_method(a, b, c, d, xc[0], -1);
			if (yc[0] * yc[1] <= 0) {
				int closer = abs(yc[0]) < abs(yc[1]) ? 0 : 1;
				x[i++] = newtons_method(a, b, c, d, xc[closer], closer == 0 ? 1 : -1);
			}
			if (yc[1] * a <= 0)
				x[i++] = newtons_method(a, b, c, d, xc[1], 1);
			return i;
		}
	}

	double newtons_method(double a, double b, double c, double d, double x0,
						  int init_dir) {
		if (init_dir != 0) {
			// quadratic approximation around x0, assuming y' = 0
			double y0 = d + x0 * (c + x0 * (b + x0 * a)),
					ddy0 = 2 * b + x0 * (6 * a);
			x0 += init_dir * sqrt(abs(2 * y0 / ddy0));
		}
		for (int iter = 0; iter < 100; iter++) {
			double y = d + x0 * (c + x0 * (b + x0 * a));
			double dy = c + x0 * (2 * b + x0 * 3 * a);
			if (dy == 0)
				return x0;
			double x1 = x0 - y / dy;
			if (abs(x0 - x1) < 1e-6)
				return x0;
			x0 = x1;
		}
		return x0;
	}

}


