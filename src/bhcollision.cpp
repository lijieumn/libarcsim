#include "bhcollision.hpp"
#include "simulation.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "optimization.hpp"
#include "timer.hpp"
#include "proximity.hpp"
#include "blockvectors.hpp"
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <algorithm>    // std::random_shuffle
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

using namespace std;


namespace arcsim {

	static const int max_iter = 100;

	static double obs_mass;
	static bool deform_obstacles;

	static vector<Vec3> xold;
	static vector<Vec3> xold_obs;

	struct ImpactZone {
		vector<Node *> nodes;
		vector<Impact> impacts;
		bool active;
	};

// returns pair of (i) is_free(vert), and
// (ii) index of mesh in ::meshes or ::obs_meshes that contains vert
	pair<bool, int> find_in_meshes(const Node *node) {
		int m = find_mesh(node, *arcsim::meshes);
		if (m != -1)
			return make_pair(true, m);
		else
			return make_pair(false, find_mesh(node, *arcsim::obs_meshes));
	}

	void update_active(const vector<AccelStruct *> &accs,
					   const vector<AccelStruct *> &obs_accs,
					   const vector<ImpactZone *> &zones);

	double get_mass(const Node *node) { return is_free(node) ? node->m : obs_mass; }

	vector<Impact> independent_impacts(const vector<Impact> &impacts);

	void add_impacts(const vector<Impact> &impacts, vector<ImpactZone *> &zones);

	void reset_zone_velocities(ImpactZone *zone);

	void apply_inelastic_projection(ImpactZone *zone, const vector<Constraint *> &cons);

	void apply_harmon_friction(ImpactZone *zone, const vector<Constraint *> &cons);


// Bridson Repulsion
	void apply_repulsion(std::vector<Impact> &impacts, Mesh *mesh, double dt) {

		bool verbose = true;

		// Cloth thickness
		double h = arcsim::magic.repulsion_thickness;

		// Repulsion spring stiffness
		double k = arcsim::magic.collision_stiffness;

		// Applying repulsion to all the proximity impacts
		for (int i = 0; i < impacts.size(); i++) {

			Impact &impact = impacts[i];

			// Node positions at beginning of time step
			Vec3 x[4];
			x[0] = impact.nodes[0]->x0;
			x[1] = impact.nodes[1]->x0;
			x[2] = impact.nodes[2]->x0;
			x[3] = impact.nodes[3]->x0;

			// Node midstep velocities
			Vec3 v[4];
			v[0] = impact.nodes[0]->vTildeMid;  // could init to vBarMid instead and then run in parallel
			v[1] = impact.nodes[1]->vTildeMid;  // bridson says if you do that, probably need to reduce the
			v[2] = impact.nodes[2]->vTildeMid;  // applied impulses
			v[3] = impact.nodes[3]->vTildeMid;

			// Node masses
			double m[4];
			m[0] = impact.nodes[0]->m;
			m[1] = impact.nodes[1]->m;
			m[2] = impact.nodes[2]->m;
			m[3] = impact.nodes[3]->m;

			// Impact normal
			Vec3 n = impact.n;

			// dv from repulsion impulse (to be computed)
			Vec3 dv[4];
			dv[0] = Vec3(0);
			dv[1] = Vec3(0);
			dv[2] = Vec3(0);
			dv[3] = Vec3(0);


			// VF or EE impact?
			switch (impact.type) {

				case Impact::VF: {

					bool inverted = false;

					bool first = impact.verts[0]->component != -1;
					bool second = impact.verts[1]->component != -1;
					bool third = impact.verts[2]->component != -1;
					bool fourth = impact.verts[3]->component != -1;

					impact.inverted = !(first && !second && !third && !fourth);

					if (!inverted) {

						// Barycentric data
						double w[4];
						w[0] = impact.w[0];
						w[1] = -impact.w[1];
						w[2] = -impact.w[2];
						w[3] = -impact.w[3];

						// Vert and point in face positions
						Vec3 xVert = x[0];
						Vec3 xFace = w[1] * x[1] + w[2] * x[2] + w[3] * x[3];

						double xN = dot(xVert - xFace, n);

						if (i == 0 && verbose) {
							std::cout << "w[0]: " << w[0] << std::endl;
							std::cout << "w[1]: " << w[1] << std::endl;
							std::cout << "w[2]: " << w[2] << std::endl;
							std::cout << "w[3]: " << w[3] << std::endl;
							std::cout << "xVert: " << xVert << std::endl;
							std::cout << "xFace: " << xFace << std::endl;
							std::cout << "xN: " << xN << std::endl;
							std::cout << "n: " << n << std::endl;
						}

						// Vert and point in face velocities
						Vec3 vVert = v[0];
						Vec3 vFace = w[1] * v[1] + w[2] * v[2] + w[3] * v[3];
						double vN = dot(vVert - vFace, n);

						if (i == 0 && verbose) {
							std::cout << "vVert: " << vVert << std::endl;
							std::cout << "vFace: " << vFace << std::endl;
							std::cout << "vN: " << vN << std::endl;
						}

						// Overlap amount
						double d = h - xN;

						if (i == 0 && verbose) {
							std::cout << "d: " << d << std::endl;
						}


						// Effective mass for the collision
						double mEffective =
								1.0 / (1.0 / m[0] + w[1] * w[1] / m[1] + w[2] * w[2] / m[2] + w[3] * w[3] / m[3]);

						if (i == 0 && verbose) {
							std::cout << "mEffective: " << mEffective << std::endl;
							std::cout << "mAlt: " << (w[1] * m[1] + w[2] * m[2] + w[3] * m[3]) /
													 (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
						}

						// Candidate repulsion impulses
						double Jr1 = dt * k * d;
						double Jr2 = mEffective * (0.1 * d / dt - vN);

						if (i == 0 && verbose) {
							std::cout << "0.1*d/dt: " << 0.1 * d / dt << std::endl;
							std::cout << "0.1*d/dt - vN: " << 0.1 * d / dt - vN << std::endl;
							std::cout << "Jr1: " << Jr1 << std::endl;
							std::cout << "Jr2: " << Jr2 << std::endl;
						}

						// If vN < 0.1*d/dt, then impulse is negative of the minimum of the two candidates. Else no impulse.
						double Jr = 0;
						bool isSpring;
						if (vN < 0.1 * d / dt) {
							if (Jr1 < Jr2) {
								Jr = -Jr1;
								isSpring = true;
							} else {
								Jr = -Jr2;
								isSpring = false;
							}
						}


						Jr = 2.0 * Jr / (1.0 + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);


						if (i == 0 && verbose) {
							std::cout << "Jr: " << Jr << std::endl;
						}

						// Computing dv from repulsion impulse

						VisualDebugger *vd = VisualDebugger::getInstance();
						if (!impact.inverted) {
							vd->addVisualPoint3(xVert, Vec3(0, 0, 1), 'v');
							vd->addVisualPoint3(xFace, Vec3(0, 0, 1), 'v');
							vd->addVisualLine3(xVert, xFace, Vec3(0, 0, 1), 'v');
							vd->addVisualLine3(xFace, xFace - 0.05 * n, Vec3(0, 1, 1), 'v');
							vd->addVisualLine3(xVert, xVert + 0.05 * n, Vec3(0, 1, 1), 'v');
						} else {
							vd->addVisualPoint3(xVert, Vec3(0, 1, 0), 'f');
							vd->addVisualPoint3(xFace, Vec3(0, 1, 0), 'f');
							vd->addVisualLine3(xVert, xFace, Vec3(0, 1, 0), 'f');
							vd->addVisualLine3(xFace, xFace - 0.05 * n, Vec3(0, 1, 1), 'f');
							vd->addVisualLine3(xVert, xVert + 0.05 * n, Vec3(0, 1, 1), 'f');
						}

						// This works for static obstacles
						Vec3 dvt = (Jr / mEffective) * n;
						dv[0] = -dvt;
						dv[1] = w[1] * dvt;
						dv[2] = w[2] * dvt;
						dv[3] = w[3] * dvt;

						// Bridson method which isn't giving good results
						//dv[0] = -(Jr/m[0])*n;
						//dv[1] = w[1]*(Jr/ m[1] )*n;
						//dv[2] = w[2]*(Jr/ m[2] )*n;
						//dv[3] = w[3]*(Jr/ m[3] )*n;

						if (!impact.inverted) {
							for (int k = 0; k < 4; k++) {
								vd->addVisualLine3(x[k], x[k] + 10.0 * dv[k], Vec3(0, 0, 1), '5');
							}
						} else {
							for (int k = 0; k < 4; k++) {
								vd->addVisualLine3(x[k], x[k] + 10.0 * dv[k], Vec3(0, 1, 0), '6');
							}
						}


						if (i == 0 && verbose) {
							std::cout << "dv[0]: " << dv[0] << std::endl;
							std::cout << "dv[1]: " << dv[1] << std::endl;
							std::cout << "dv[2]: " << dv[2] << std::endl;
							std::cout << "dv[3]: " << dv[3] << std::endl;
						}

						// Relative change in velocity from repulsion forces, in the normal direction
						impact.dvN = dot(dv[0] - (w[1] * dv[1] + w[2] * dv[2] + w[3] * dv[3]), n);

						if (i == 0 && verbose) {
							std::cout << "impact.dvN: " << impact.dvN << std::endl;
							std::cout << "-----------\n";
						}

					}


					break; // break Impact::VF
				}

				case Impact::EE: {

					// Barycentric data
					double w[4];
					w[0] = impact.w[0];
					w[1] = impact.w[1];
					w[2] = -impact.w[2];
					w[3] = -impact.w[3];

					double a = w[1];
					double b = w[3];

					// Point in edge0 and point in edge1 positions
					Vec3 xEdge0 = x[0] + a * (x[1] - x[0]);
					Vec3 xEdge1 = x[2] + b * (x[3] - x[2]);
					double xN = dot(xEdge0 - xEdge1, n);

					if (verbose) {

						std::cout << "xEdge0: " << xEdge0 << std::endl;
						std::cout << "xEdge1: " << xEdge1 << std::endl;
						std::cout << "n: " << n << std::endl;
						std::cout << "xN: " << xN << std::endl;

					}

					// Point in edge0 and point in edge1 velocities
					Vec3 vEdge0 = v[0] + a * (v[1] - v[0]);
					Vec3 vEdge1 = v[2] + b * (v[3] - v[2]);
					double vN = dot(vEdge0 - vEdge1, n);

					if (verbose) {

						std::cout << "vEdge0: " << vEdge0 << std::endl;
						std::cout << "vEdge1: " << vEdge1 << std::endl;
						std::cout << "vN: " << vN << std::endl;

					}

					// Overlap amount
					double d = h - xN;

					if (verbose) {

						std::cout << "d: " << d << std::endl;
					}

					// Effective mass for the collision
					double mEffective = 1.0 / (((1.0 - a) * (1.0 - a)) / m[0] + (a * a) / m[1] +
											   ((1.0 - b) * (1.0 - b)) / m[2] + (b * b) / m[3]);


					if (verbose) {
						std::cout << "mEffective: " << mEffective << std::endl;
					}

					// Candidate repulsion impulses
					double Jr1 = dt * k * d;
					double Jr2 = 0.25 * mEffective * (0.1 * d / dt - vN);

					if (verbose) {

						std::cout << "0.1*d/dt: " << 0.1 * d / dt << std::endl;
						std::cout << "0.1*d/dt - vN: " << 0.1 * d / dt - vN << std::endl;
						std::cout << "Jr1: " << Jr1 << std::endl;
						std::cout << "Jr2: " << Jr2 << std::endl;

					}

					// Impulse is negative of the minimum of the two candidates
					//double Jr = -1.0 * std::min(Jr1,Jr2);

					double Jr = 0;
					bool isSpring;
					if (vN < 0.1 * d / dt) {
						if (Jr1 < Jr2) {
							Jr = -Jr1;
							isSpring = true;
						} else {
							Jr = -Jr2;
							isSpring = false;
						}
					}

					if (verbose) {

						std::cout << "Jr: " << Jr << std::endl;

					}

					Jr = 2.0 * Jr / (w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);

					mesh->bridson_ee_impacts.push_back(std::make_tuple(xEdge0, xEdge1, n, -Jr, isSpring));
					VisualDebugger *vd = VisualDebugger::getInstance();
					vd->addVisualPoint3(xEdge0, Vec3(1, 0, 0), 'e');
					vd->addVisualPoint3(xEdge1, Vec3(1, 0, 0), 'e');
					vd->addVisualLine3(xEdge0, xEdge1, Vec3(1, 0, 0), 'e');
					vd->addVisualLine3(xEdge1, xEdge1 - 0.05 * n, Vec3(0, 1, 1), 'e');
					vd->addVisualLine3(xEdge0, xEdge0 + 0.05 * n, Vec3(0, 1, 1), 'e');


					// This works for static obstacles
					//Vec3 dvt = (Jr / mEffective) * n;
					//dv[0] = -(1.0-a)*dvt;
					//dv[1] = -(a)*dvt;
					//dv[2] = (1.0-b)*dvt;
					//dv[3] = (b)*dvt;


					// Bridson strategy
					dv[0] = -(1.0 - a) * (Jr / m[0]) * n;
					dv[1] = -(a) * (Jr / m[1]) * n;
					dv[2] = (1.0 - b) * (Jr / m[2]) * n;
					dv[3] = (b) * (Jr / m[3]) * n;

					if (verbose) {

						std::cout << "dv[0]: " << dv[0] << std::endl;
						std::cout << "dv[1]: " << dv[1] << std::endl;
						std::cout << "dv[2]: " << dv[2] << std::endl;
						std::cout << "dv[3]: " << dv[3] << std::endl;

					}

					// Relative change in velocity from repulsion forces, in the normal direction
					impact.dvN = dot((dv[0] + a * (dv[1] - dv[0])) - (dv[2] + b * (dv[3] - dv[2])), n);

					if (verbose) {

						std::cout << "impact.dvN: " << impact.dvN << std::endl;

					}

					break;

				} // end of case Impact::EE

			}

			// Applying impulse to the nodes
			impact.nodes[0]->vTildeMid += dv[0];
			impact.nodes[1]->vTildeMid += dv[1];
			impact.nodes[2]->vTildeMid += dv[2];
			impact.nodes[3]->vTildeMid += dv[3];

			if (verbose && i == 0) {

				std::cout << "after applying impulse to nodes..\n";
				std::cout << "vMid[0]: " << impact.nodes[0]->vTildeMid << std::endl;
				std::cout << "vMid[1]: " << impact.nodes[1]->vTildeMid << std::endl;
				std::cout << "vMid[2]: " << impact.nodes[2]->vTildeMid << std::endl;
				std::cout << "vMid[3]: " << impact.nodes[3]->vTildeMid << std::endl;
				std::cout << "---------------------------------\n";

			}


		}
	}

// Bridson Friction
	void apply_friction(std::vector<Impact> &impacts) {

		double mu = arcsim::magic.friction_coeff;

		for (int i = 0; i < impacts.size(); i++) {

			Impact &impact = impacts[i];

			// Node midstep velocities
			Vec3 v[4];
			v[0] = impact.nodes[0]->vTildeMid;
			v[1] = impact.nodes[1]->vTildeMid;
			v[2] = impact.nodes[2]->vTildeMid;
			v[3] = impact.nodes[3]->vTildeMid;

			// Node masses
			double m[4];
			m[0] = impact.nodes[0]->m;
			m[1] = impact.nodes[1]->m;
			m[2] = impact.nodes[2]->m;
			m[3] = impact.nodes[3]->m;

			// Impact normal
			Vec3 n = impact.n;

			// dv from friction impulse (to be computed)
			Vec3 dv[4];
			dv[0] = Vec3(0);
			dv[1] = Vec3(0);
			dv[2] = Vec3(0);
			dv[3] = Vec3(0);

			// Change in relative velocity in normal direction from the associated repulsion impulse
			double dvN = impact.dvN;

			switch (impact.type) {

				case Impact::VF: {

					// Barycentric data
					double w[4];
					w[0] = impact.w[0];
					w[1] = -impact.w[1];
					w[2] = -impact.w[2];
					w[3] = -impact.w[3];

					// Midstep updated velocities of the vert and point in face (post-repulsion)
					Vec3 vVert = v[0];
					Vec3 vFace = w[1] * v[1] + w[2] * v[2] + w[3] * v[3];

					// Midstep velocity in the transverse direction
					Vec3 vT = (vVert - vFace) - dot(vVert - vFace, n) * n;

					// Negative unit transverse direction (this is the direction friction is applied in)
					Vec3 u = -1.0 * normalize(vT);

					// Effective mass for the collision
					double mEffective =
							1.0 / (1.0 / m[0] + w[1] * w[1] / m[1] + w[2] * w[2] / m[2] + w[3] * w[3] / m[3]);

					// Candidate friction impulses
					double Jf1 = mEffective * mu * dvN;
					double Jf2 = mEffective * norm(vT);
					double Jf = -1.0 * std::min(Jf1, Jf2);

					//Jf =  2.0*Jf / ( 1.0 + w[1]*w[1] + w[2]*w[2] + w[3]*w[3] );

					// This works for static obstacles
					//Vec3 dvt = (Jf / mEffective) * u;
					//dv[0] = -dvt;
					//dv[1] = w[1]*dvt;
					//dv[2] = w[2]*dvt;
					//dv[3] = w[3]*dvt;

					// Computing dv from friction impulse
					dv[0] = -(Jf / m[0]) * u;
					dv[1] = w[1] * (Jf / m[1]) * u;
					dv[2] = w[2] * (Jf / m[2]) * u;
					dv[3] = w[3] * (Jf / m[3]) * u;

					break;

				} // end of case Impact::VF

				case Impact::EE: {

					// Barycentric data
					/*double w[4];
                    w[0] = impact.w[0];
                    w[1] = impact.w[1];
                    w[2] = -impact.w[2];
                    w[3] = -impact.w[3];

                    // Midstep updated velocities of the points in edges (post-repulsion)
                    Vec3 vEdge0 = w[0]*v[0] + w[1]*v[1];
                    Vec3 vEdge1 = w[2]*v[2] + w[3]*v[3];

                    // Midstep velocity in the transverse direction
                    Vec3 vT = (vEdge0 - vEdge1) - dot(vEdge0 - vEdge1,n)*n;

                    // Negative unit transverse direction (this is the direction friction is applied in)
                    Vec3 u = -1.0 * normalize( vT );

                    // Effective mass for the collision
                    double mEffective = 1.0 / ( w[0]*w[0]/m[0] + w[1]*w[1]/m[1] + w[2]*w[2]/m[2] + w[3]*w[3]/m[3] );

                    // Candidate friction impulses
                    double Jf1 = mEffective*mu*dvN;
                    double Jf2 = mEffective*norm(vT);
                    double Jf = -1.0*min(Jf1,Jf2);


                    Jf =  2.0*Jf / ( 1.0 + w[1]*w[1] + w[2]*w[2] + w[3]*w[3] );

                    // This works for static obstacles
                    Vec3 dvt = (Jf / mEffective) * u;
                    dv[0] = -dvt;
                    dv[1] = w[1]*dvt;
                    dv[2] = w[2]*dvt;
                    dv[3] = w[3]*dvt;

                    // Computing dv from friction impulse
                    dv[0] = -w[0]*(Jf/m[0])*u;
                    dv[1] = -w[1]*(Jf/m[1])*u;
                    dv[2] = w[2]*(Jf/m[2])*u;
                    dv[3] = w[3]*(Jf/m[3])*u;

                    break; */

				} // end of case Impact::EE

			}

			// Updating midstep velocities based on friction impulses
			impact.nodes[0]->vTildeMid += dv[0];
			impact.nodes[1]->vTildeMid += dv[1];
			impact.nodes[2]->vTildeMid += dv[2];
			impact.nodes[3]->vTildeMid += dv[3];

		}

	}


	enum CompareType {
		DIFFERENT,
		EQUAL,
		EARLIER,
		LATER
	};

	typedef CompareType (*CompareCallback)(const Impact &imp1, const Impact &imp2);


	CompareType impacts_equal2(const Impact &imp1, const Impact &imp2) {
		if (imp1.type != imp2.type) {
			return DIFFERENT;
		}
		if (imp1.type == Impact::EE) {
			bool edge1_equal = false, edge2_equal = false;
			if (imp1.nodes[0] == imp2.nodes[0] && imp1.nodes[1] == imp2.nodes[1]) {
				edge1_equal = true;
			} else if (imp1.nodes[0] == imp2.nodes[1] && imp1.nodes[1] == imp2.nodes[0]) {
				edge1_equal = true;
			}
			if (!edge1_equal) {
				return DIFFERENT;
			}
			if (imp1.nodes[2] == imp2.nodes[2] && imp1.nodes[3] == imp2.nodes[3]) {
				edge2_equal = true;
			} else if (imp1.nodes[2] == imp2.nodes[3] && imp1.nodes[3] == imp2.nodes[2]) {
				edge2_equal = true;
			}
			if (edge2_equal) {
				return EQUAL;
			}
		} else {
			if (imp1.nodes[0] != imp2.nodes[0]) {
				return DIFFERENT;
			}
			bool nodes_equal[3] = {false, false, false};
			for (int n = 1; n <= 3; n++) {
				if (imp1.nodes[n] == imp2.nodes[n] || imp1.nodes[n] == imp2.nodes[n % 3 + 1]
					|| imp1.nodes[n] == imp2.nodes[(n + 1) % 3 + 1]) {
					nodes_equal[n - 1] = true;
				} else {
					return DIFFERENT;
				}
			}
			if (nodes_equal[0] && nodes_equal[1] && nodes_equal[2]) {
				return EQUAL;
			}
		}
		return DIFFERENT;
	}

	CompareType vf_vertex_equal2(const Impact &imp1, const Impact &imp2) {
		if (imp1.type == Impact::VF && imp2.type == Impact::VF) {
			if (imp1.nodes[0] == imp2.nodes[0]) {
				return imp1.t < imp2.t ? EARLIER : LATER;
			}
		}
		return DIFFERENT;
	}


/* Removes duplicate impacts */
	void prune_impacts(std::vector<Impact> &impacts, CompareCallback compare) {

		std::vector<int> markToDelete;
		vector<int> markToRemain;

		for (int i = 0; i < impacts.size(); i++) {
			Impact &imp1 = impacts[i];
			CompareType type = DIFFERENT;
			for (int j = 0; j < markToRemain.size(); j++) {
				Impact &imp2 = impacts[markToRemain[j]];
				type = compare(imp1, imp2);
				if (type == EARLIER) {
					markToRemain.erase(markToRemain.begin() + j);
					break;
				} else if (type == LATER || type == EQUAL) {
					break;
				}
			}
			if (type == DIFFERENT || type == EARLIER) {
				markToRemain.push_back(i);
			}
		}

		int delete_count = 0;
		for (int i = 0; i < impacts.size(); i++) {
			if (markToRemain.size() <= i - delete_count || markToRemain[i - delete_count] != i) {
				delete_count++;
				markToDelete.push_back(i);
			}
		}

		for (int i = int(markToDelete.size()) - 1; i >= 0; i--) {
			impacts.erase(impacts.begin() + markToDelete[i]);
		}

		markToDelete.clear();

	}


	void collision_response_arcsim(vector<Mesh *> &meshes, const vector<Constraint *> &cons,
								   const vector<Mesh *> &obs_meshes) {
		arcsim::meshes = &meshes;
		arcsim::obs_meshes = &obs_meshes;
		arcsim::xold = node_positions(meshes);
		arcsim::xold_obs = node_positions(obs_meshes);
		vector<AccelStruct *> accs = create_accel_structs(meshes, true),
				obs_accs = create_accel_structs(obs_meshes, true);
		vector<ImpactZone *> zones;
		arcsim::obs_mass = 1e3;
		int iter;
		for (int deform = 0; deform <= 1; deform++) {
			arcsim::deform_obstacles = deform;
			zones.clear();
			for (iter = 0; iter < max_iter; iter++) {
				if (!zones.empty())
					update_active(accs, obs_accs, zones);
				vector<Impact> impacts = find_continuous_impacts(accs, obs_accs);
				impacts = independent_impacts(impacts);
				if (impacts.empty())
					break;
				add_impacts(impacts, zones);
				for (int z = 0; z < zones.size(); z++) {
					ImpactZone *zone = zones[z];
					apply_inelastic_projection(zone, cons);
				}
				for (int a = 0; a < accs.size(); a++)
					update_accel_struct(*accs[a]);
				for (int a = 0; a < obs_accs.size(); a++)
					update_accel_struct(*obs_accs[a]);
				if (deform_obstacles)
					arcsim::obs_mass /= 2;
			}
			if (iter < max_iter) // success!
				break;
		}
		if (iter == max_iter) {
			cerr << "Collision resolution failed to converge!" << endl;
			debug_save_meshes(meshes, "meshes");
			debug_save_meshes(obs_meshes, "obsmeshes");
			exit(1);
		}
		for (int m = 0; m < meshes.size(); m++) {
			compute_ws_data(*meshes[m]);
			update_x0(*meshes[m]);
		}
		for (int o = 0; o < obs_meshes.size(); o++) {
			compute_ws_data(*obs_meshes[o]);
			update_x0(*obs_meshes[o]);
		}
		for (int z = 0; z < zones.size(); z++)
			delete zones[z];
		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);
	}


// Bridson-Harmon Collision Response
	bool collision_response_bh(vector<Mesh *> &meshes, const vector<Constraint *> &cons,
							   const vector<Mesh *> &obs_meshes, double dt) {

		bool resolvedAnything = false;

		// Old positions of all the nodes in the cloth and obstacle meshes
		arcsim::xold = node_positions(meshes);
		arcsim::xold_obs = node_positions(obs_meshes);

		// Cloth and obstacle meshes
		arcsim::meshes = &meshes;
		arcsim::obs_meshes = &obs_meshes;

		std::cout << " -- Building accel structs\n";

		// Acceleration structures for cloth and obstacle meshes
		vector<AccelStruct *> accs = create_accel_structs(meshes, true);
		vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, true);

		// 3. Detect proximity impacts using start of timestep positions
		std::cout << " -- Finding prox impacts \n";
		std::vector<Impact> proxImpacts = find_proximity_impacts(accs, obs_accs);

		std::cout << " -- Pruning prox impacts \n";
		prune_impacts(proxImpacts, impacts_equal2);

		if (proxImpacts.size() > 0) {
			resolvedAnything = true;
		}

		// 4. Apply repulsion impulses
		std::cout << " -- Applying repulsion to " << proxImpacts.size() << " prox impacts\n";
		//apply_repulsion(proxImpacts,meshes[0],dt);

		// 5. Apply friction
		std::cout << " -- Applying friction \n";
		if (arcsim::magic.friction_coeff > 1e-6) {
			apply_friction(proxImpacts);
		}

		// Compute midstep positions (x0 + vMid*dt)
		//std::cout << " -- Computing midstep positions \n";
		//for(int m = 0; m < meshes.size(); m++){
		//	Mesh* mesh = meshes[m];
		//	for(int n = 0; n < mesh->nodes.size(); n++){
		//		mesh->nodes[n]->x = mesh->nodes[n]->x0 + mesh->nodes[n]->vMid * dt;
		//	}
		//}

		std::cout << " -- Continuous Collision / Impact Zones \n";

		if (arcsim::magic.inelastic_projection) {

			// Impact Zones
			std::vector<ImpactZone *> zones;

			// Initializing mass of obstacles to large value. Can decrease it when attempting to solve collision problem
			arcsim::obs_mass = 1e3;

			// Keep track of number of iterations in the inner loop
			int iter;

			// First pass, don't allow obstacles to deform. Second pass, let them deform.
			for (int deform = 0; deform <= 1; deform++) {
				arcsim::deform_obstacles = deform;
				zones.clear();

				std::cout << " -- deform: " << deform << std::endl;

				// Inner loop - iterate until no more impacts are left to be resolved or until max iters is reached
				for (iter = 0; iter < max_iter; iter++) {

					std::cout << " -- ~~~ iter : " << iter << std::endl;

					// Update the active zones, if there are any
					if (!zones.empty()) {
						update_active(accs, obs_accs, zones);
					}

					// Continuous collision detection using midstep velocities
					std::cout << " -- ~~~ Finding continuous impacts\n";
					vector<Impact> continuousImpacts = find_continuous_impacts(accs, obs_accs);
					std::cout << " -- ~~~ Getting independent impacts\n";
					continuousImpacts = independent_impacts(continuousImpacts);

					if (continuousImpacts.size() > 0) {
						resolvedAnything = true;
					}

					// If no collisions were found, then nothing left to resolve
					if (continuousImpacts.empty()) {
						break;
					}

					// Add the found impacts to the zones
					std::cout << " -- ~~~ Adding " << continuousImpacts.size() << " impacts to zones\n";
					add_impacts(continuousImpacts, zones);

					// Apply inelastic projection on each zone
					std::cout << "-- ~~~ Applying inelastic projection on " << zones.size() << " zones\n";
					for (int z = 0; z < zones.size(); z++) {
						ImpactZone *zone = zones[z];
						//reset_zone_velocities(zone);
						apply_inelastic_projection(zone, cons);
						//apply_harmon_friction(zone, cons);
					}

					// Update the acceleration structures
					std::cout << "-- ~~~ Updating acceleration structures\n";
					for (int a = 0; a < accs.size(); a++) {
						update_accel_struct(*accs[a]);
					}
					for (int a = 0; a < obs_accs.size(); a++) {
						update_accel_struct(*obs_accs[a]);
					}

					// Halve the obstacle mass to allow for more deformations in next attempt
					if (deform_obstacles) {
						arcsim::obs_mass /= 2;
					}
				}

			}

			if (iter == max_iter) {
				cerr << "Collision resolution failed to converge!" << endl;
				debug_save_meshes(meshes, "meshes");
				debug_save_meshes(obs_meshes, "obsmeshes");
				exit(1);
			}

			// Delete all the zones and destroy acceleration structures
			for (int z = 0; z < zones.size(); z++)
				delete zones[z];

		}

		std::cout << " -- Computing ws data and updating cloth and obs meshes \n";





		// Recompute world space data and buffered positions in the cloth and obstacle meshes
		//for (int m = 0; m < meshes.size(); m++) {
		//    compute_ws_data(*meshes[m]);
		//    update_x0(*meshes[m]);
		//}
		//for (int o = 0; o < obs_meshes.size(); o++) {
		//    compute_ws_data(*obs_meshes[o]);
		//    update_x0(*obs_meshes[o]);
		//}


		std::cout << " -- Destroying accel structs \n";

		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);

		return resolvedAnything;

	}

// Independent impacts

	bool operator<(const Impact &impact0, const Impact &impact1) {
		return impact0.t < impact1.t;
	}

	bool conflict(const Impact &impact0, const Impact &impact1);

	vector<Impact> independent_impacts(const vector<Impact> &impacts) {
		vector<Impact> sorted = impacts;
		sort(sorted.begin(), sorted.end());
		vector<Impact> indep;
		for (int e = 0; e < sorted.size(); e++) {
			const Impact &impact = sorted[e];
			bool con = false;
			for (int e1 = 0; e1 < indep.size(); e1++)
				if (conflict(impact, indep[e1]))
					con = true;
			if (!con)
				indep.push_back(impact);
		}
		return indep;
	}

	bool conflict(const Impact &i0, const Impact &i1) {
		return (is_free(i0.nodes[0]) && is_in(i0.nodes[0], i1.nodes, 4))
			   || (is_free(i0.nodes[1]) && is_in(i0.nodes[1], i1.nodes, 4))
			   || (is_free(i0.nodes[2]) && is_in(i0.nodes[2], i1.nodes, 4))
			   || (is_free(i0.nodes[3]) && is_in(i0.nodes[3], i1.nodes, 4));
	}

	void update_active(const vector<AccelStruct *> &accs,
					   const vector<AccelStruct *> &obs_accs,
					   const vector<ImpactZone *> &zones) {
		for (int a = 0; a < accs.size(); a++)
			mark_all_inactive(*accs[a]);
		for (int a = 0; a < obs_accs.size(); a++)
			mark_all_inactive(*obs_accs[a]);
		for (int z = 0; z < zones.size(); z++) {
			const ImpactZone *zone = zones[z];
			if (!zone->active)
				continue;
			for (int n = 0; n < zone->nodes.size(); n++) {
				const Node *node = zone->nodes[n];
				pair<bool, int> mi = find_in_meshes(node);
				AccelStruct *acc = (mi.first ? accs : obs_accs)[mi.second];
				for (int v = 0; v < node->verts.size(); v++)
					for (int f = 0; f < node->verts[v]->adjf.size(); f++)
						mark_active(*acc, node->verts[v]->adjf[f]);
			}
		}
	}



// Impact zones

	ImpactZone *find_or_create_zone(const Node *node, vector<ImpactZone *> &zones);

	void merge_zones(ImpactZone *zone0, ImpactZone *zone1,
					 vector<ImpactZone *> &zones);

	void add_impacts(const vector<Impact> &impacts, vector<ImpactZone *> &zones) {
		for (int z = 0; z < zones.size(); z++)
			zones[z]->active = false;
		for (int i = 0; i < impacts.size(); i++) {
			const Impact &impact = impacts[i];
			Node *node = impact.nodes[is_free(impact.nodes[0]) ? 0 : 3];
			ImpactZone *zone = find_or_create_zone(node, zones);
			for (int n = 0; n < 4; n++)
				if (is_free(impact.nodes[n]) || arcsim::deform_obstacles)
					merge_zones(zone, find_or_create_zone(impact.nodes[n], zones),
								zones);
			zone->impacts.push_back(impact);
			zone->active = true;
		}
	}

	ImpactZone *find_or_create_zone(const Node *node, vector<ImpactZone *> &zones) {
		for (int z = 0; z < zones.size(); z++)
			if (is_in((Node *) node, zones[z]->nodes))
				return zones[z];
		ImpactZone *zone = new ImpactZone;
		zone->nodes.push_back((Node *) node);
		zones.push_back(zone);
		return zone;
	}

	void merge_zones(ImpactZone *zone0, ImpactZone *zone1,
					 vector<ImpactZone *> &zones) {
		if (zone0 == zone1)
			return;
		append(zone0->nodes, zone1->nodes);
		append(zone0->impacts, zone1->impacts);
		exclude(zone1, zones);
		delete zone1;
	}







// Response

	struct NormalOpt : public NLConOpt {
		ImpactZone *zone;
		double inv_m;

		NormalOpt() : zone(NULL), inv_m(0) { nvar = ncon = 0; }

		NormalOpt(ImpactZone *zone) : zone(zone), inv_m(0) {
			nvar = zone->nodes.size() * 3;
			ncon = zone->impacts.size();
			for (int n = 0; n < zone->nodes.size(); n++)
				inv_m += 1 / get_mass(zone->nodes[n]);
			inv_m /= zone->nodes.size();
		}

		void initialize(double *x) const;

		void precompute(const double *x) const;

		double objective(const double *x) const;

		void obj_grad(const double *x, double *grad) const;

		double constraint(const double *x, int i, int &sign) const;

		void con_grad(const double *x, int i, double factor, double *grad) const;

		void finalize(const double *x) const;
	};


	void reset_zone_velocities(ImpactZone *zone) {

		for (int i = 0; i < zone->nodes.size(); i++) {
			Node *node = zone->nodes[i];
			//node->vPrimeMid = node->vTildeMid;
		}

	}


	void apply_inelastic_projection(ImpactZone *zone,
									const vector<Constraint *> &cons) {
		if (!zone->active)
			return;

		if (arcsim::magic.harmon_inelastic_projection) {

			for (int i = 0; i < zone->impacts.size(); i++) {
				for (int j = 0; j < 4; j++) {
					include(zone->impacts[i].nodes[j], zone->nodes);
				}
			}

			const int k = zone->impacts.size();
			const int n = zone->nodes.size();


			/*std::cout << "PREPROCESSING.\n";
            std::cout << "Printing the indices of all nodes in the impact zone.\n";
            for(int i = 0; i < zone->nodes.size(); i++){
                std::cout << zone->nodes[i]->index << " ";
            }
            std::cout << std::endl;

            std::cout << "Printing the indices of all nodes in all the impact signatures.\n";
            for(int i = 0; i < zone->impacts.size(); i++){
                for(int j = 0; j < 4; j++){
                    std::cout << zone->impacts[i].nodes[j]->index << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;*/

			// Constraint gradient matrix for entire impact zone
			Eigen::MatrixXd gradC(k, 3 * n);

			// Looping through each impact in the impact zone
			// Computing the gradient of the constraint, then
			// adding it to the total gradC matrix
			//std::cout << "adding local gradCi to global gradC\n";
			for (int i = 0; i < k; i++) {

				//std::cout << "iter: " << i << std::endl;

				Impact &impact = zone->impacts[i];
				std::vector<Vec3> gradCi(4, Vec3(0));
				Vec3 n = impact.n;

				if (impact.type == Impact::VF) {
					//std::cout << "VF impact\n";
					for (int j = 0; j < 4; j++) {
						gradCi[j] = impact.w[j] * n;

					}
				} else if (impact.type == Impact::EE) {
					//std::cout << "EE impact\n";
					for (int j = 0; j < 4; j++) {
						gradCi[j] = -impact.w[j] * n;
					}
				} else {
					std::cout << "Error: Impact does not have proper type in inelastic projection.\n";
					exit(0);
				}

				//std::cout << "Done computing the local 4. Going to loop up and add to global\n";

				//std::cout << "zone->nodes size is " << zone->nodes.size() << std::endl;

				for (int j = 0; j < 4; j++) {
					int index = find(impact.nodes[j], zone->nodes);
					if (index >= 0) {
						gradC(i, 3 * index) = gradCi[j][0];
						gradC(i, 3 * index + 1) = gradCi[j][1];
						gradC(i, 3 * index + 2) = gradCi[j][2];
					}
				}
			}


			//std::cout << "gradC: " << gradC << std::endl;

			// Setting the inverse mass matrix
			//std::cout << "Filling inv mass matrix ...\n";
			Eigen::MatrixXd Minv(3 * n, 3 * n);
			for (int i = 0; i < n; i++) {
				Minv(3 * i, 3 * i) = 1.0 / zone->nodes[i]->m;
				Minv(3 * i + 1, 3 * i + 1) = 1.0 / zone->nodes[i]->m;
				Minv(3 * i + 2, 3 * i + 2) = 1.0 / zone->nodes[i]->m;
			}

			//std::cout << "Minv: " << Minv << std::endl;

			// Setting configurational velocity vector
			//std::cout << "Filling qdot ... \n";
			Eigen::VectorXd qdot(3 * n);
			for (int i = 0; i < n; i++) {
				Vec3 vel = zone->nodes[i]->vTildeMid;
				qdot(3 * i) = vel[0];
				qdot(3 * i + 1) = vel[1];
				qdot(3 * i + 2) = vel[2];
			}

			// Setting up and solving the linear system for the lagrange multipliers
			Eigen::MatrixXd A = gradC * Minv * gradC.transpose();
			Eigen::VectorXd b = gradC * qdot;
			Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
			Eigen::VectorXd lambda = solver.compute(A.sparseView()).solve(b);

			// Computing the change in momentum for inelastic projection
			Eigen::VectorXd dp = -gradC.transpose() * lambda;

			// Applying the momentum change to the nodes in the zone
			for (int i = 0; i < n; i++) {
				Node *node = zone->nodes[i];
				Vec3 dpNode(dp[3 * i], dp[3 * i + 1], dp[3 * i + 2]);
				node->vPrimeMid += (dpNode / node->m);
				node->vFinalMid = node->vPrimeMid;
			}


		} else {

			// Old inelastic projection code

			augmented_lagrangian_method(NormalOpt(zone));

		}

	}


	void apply_harmon_friction(ImpactZone *zone, const vector<Constraint *> &cons) {





		// Applying velocity change to the nodes in the zone
		//for(int i = 0; i < n; i++){
		//	Node* node = zone->nodes[i];
		//	Vec3 dqdotNode( dqdot[3*i] , dqdot[3*i+1] , dqdot[3*i+2] );
		//	node->vFinalMid = node->vTildeMid + dqdotNode;
		//}

	}


	void NormalOpt::initialize(double *x) const {
		for (int n = 0; n < zone->nodes.size(); n++)
			set_subvec(x, n, zone->nodes[n]->x);
	}

	void NormalOpt::precompute(const double *x) const {
		for (int n = 0; n < zone->nodes.size(); n++)
			zone->nodes[n]->x = get_subvec(x, n);
	}

	const Vec3 &get_xold(const Node *node);

	double NormalOpt::objective(const double *x) const {
		double e = 0;
		for (int n = 0; n < zone->nodes.size(); n++) {
			const Node *node = zone->nodes[n];
			Vec3 dx = node->x - get_xold(node);
			e += inv_m * get_mass(node) * norm2(dx) / 2;
		}
		return e;
	}

	void NormalOpt::obj_grad(const double *x, double *grad) const {
		for (int n = 0; n < zone->nodes.size(); n++) {
			const Node *node = zone->nodes[n];
			Vec3 dx = node->x - get_xold(node);
			set_subvec(grad, n, inv_m * get_mass(node) * dx);
		}
	}

	double NormalOpt::constraint(const double *x, int j, int &sign) const {
		sign = 1;
		//double c = -::thickness;
		double c = -arcsim::magic.projection_thickness;
		const Impact &impact = zone->impacts[j];
		for (int n = 0; n < 4; n++)
			c += impact.w[n] * dot(impact.n, impact.nodes[n]->x);
		return c;
	}

	void NormalOpt::con_grad(const double *x, int j, double factor,
							 double *grad) const {
		const Impact &impact = zone->impacts[j];
		for (int n = 0; n < 4; n++) {
			int i = find(impact.nodes[n], zone->nodes);
			if (i != -1)
				add_subvec(grad, i, factor * impact.w[n] * impact.n);
		}
	}

	void NormalOpt::finalize(const double *x) const {
		precompute(x);
	}

	const Vec3 &get_xold(const Node *node) {
		pair<bool, int> mi = find_in_meshes(node);
		int ni = get_index(node, mi.first ? *arcsim::meshes : *arcsim::obs_meshes);
		return (mi.first ? arcsim::xold : arcsim::xold_obs)[ni];
	}

}
