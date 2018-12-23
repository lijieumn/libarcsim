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

#include "simulation.hpp"

#include "bhcollision.hpp"
#include "collision.hpp"
#include "dynamicremesh.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "nearobs.hpp"
#include "physics.hpp"
#include "plasticity.hpp"
#include "popfilter.hpp"
#include "proximity.hpp"
#include "separate.hpp"
#include "strainlimiting.hpp"
#include "io.hpp"
#include "mergeimpacts.hpp"
#include "mergehelper.hpp"
#include "visualdebug.hpp"
#include "log.hpp"

#include "sparse_solver.hpp"
#include <iostream>
#include <fstream>

using namespace std;

namespace arcsim {

    static const bool verbose = false;
    static const int proximity = Simulation::Proximity,
            physics = Simulation::Physics,
            strainlimiting = Simulation::StrainLimiting,
            collision = Simulation::Collision,
            remeshing = Simulation::Remeshing,
            separation = Simulation::Separation,
            popfilter = Simulation::PopFilter,
            plasticity = Simulation::Plasticity;

    Simulation::~Simulation() {
        delete_simulation(*this);
    }

    void plasticity_step(Simulation &sim);

    void strainlimiting_step(Simulation &sim, const vector<Constraint *> &cons);

    void strainzeroing_step(Simulation &sim);

    void equilibration_step(Simulation &sim);

    void collision_step_bh(Simulation &sim);

    void remeshing_step(Simulation &sim, bool initializing = false);

    void validate_handles(const Simulation &sim);

    void prepare(Simulation &sim) {

#ifndef SILENCE_ARGUS
        std::cout << "Preparing Simulation...\n";
#endif

        sim.bufferedMeshes.resize(Simulation::nStages);
        for (int i = 0; i < Simulation::nStages; i++) {
            if (sim.bufferedMeshes[i]) {
                delete_mesh(*sim.bufferedMeshes[i]);
                delete sim.bufferedMeshes[i];
            }
            sim.bufferedMeshes[i] = new Mesh();
        }

        sim.cloth_meshes.resize(sim.cloths.size());
        for (int c = 0; c < sim.cloths.size(); c++) {
            compute_masses(sim.cloths[c]);
            sim.cloth_meshes[c] = &sim.cloths[c].mesh;
            update_x0(*sim.cloth_meshes[c]);
            update_n0(*sim.cloth_meshes[c]);
        }
        sim.obstacle_meshes.resize(sim.obstacles.size());
        for (int o = 0; o < sim.obstacles.size(); o++) {
            sim.obstacle_meshes[o] = &sim.obstacles[o].get_mesh();
            update_x0(*sim.obstacle_meshes[o]);
            update_n0(*sim.obstacle_meshes[o]);

            // setting all the obstacle node masses to inf
            for (int j = 0; j < sim.obstacle_meshes[o]->nodes.size(); j++) {
                sim.obstacle_meshes[o]->nodes[j]->m = 1.e8;
            }

        }

        sim.liveMesh = sim.cloth_meshes[0];
        sim.bufferId = -1;

    }

    void relax_initial_state(Simulation &sim) {
        sim.cloths[0].mesh.reset_face_size_min(sim.cloths[0].remeshing.size_min);

        if (arcsim::magic.relax_initial_state) {

#ifndef SILENCE_ARGUS
            std::cout << "Relaxing initial state...\n";
#endif

            validate_handles(sim);
            if (arcsim::magic.preserve_creases)
                for (int c = 0; c < sim.cloths.size(); c++)
                    reset_plasticity(sim.cloths[c]);
            bool equilibrate = false;
            if (equilibrate) {
#ifndef SILENCE_ARGUS
                std::cout << "-- Equilibration step\n";
#endif
                equilibration_step(sim);
#ifndef SILENCE_ARGUS
                std::cout << "-- Remeshing step\n";
#endif
                remeshing_step(sim, true);
#ifndef SILENCE_ARGUS
                std::cout << "-- Equilibration step\n";
#endif
                equilibration_step(sim);
            } else {
                remeshing_step(sim, true);
                strainzeroing_step(sim);
                remeshing_step(sim, true);
                strainzeroing_step(sim);
            }
            if (arcsim::magic.preserve_creases)
                for (int c = 0; c < sim.cloths.size(); c++)
                    reset_plasticity(sim.cloths[c]);
            arcsim::magic.preserve_creases = false;
            if (arcsim::magic.fixed_high_res_mesh)
                sim.enabled[remeshing] = false;

        }

        // std::cout << "before creating kd tree, sim.cloth_meshes size is " << sim.cloth_meshes.size() << std::endl;
        // std::cout << "also, the first mesh has num verts size of " << sim.cloth_meshes[0]->verts.size() << std::endl;

        for (int m = 0; m < sim.cloth_meshes.size(); m++) {

            // std::cout << "m = ... " << m << std::endl;

            Mesh *mesh = sim.cloth_meshes[m];
            mesh->kdTree = new KDTree(mesh->verts);

            for (int n = 0; n < mesh->nodes.size(); n++) {
                mesh->nodes[n]->inMesh = true;
            }

            // if (::magic.merge_radius < sim.cloths[m].remeshing.size_min) {
            // 	sim.cloths[m].remeshing.size_min = ::magic.merge_radius;
            // } else {
            // 	::magic.merge_radius = sim.cloths[m].remeshing.size_min;
            // }
        }
        for (int m = 0; m < sim.obstacle_meshes.size(); m++) {
            Mesh *mesh = sim.obstacle_meshes[m];
            for (int n = 0; n < mesh->nodes.size(); n++) {
                mesh->nodes[n]->inMesh = false;
            }
        }

    }

    void validate_handles(const Simulation &sim) {
        for (int h = 0; h < sim.handles.size(); h++) {
            vector<Node *> nodes = sim.handles[h]->get_nodes();
            for (int n = 0; n < nodes.size(); n++) {
                if (!nodes[n]->preserve) {
                    cout << "Constrained node " << nodes[n]->index << " will not be preserved by remeshing" << endl;
                    abort();
                }
            }
        }
    }

    vector<Constraint *> get_constraints(Simulation &sim, bool include_proximity);

    void delete_constraints(const vector<Constraint *> &cons);

    void update_obstacles(Simulation &sim, bool update_positions = true);

    void delete_simulation(Simulation &sim) {
        for (int m = 0; m < sim.cloth_meshes.size(); m++) {
            Mesh *mesh = sim.cloth_meshes[m];
            delete_mesh(*mesh);
            // delete mesh; // don't do this since cloth_meshes was not created by new operator.
        }
        for (int o = 0; o < sim.obstacles.size(); o++) {
            Obstacle &obs = sim.obstacles[o];
            delete_mesh(obs.curr_state_mesh);
            delete_mesh(obs.base_mesh);
            if (obs.prev) {
                delete_mesh(*obs.prev);
                delete obs.prev;
            }
            if (obs.prevprev) {
                delete_mesh(*obs.prevprev);
                delete obs.prevprev;
            }
            if (obs.next) {
                delete_mesh(*obs.next);
                delete obs.next;
            }
            if (obs.nextnext) {
                delete_mesh(*obs.nextnext);
                delete obs.nextnext;
            }
        }
        // for (int m = 0; m < sim.obstacle_meshes.size(); m++) {
        // 	Mesh* mesh = sim.obstacle_meshes[m];
        // 	delete_mesh(*mesh);
        // 	// delete mesh;	// don't do this since obstacle_meshes was not created by new operator.
        // }
        for (int m = 0; m < sim.bufferedMeshes.size(); m++) {
            Mesh *mesh = sim.bufferedMeshes[m];
            delete_mesh(*mesh);
            delete mesh;
        }
        // delete_mesh(*sim.liveMesh);
        // delete sim.liveMesh;
        for (int c = 0; c < sim.cloths.size(); c++) {
            Cloth &cloth = sim.cloths[c];
            for (int m = 0; m < cloth.materials.size(); m++) {
                delete cloth.materials[m];
            }
        }

    }

    vector<Constraint *> get_constraints(Simulation &sim, bool include_proximity) {
        vector<Constraint *> cons;
        for (int h = 0; h < sim.handles.size(); h++)
            append(cons, sim.handles[h]->get_constraints(sim.time));
        for (int m = 0; m < sim.cloth_meshes.size(); m++) {
            mark_nodes_to_preserve(*sim.cloth_meshes[m]);
        }
        if (include_proximity && sim.enabled[proximity]) {
            sim.timers[proximity].tick();
            append(cons, proximity_constraints(sim.cloth_meshes,
                                               sim.obstacle_meshes,
                                               sim.friction, sim.obs_friction));
            sim.timers[proximity].tock();
        }
        return cons;
    }

    void delete_constraints(const vector<Constraint *> &cons) {
        for (int c = 0; c < cons.size(); c++)
            delete cons[c];
    }

// Steps

    void update_velocities(vector<Mesh *> &meshes, vector<Vec3> &xold, double dt);

    void compute_avg_velocities(std::vector<Mesh *> &meshes, double dt);

    void step_mesh(Mesh &mesh, double dt);


    void add_ext_and_morphs(Simulation &sim, const vector<Constraint *> &cons, std::vector<Vec3> &fext,
                            std::vector<Mat3x3> &Jext, int c) {

        int nn = sim.cloths[c].mesh.nodes.size();
        vector<Vec3> fext2(nn, Vec3(0));
        vector<Mat3x3> Jext2(nn, Mat3x3(0));
        fext = fext2;
        Jext = Jext2;
        add_external_forces(sim.cloths[c], sim.gravity, sim.wind, fext, Jext);
        for (int m = 0; m < sim.morphs.size(); m++) {
            if (sim.morphs[m].mesh == &sim.cloths[c].mesh) {
                add_morph_forces(sim.cloths[c], sim.morphs[m], sim.time, sim.step_time, fext, Jext);
            }
        }


    }

    enum CompareType {
        DIFFERENT,
        EQUAL,
        EARLIER,
        LATER
    };

    typedef CompareType (*CompareCallback)(const Impact &imp1, const Impact &imp2);


    CompareType impacts_equal(const Impact &imp1, const Impact &imp2) {
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

    CompareType vf_vertex_equal(const Impact &imp1, const Impact &imp2) {
        if (imp1.type == Impact::VF && imp2.type == Impact::VF) {
            if (imp1.nodes[0] == imp2.nodes[0]) {
                return imp1.t < imp2.t ? EARLIER : LATER;
            }
        }
        return DIFFERENT;
    }


    bool find_in_mesh3(const Node *node, const Mesh &mesh) {
        for (int n = 0; n < mesh.nodes.size(); n++) {
            if (node == mesh.nodes[n]) {
                return true;
            }
        }
        return false;
    }

    struct ImpactCompare {
        bool operator()(const Impact &lhs, const Impact &rhs) const {
            if (lhs.type != rhs.type) {
                return lhs.type < rhs.type;
            }
            for (int n = 0; n < 4; n++) {
                if (lhs.nodes[n] != rhs.nodes[n]) {
                    return lhs.nodes[n] < rhs.nodes[n];
                }
            }
            return false;
        }
    };

    vector<Impact> faster_pruning(vector<Impact> &impacts) {
        map<Impact, int, ImpactCompare> iMap;
        for (int i = 0; i < impacts.size(); i++) {
            Impact &impact = impacts[i];
            auto search = iMap.find(impact);
            if (search != iMap.end()) {    // found
                // const Impact& impact2 = search->first;
                // if (impact.d )
            } else {
                iMap[impact] = i;
            }
        }
        vector<Impact> pruned;
        for (auto &search : iMap) {
            pruned.push_back(impacts[search.second]);
        }
        return pruned;
    }

/* Removes duplicate impacts */
    void prune_impacts(Simulation &sim, std::vector<std::vector<Impact> > &impacts, CompareCallback compare) {

        std::vector<int> markToDelete;
        vector<int> markToRemain;

        for (int i = 0; i < impacts[0].size(); i++) {
            Impact &imp1 = impacts[0][i];
            CompareType type = DIFFERENT;
            for (int j = 0; j < markToRemain.size(); j++) {
                Impact &imp2 = impacts[0][markToRemain[j]];
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
        for (int i = 0; i < impacts[0].size(); i++) {
            if (markToRemain.size() <= i - delete_count || markToRemain[i - delete_count] != i) {
                delete_count++;
                markToDelete.push_back(i);
            }
        }

        for (int i = int(markToDelete.size()) - 1; i >= 0; i--) {
            impacts[0].erase(impacts[0].begin() + markToDelete[i]);
        }

        markToDelete.clear();

    }


/* VF and EE refinement at impact sites to generate ArgusImpacts */
    void produce_argus_impacts(Simulation &sim, const std::vector<std::vector<Impact> > &impacts,
                               std::vector<std::vector<ArgusImpact> > &argusImpacts) {

        // copy old meshes
        vector<Mesh> old_meshes(sim.cloths.size());
        vector<Mesh *> old_meshes_p(sim.cloths.size()); // for symmetry in separate()
        for (int c = 0; c < sim.cloths.size(); c++) {
            old_meshes[c] = deep_copy(sim.cloths[c].mesh);
            old_meshes_p[c] = &old_meshes[c];
        }

        // init the meshes
        init_meshes(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes);

        // remesh
        for (int c = 0; c < sim.cloths.size(); c++) {
            vector<Plane> planes = nearest_obstacle_planes(sim.cloths[c].mesh,
                                                           sim.obstacle_meshes);

            std::vector<ArgusImpact> aImps = collision_refine(sim.cloths[c], sim.obstacle_meshes, planes, impacts[c],
                                                              true);
            argusImpacts.push_back(aImps);

        }

        // delete old meshes
        for (int c = 0; c < sim.cloths.size(); c++)
            delete_mesh(old_meshes[c]);

    }


    void coarsen_impacts(Simulation &sim, std::vector<std::vector<Impact> > &impacts) {

        // copy old meshes
        vector<Mesh> old_meshes(sim.cloths.size());
        vector<Mesh *> old_meshes_p(sim.cloths.size()); // for symmetry in separate()
        for (int c = 0; c < sim.cloths.size(); c++) {
            old_meshes[c] = deep_copy(sim.cloths[c].mesh);
            old_meshes_p[c] = &old_meshes[c];
        }
        // remesh
        for (int c = 0; c < sim.cloths.size(); c++) {
            vector<Plane> planes = nearest_obstacle_planes(sim.cloths[c].mesh,
                                                           sim.obstacle_meshes);

            init_meshes(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes);
            collision_coarsen(sim.cloths[c], sim.obstacle_meshes, planes, true);

        }
        // delete old meshes
        for (int c = 0; c < sim.cloths.size(); c++)
            delete_mesh(old_meshes[c]);

    }

    void prune_argus_impacts(Simulation &sim, std::vector<std::vector<ArgusImpact> > &argusImpacts) {

        for (int i = 0; i < sim.cloths[0].mesh.nodes.size(); i++) {

            Node *node = sim.cloths[0].mesh.nodes[i];
            std::vector<int> matches;
            double dmin = 1.e8;
            int jmin = 10000000;

            for (int j = 0; j < argusImpacts[0].size(); j++) {

                if (argusImpacts[0][j].nodeA == node) {
                    double d = norm(argusImpacts[0][j].nodeA->x - argusImpacts[0][j].posB);
                    if (d < dmin) {
                        dmin = d;
                        jmin = j;
                    }
                    matches.push_back(j);
                }

            }

            if (matches.size() > 0) {
                for (int k = matches.size() - 1; k >= 0; k--) {
                    if (matches[k] != jmin) {
                        argusImpacts[0].erase(argusImpacts[0].begin() + matches[k]);
                    }
                }
            }


        }
    }


    void setPreMergeNormals(Simulation &sim, std::vector<std::vector<Impact> > &impacts) {

        VisualDebugger *vd = VisualDebugger::getInstance();
        for (int i = 0; i < impacts[0].size(); i++) {
            Impact impact = impacts[0][i];
            Vec3 p;
            Vec3 p2;
            if (impact.type == Impact::EE) {
                p = impact.nodes[2]->x * (-impact.w[2]) + impact.nodes[3]->x * (-impact.w[3]);
                p2 = impact.nodes[0]->x * impact.w[0] + impact.nodes[1]->x * impact.w[1];
            } else if (impact.type == Impact::VF) {
                p = impact.nodes[1]->x * (-impact.w[1]) + impact.nodes[2]->x * (-impact.w[2]) +
                    impact.nodes[3]->x * (-impact.w[3]);
                p2 = impact.nodes[0]->x;
            } else if (impact.type == Impact::VE) {
                p = impact.nodes[1]->x * (impact.w[1]) + impact.nodes[2]->x * (impact.w[2]);
                p2 = impact.nodes[0]->x;
            } else {
                p = impact.nodes[1]->x;
                p2 = impact.nodes[0]->x;
            }
            vd->addVisualPoint3(p, Vec3(1, 0, 0), 'b');
            vd->addVisualLine3(p, p + impact.n * .05, Vec3(1, 0, 0), 'b');
            vd->addVisualPoint3(p2, Vec3(1, 0, 0), 'b');

        }

    }


    void setPostMergeNormals(Simulation &sim, std::vector<std::vector<Impact> > &impacts) {

        VisualDebugger *vd = VisualDebugger::getInstance();
        for (int i = 0; i < impacts[0].size(); i++) {
            Impact impact = impacts[0][i];
            Vec3 p;
            Vec3 p2;
            if (impact.type == Impact::EE) {
                p = impact.nodes[2]->x * (-impact.w[2]) + impact.nodes[3]->x * (-impact.w[3]);
                p2 = impact.nodes[0]->x * impact.w[0] + impact.nodes[1]->x * impact.w[1];
            } else if (impact.type == Impact::VF) {
                p = impact.nodes[1]->x * (-impact.w[1]) + impact.nodes[2]->x * (-impact.w[2]) +
                    impact.nodes[3]->x * (-impact.w[3]);
                p2 = impact.nodes[0]->x;
            } else if (impact.type == Impact::VE) {
                p = impact.nodes[1]->x * (impact.w[1]) + impact.nodes[2]->x * (impact.w[2]);
            } else {
                p = impact.nodes[1]->x;
            }
            vd->addVisualPoint3(p, Vec3(0, 0, 1), 'a');
            vd->addVisualLine3(p, p + impact.n * .05, Vec3(0, 0, 1), 'a');

            char switchChar;
            Vec3 color;

            vd->addVisualPoint3(p2, color, switchChar);
            vd->addVisualPoint3(p, color, switchChar);
            //vd->addVisualLine3(p2,p,color,switchChar);
            vd->addVisualLine3(p2, p2 + 0.05 * impact.n, color, switchChar);
            //vd->addVisualLine3(p,p-0.05*impact.n, color, switchChar);



        }
    }


    void projectArgusImpacts(Simulation &sim, std::vector<std::vector<ArgusImpact> > &argusImpacts) {

        VisualDebugger *vd = VisualDebugger::getInstance();

        // Project impacts before passing into the solver
        for (int i = 0; i < argusImpacts[0].size(); i++) {
            ArgusImpact &argusImpact = argusImpacts[0][i];
            Vec3 posA = argusImpact.nodeA->x;
            Vec3 posB = argusImpact.nodeB ? argusImpact.nodeB->x : argusImpact.posB;
            Vec3 normal = argusImpact.normal;
            // double d = dot(posA - posB, normal);
            // double thickness = ::magic.projection_thickness;
            // if (d < thickness) {
            // 	double difference = thickness - d;
            // 	argusImpact.nodeA->x += difference*normal;
            // }

            vd->addVisualPoint3(posA, Vec3(0, 1, 1), 'n');
            vd->addVisualPoint3(posB, Vec3(0, 1, 1), 'n');
            vd->addVisualLine3(posA, posA + normal * .05, Vec3(0, 1, 1), 'n');

        }
    }

    void projectImpacts(vector<Impact> &impacts) {
        for (int i = 0; i < impacts.size(); i++) {
            Impact &impact = impacts[i];
            if (impact.self) {    // Not clear how to project self contact
                continue;
            }
            if (impact.type == Impact::VF) {
                if (!impact.inverted) {
                    Node *n = impact.nodes[0];
                    Vec3 obs_x = impact.nodes[1]->x * impact.w[1] + impact.nodes[2]->x * impact.w[2] +
                                 impact.nodes[3]->x * impact.w[3];
                    Vec3 normal = impact.n;
                    double d = dot(n->x + obs_x, normal);
                    double thickness = arcsim::magic.projection_thickness;
                    if (d < thickness) {
                        double difference = thickness - d;
                        n->x += difference * normal;
                    }
                }
            }
        }
    }

    void cacheArgusImpacts(Simulation &sim, std::vector<std::vector<ArgusImpact> > &argusImpacts) {

        sim.cloths[0].mesh.impactsCache.clear();
        for (int i = 0; i < argusImpacts[0].size(); i++) {
            Mesh::ImpactCache impactCache;
            impactCache.pos = argusImpacts[0][i].nodeA->x;
            impactCache.pos2 = argusImpacts[0][i].posB;
            impactCache.normal = argusImpacts[0][i].normal;
            if (argusImpacts[0][i].type == ArgusImpact::VF) {
                impactCache.type = Mesh::ImpactCache::VF;
            } else if (argusImpacts[0][i].type == ArgusImpact::EE) {
                impactCache.type = Mesh::ImpactCache::EE;
            } else if (argusImpacts[0][i].type == ArgusImpact::VE) {
                impactCache.type = Mesh::ImpactCache::VE;
            } else if (argusImpacts[0][i].type == ArgusImpact::VV) {
                impactCache.type = Mesh::ImpactCache::VV;
            }
            impactCache.inverted = argusImpacts[0][i].inverted;
            impactCache.self = argusImpacts[0][i].self;
            impactCache.debug = argusImpacts[0][i].debug;
            sim.cloths[0].mesh.impactsCache.push_back(impactCache);
        }

    }

    double smallest(Mesh &mesh) {
        double min = 100;
        Vec3 vi, vj;
        int di, dj;
        for (int i = 0; i < mesh.verts.size(); i++) {
            for (int j = i + 1; j < mesh.verts.size(); j++) {
                Vert *v1 = mesh.verts[i];
                Vert *v2 = mesh.verts[j];
                double dis = norm(v1->u - v2->u);
                if (dis < min) {
                    min = dis;
                    di = i;
                    dj = j;
                    vi = v1->node->x;
                    vj = v2->node->x;
                }
            }
        }
        // VisualDebugger *vd = VisualDebugger::getInstance();
        // vd->addVisualPoint3(mesh.verts[di]->node->x, Vec3(1, 0, 0), 'l');
        // vd->addVisualPoint3(mesh.verts[dj]->node->x, Vec3(1, 0, 0), 'l');
        return min;
    }

    double smallest_ve(Mesh &mesh) {
        double min = 100;
        for (int f = 0; f < mesh.faces.size(); f++) {
            Face *face = mesh.faces[f];
            for (int i = 0; i < 3; i++) {
                Vec2 proj;
                double d = material_ve_projection(face->v[i]->u, face->v[(i + 1) % 3]->u, face->v[(i + 2) % 3]->u,
                                                  proj);
                if (d >= 0 && d < min) {
                    min = d;
                }
            }
        }
        return min;
    }

    void print_mesh(Mesh &mesh) {
        for (int i = 0; i < mesh.nodes.size(); i++) {
            Node *node = mesh.nodes[i];
            cout << "x:" << node->x << endl;
            cout << "v:" << node->v << endl;
            cout << "u:" << node->verts[0]->u << endl;
            cout << "y:" << node->y << endl;
            cout << "temp:" << node->temp << endl;
            cout << "temp2:" << node->temp2 << endl;
            cout << "adje size:" << node->adje.size() << endl;
            for (int j = 0; j < node->adje.size(); j++) {
                Edge *e = node->adje[j];
                cout << "\t" << e->index << endl;
            }
            cout << "n:" << node->n << endl;
        }
        for (int i = 0; i < mesh.faces.size(); i++) {
            Face *face = mesh.faces[i];
            cout << "node index:";
            for (int i = 0; i < 3; i++) {
                cout << "\t" << face->v[i]->node->index;
            }
            cout << endl;
            for (int i = 0; i < 3; i++) {
                cout << "\t" << face->v[i]->index;
            }
            cout << endl;
            // cout << "adje:" << endl;
            // cout << "\t" << (face->adje[0] ? face->adje[0]->index : -1) << endl;
            // cout << "\t" << (face->adje[1] ? face->adje[1]->index : -1) << endl;
            // cout << "\t" << (face->adje[2] ? face->adje[2]->index : -1) << endl;
        }
        // for (int i = 0; i < mesh.edges.size(); i++) {
        // 	Edge* edge = mesh.edges[i];
        // 	cout << "adjf:" << endl;
        // 	cout << "\t" << (edge->adjf[0] ? edge->adjf[0]->index : -1) << endl;
        // 	cout << "\t" << (edge->adjf[1] ? edge->adjf[1]->index : -1) << endl;
        // }
        for (int i = 0; i < mesh.verts.size(); i++) {
            Vert *vert = mesh.verts[i];
            cout << "adjf size:" << vert->adjf.size() << endl;
            for (int j = 0; j < vert->adjf.size(); j++) {
                Face *f = vert->adjf[j];
                cout << "\t" << f->index << endl;
            }
        }
        cout << endl;
    }

    bool impact_based_remesh(const vector<Impact> &impacts, Simulation &sim) {
        bool *impactFaces = new bool[sim.cloths[0].mesh.faces.size()];
        for (int i = 0; i < sim.cloths[0].mesh.faces.size(); i++) {
            impactFaces[i] = false;
        }
        for (int i = 0; i < sim.cloths[0].mesh.verts.size(); i++) {
            sim.cloths[0].mesh.verts[i]->impactVert = false;
        }
        for (int i = 0; i < impacts.size(); i++) {
            const Impact &impact = impacts[i];
            Vec2 impact_u[2];
            if (impact.type == Impact::VF) {
                if (impact.self) {    // self
                    impact_u[0] = impact.verts[0]->u;
                    impact_u[1] = impact.verts[1]->u * (-impact.w[1]) + impact.verts[2]->u * (-impact.w[2])
                                  + impact.verts[3]->u * (-impact.w[3]);
                } else if (impact.inverted) {    // inverted
                    impact_u[0] = impact.verts[1]->u * (-impact.w[1]) + impact.verts[2]->u * (-impact.w[2])
                                  + impact.verts[3]->u * (-impact.w[3]);
                } else {
                    // Standard
                    impact_u[0] = impact.verts[0]->u;
                }
            } else if (impact.type == Impact::EE) {
                if (impact.self) {
                    impact_u[0] = impact.verts[0]->u * (impact.w[0]) + impact.verts[1]->u * (impact.w[1]);
                    impact_u[1] = impact.verts[2]->u * (-impact.w[2]) + impact.verts[3]->u * (-impact.w[3]);
                } else {
                    impact_u[0] = impact.verts[0]->u * (impact.w[0]) + impact.verts[1]->u * (impact.w[1]);
                }
            } else if (impact.type == Impact::VE) {
                if (impact.self) {
                    impact_u[0] = impact.verts[0]->u;
                    impact_u[1] = impact.verts[1]->u * (-impact.w[1]) + impact.verts[2]->u * (-impact.w[2]);
                } else if (impact.inverted) {
                    impact_u[0] = impact.verts[1]->u * (-impact.w[1]) + impact.verts[2]->u * (-impact.w[2]);
                } else {
                    impact_u[0] = impact.verts[0]->u;
                }
            } else if (impact.type == Impact::VV) {
                if (impact.self) {
                    impact_u[0] = impact.verts[0]->u;
                    impact_u[1] = impact.verts[1]->u;
                } else {
                    impact_u[0] = impact.verts[0]->u;
                }
            }
            Face *face1 = get_enclosing_face(sim.cloths[0].mesh, impact_u[0]);
            double min = _MAX;
            int index = -1;
            for (int v = 0; v < 3; v++) {
                double d = norm(face1->v[v]->u - impact_u[0]);
                if (d < min) {
                    min = d;
                    index = v;
                }
            }
            if (min > arcsim::magic.merge_radius * sqrt(3)) {
                impactFaces[face1->index] = true;
                // face1->v[0]->impactVert = true;
                // face1->v[1]->impactVert = true;
                // face1->v[2]->impactVert = true;
            }
            if (impact.self) {
                Face *face2 = get_enclosing_face(sim.cloths[0].mesh, impact_u[1]);
                assert(sim.cloths[0].mesh.faces[face2->index] == face2);
                min = _MAX;
                index = -1;
                for (int v = 0; v < 3; v++) {
                    double d = norm(face2->v[v]->u - impact_u[1]);
                    if (d < min) {
                        min = d;
                        index = v;
                    }
                }
                if (min > arcsim::magic.merge_radius * sqrt(3)) {
                    impactFaces[index] = true;
                    // face2->v[0]->impactVert = true;
                    // face2->v[1]->impactVert = true;
                    // face2->v[2]->impactVert = true;
                }
            }
        }
        bool result = false;
        for (int i = 0; i < sim.cloths[0].mesh.faces.size(); i++) {
            if (impactFaces[i]) {
                Face *face = sim.cloths[0].mesh.faces[i];
                VisualDebugger *vd = VisualDebugger::getInstance();
                vd->addVisualLine3(face->v[0]->node->x, face->v[1]->node->x, Vec3(1, 1, 0), '9');
                vd->addVisualLine3(face->v[1]->node->x, face->v[2]->node->x, Vec3(1, 1, 0), '9');
                vd->addVisualLine3(face->v[2]->node->x, face->v[0]->node->x, Vec3(1, 1, 0), '9');
                if (face->size_min > arcsim::magic.merge_radius * sqrt(3)) {
                    face->size_min = std::max(face->size_min / 2, arcsim::magic.merge_radius * sqrt(3));
                    // face->v[0]->size_min = ::magic.merge_radius*sqrt(3);
                    // face->v[1]->size_min = ::magic.merge_radius*sqrt(3);
                    // face->v[2]->size_min = ::magic.merge_radius*sqrt(3);
                    result = true;
                }
            }
        }
        // for (int i = 0; i < sim.cloths[0].mesh.verts.size(); i++) {
        // 	Vert *vert = sim.cloths[0].mesh.verts[i];
        // 	if (vert->impactVert) {
        // 		cout << "here" << endl;
        // 		if (vert->size_min > ::magic.merge_radius*sqrt(3)) {
        // 			vert->size_min = max(vert->size_min/2, ::magic.merge_radius*sqrt(3));
        // 		}
        // 	}
        // }
        return result;
    }

    void collapse_edges(Simulation &sim) {
        vector<Mesh> old_meshes(sim.cloths.size());
        vector<Mesh *> old_meshes_p(sim.cloths.size()); // for symmetry in separate()
        for (int c = 0; c < sim.cloths.size(); c++) {
            old_meshes[c] = deep_copy(sim.cloths[c].mesh);
            old_meshes_p[c] = &old_meshes[c];
        }
        init_meshes(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes);
        vector<AccelStruct *> obs_accs = create_accel_structs(sim.obstacle_meshes, false);
        vector<Plane> planes = nearest_obstacle_planes(sim.cloths[0].mesh,
                                                       sim.obstacle_meshes);
        collapse_edges(sim.cloths[0], planes, obs_accs);
        destroy_accel_structs(obs_accs);
        for (int c = 0; c < sim.cloths.size(); c++)
            delete_mesh(old_meshes[c]);
    }

    void step_mesh(Mesh &mesh, double dt) {
        for (int n = 0; n < mesh.nodes.size(); n++)
            mesh.nodes[n]->x += mesh.nodes[n]->v * dt;
    }


    void plasticity_step(Simulation &sim) {
        if (!sim.enabled[plasticity])
            return;
        sim.timers[plasticity].tick();
        for (int c = 0; c < sim.cloths.size(); c++) {
            plastic_update(sim.cloths[c]);
            optimize_plastic_embedding(sim.cloths[c]);
        }
        sim.timers[plasticity].tock();
    }

    void strainlimiting_step(Simulation &sim, const vector<Constraint *> &cons) {
        if (!sim.enabled[strainlimiting])
            return;
        sim.timers[strainlimiting].tick();
        vector<Vec3> xold = node_positions(sim.cloth_meshes);
        strain_limiting(sim.cloth_meshes, get_strain_limits(sim.cloths), cons);
        update_velocities(sim.cloth_meshes, xold, sim.step_time);
        sim.timers[strainlimiting].tock();
    }


    void equilibration_step(Simulation &sim) {
        sim.timers[remeshing].tick();
        vector<Constraint *> cons;// = get_constraints(sim, true);
        // double stiff = 1;
        // swap(stiff, ::magic.handle_stiffness);
        for (int c = 0; c < sim.cloths.size(); c++) {
            Mesh &mesh = sim.cloths[c].mesh;
            for (int n = 0; n < mesh.nodes.size(); n++)
                mesh.nodes[n]->acceleration = Vec3(0);
            apply_pop_filter(sim.cloths[c], cons, 1);
        }
        // swap(stiff, ::magic.handle_stiffness);
        sim.timers[remeshing].tock();
        delete_constraints(cons);
        cons = get_constraints(sim, false);

        if (sim.enabled[collision]) {
            sim.timers[collision].tick();
            collision_response_arcsim(sim.cloth_meshes, cons, sim.obstacle_meshes);
            sim.timers[collision].tock();
        }


        delete_constraints(cons);
    }

    void strainzeroing_step(Simulation &sim) {
        sim.timers[strainlimiting].tick();
        vector<Vec2> strain_limits(size<Face>(sim.cloth_meshes), Vec2(1, 1));
        vector<Constraint *> cons =
                proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes,
                                      sim.friction, sim.obs_friction);
        strain_limiting(sim.cloth_meshes, strain_limits, cons);
        delete_constraints(cons);
        sim.timers[strainlimiting].tock();

        // commented this out because collision_response isn't used anymore, we have bh_collision_response, and argus style now ~ george
        if (sim.enabled[collision]) {
            sim.timers[collision].tock();
            collision_response_arcsim(sim.cloth_meshes, vector<Constraint *>(),
                                      sim.obstacle_meshes);
            sim.timers[collision].tock();
        }
    }


// Bridson-Harmon collision step
    void collision_step_bh(Simulation &sim) {

        if (!sim.enabled[collision] || magic.sim_type != "bridsonharmon")
            return;
        sim.timers[collision].tick();


        // Store x and v as they were computed from internal dynamics before collision resolution
        // Also compute the average velocity for the timestep and store it in vBarMid
        for (int m = 0; m < sim.cloth_meshes.size(); m++) {
            Mesh *mesh = sim.cloth_meshes[m];
            for (int n = 0; n < mesh->nodes.size(); n++) {
                Node *node = mesh->nodes[n];
                node->xBar = node->x;
                node->vBar = node->v;
                node->vBarMid = (node->xBar - node->x0) / sim.step_time;
                node->vTildeMid = node->vBarMid;
                node->vPrimeMid = node->vTildeMid;
                node->vFinalMid = node->vPrimeMid;
            }
        }
        for (int m = 0; m < sim.obstacle_meshes.size(); m++) {
            Mesh *mesh = sim.obstacle_meshes[m];
            for (int n = 0; n < mesh->nodes.size(); n++) {
                Node *node = mesh->nodes[n];
                node->xBar = node->x;
                node->vBar = node->v;
                node->vBarMid = (node->xBar - node->x0) / sim.step_time;
                node->vTildeMid = node->vBarMid;
                node->vPrimeMid = node->vTildeMid;
                node->vFinalMid = node->vPrimeMid;
            }
        }

        // Getting constraints
        vector<Constraint *> cons = get_constraints(sim, false);

        // Apply impulses to get vFinalMid for each node
        bool resolvedAnything = collision_response_bh(sim.cloth_meshes, cons, sim.obstacle_meshes, sim.step_time);
        //bool resolvedAnything = false;

        // If the collision response changed anything, then we need to update the final node positions and velocities
        // Otherwise, positions and velocities computed in the physics step are the official end of time step positions and velocities
        if (resolvedAnything) {

            // Compute final positions and velocities
            for (int m = 0; m < sim.cloth_meshes.size(); m++) {
                Mesh *mesh = sim.cloth_meshes[m];
                for (int n = 0; n < mesh->nodes.size(); n++) {
                    Node *node = mesh->nodes[n];
                    node->x = node->x0 + sim.step_time * node->vFinalMid;
                    node->v = node->vBar + (node->x - node->xBar) / sim.step_time;
                    if (n == 0) {
                        std::cout << "node->x: " << node->x << std::endl;
                        std::cout << "node->v: " << node->v << std::endl;
                    }
                }
            }
        }

        // Deleting constraints
        delete_constraints(cons);

        sim.timers[collision].tock();


    }

    void compute_internal_force_density(Cloth &cloth) {
        vector<Vec3> forceVec(cloth.mesh.nodes.size(), Vec3(0, 0, 0));
        add_internal_forces<WS>(cloth, forceVec);
        for (int n = 0; n < cloth.mesh.nodes.size(); n++) {
            Node *node = cloth.mesh.nodes[n];
            node->internal_force_density = forceVec[n] / node->a;
        }
    }

    void remeshing_step(Simulation &sim, bool initializing) {
        if (!sim.enabled[remeshing])
            return;
        // copy old meshes
        vector<Mesh> old_meshes(sim.cloths.size());
        vector<Mesh *> old_meshes_p(sim.cloths.size()); // for symmetry in separate()
        for (int c = 0; c < sim.cloths.size(); c++) {
            old_meshes[c] = deep_copy(sim.cloths[c].mesh);
            old_meshes_p[c] = &old_meshes[c];
        }
        // back up residuals
        // typedef vector<Residual> MeshResidual;
        // vector<MeshResidual> res;
        // if (sim.enabled[plasticity] && !initializing) {
        //     sim.timers[plasticity].tick();
        //     res.resize(sim.cloths.size());
        //     for (int c = 0; c < sim.cloths.size(); c++)
        //         res[c] = back_up_residuals(sim.cloths[c].mesh);
        //     sim.timers[plasticity].tock();
        // }
        // remesh
        sim.timers[remeshing].tick();
        // init_meshes(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes);
        vector<AccelStruct *> obs_accs = create_accel_structs(sim.obstacle_meshes, false);
        for (int c = 0; c < sim.cloths.size(); c++) {
            compute_internal_force_density(sim.cloths[c]);
            deactivate_nodes(sim.cloths[c].mesh.nodes);
            if (magic.fixed_high_res_mesh) {
                static_remesh(sim.cloths[c], obs_accs, true);
            } else {

                // vector<Plane> obs_planes = nearest_obstacle_planes(sim.cloths[c].mesh,
                //                                                sim.obstacle_meshes);
                NearestSearch ns;
                // vector<Plane> new_planes = ns.getNearestPlanes(sim.cloth_meshes, sim.obstacle_meshes);

                if (magic.use_proximity_metric) {
                    CuttingPlaneSet cutting_planes = ns.getCuttingPlanes(sim.cloth_meshes, sim.obstacle_meshes);
                    setPlaneSet(cutting_planes);
                } else {
                    vector<Plane> obs_planes = nearest_obstacle_planes(sim.cloths[c].mesh, sim.obstacle_meshes);
                    setPlanes(obs_planes);
                }
                dynamic_remesh(sim.cloths[c], /*cutting_planes, */sim.enabled[plasticity], obs_accs, true);
            }
            activate_nodes(sim.cloths[c].mesh.nodes);
        }
        ICM_separate(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes, false);
        destroy_accel_structs(obs_accs);
        sim.timers[remeshing].tock();
        // restore residuals
        // if (sim.enabled[plasticity] && !initializing) {
        //     sim.timers[plasticity].tick();
        //     for (int c = 0; c < sim.cloths.size(); c++)
        //         restore_residuals(sim.cloths[c].mesh, old_meshes[c], res[c]);
        //     sim.timers[plasticity].tock();
        // }
        // // separate
        // if (sim.enabled[separation]) {
        //     sim.timers[separation].tick();
        //     separate(sim.cloth_meshes, old_meshes_p, sim.obstacle_meshes);
        //     sim.timers[separation].tock();
        // }
        // // apply pop filter
        // if (sim.enabled[popfilter] && !initializing) {
        //     sim.timers[popfilter].tick();
        //     vector<Constraint*> cons = get_constraints(sim, true);
        //     for (int c = 0; c < sim.cloths.size(); c++)
        //         apply_pop_filter(sim.cloths[c], cons);
        //     delete_constraints(cons);
        //     sim.timers[popfilter].tock();
        // }
        // delete old meshes
        for (int c = 0; c < sim.cloths.size(); c++)
            delete_mesh(old_meshes[c]);
    }

    void update_velocities(vector<Mesh *> &meshes, vector<Vec3> &xold, double dt) {
        double inv_dt = 1 / dt;
#pragma omp parallel for
        for (int n = 0; n < xold.size(); n++) {
            Node *node = get<Node>(n, meshes);
            node->v += (node->x - xold[n]) * inv_dt;
        }
    }

    void compute_avg_velocities(std::vector<Mesh *> &meshes, double dt) {

/*	double inv_dt = 1.0/dt;
#pragma omp parallel for
	for(int m = 0; m < meshes.size(); m++){
		for(int n = 0; n < meshes[m]->nodes.size(); n++){
			Node* node = meshes[m]->nodes[n];
			node->vMid = inv_dt * ( node->x - node->x0 );
		}
	}
	*/

    }

    void update_obstacles(Simulation &sim, bool update_positions) {
        double decay_time = 0.1,
                blend = sim.step_time / decay_time;
        blend = blend / (1 + blend);
        for (int o = 0; o < sim.obstacles.size(); o++) {
            sim.obstacles[o].get_mesh(sim.time);
            // sim.obstacles[o].blend_with_previous(sim.time, sim.step_time, blend);
            if (!update_positions) {
                // put positions back where they were
                Mesh &mesh = sim.obstacles[o].get_mesh();
                for (int n = 0; n < mesh.nodes.size(); n++) {
                    Node *node = mesh.nodes[n];
                    node->v = (node->x - node->x0) / sim.step_time;
                    node->x = node->x0;
                }
            }
        }
    }

// Helper functions

    template<typename Prim>
    int size(const vector<Mesh *> &meshes) {
        int np = 0;
        for (int m = 0; m < meshes.size(); m++) np += get<Prim>(*meshes[m]).size();
        return np;
    }

    template int size<Vert>(const vector<Mesh *> &);

    template int size<Node>(const vector<Mesh *> &);

    template int size<Edge>(const vector<Mesh *> &);

    template int size<Face>(const vector<Mesh *> &);

    template<typename Prim>
    int get_index(const Prim *p,
                  const vector<Mesh *> &meshes) {
        int i = 0;
        for (int m = 0; m < meshes.size(); m++) {
            const vector<Prim *> &ps = get<Prim>(*meshes[m]);
            if (p->index < ps.size() && p == ps[p->index])
                return i + p->index;
            else i += ps.size();
        }
        return -1;
    }

    template int get_index(const Vert *, const vector<Mesh *> &);

    template int get_index(const Node *, const vector<Mesh *> &);

    template int get_index(const Edge *, const vector<Mesh *> &);

    template int get_index(const Face *, const vector<Mesh *> &);

    template<typename Prim>
    Prim *get(int i, const vector<Mesh *> &meshes) {
        for (int m = 0; m < meshes.size(); m++) {
            const vector<Prim *> &ps = get<Prim>(*meshes[m]);
            if (i < ps.size())
                return ps[i];
            else
                i -= ps.size();
        }
        return NULL;
    }

    template Vert *get(int, const vector<Mesh *> &);

    template Node *get(int, const vector<Mesh *> &);

    template Edge *get(int, const vector<Mesh *> &);

    template Face *get(int, const vector<Mesh *> &);

    vector<Vec3> node_positions(const vector<Mesh *> &meshes) {
        vector<Vec3> xs(size<Node>(meshes));
        for (int n = 0; n < xs.size(); n++)
            xs[n] = get<Node>(n, meshes)->x;
        return xs;
    }

    vector<Vec3> node_avg_velocities(const vector<Mesh *> &meshes, double dt) {
        vector<Vec3> vs(size<Node>(meshes));
        for (int n = 0; n < vs.size(); n++) {
            vs[n] = (get<Node>(n, meshes)->x - get<Node>(n, meshes)->x0) / dt;
        }
        return vs;
    }

}
