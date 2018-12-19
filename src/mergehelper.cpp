#include "mergehelper.hpp"
#include <assert.h>
#include <algorithm>
#include "util.hpp"
#include "magic.hpp"
#include "geometry.hpp"
#include "dynamicremesh.hpp"

namespace arcsim {
	void get_VF_impactPos(const Impact &impact, NodalImpact *nodalImpact) {
		assert(impact.type == Impact::VF);
		Vec3 posV = impact.nodes[0]->x;
		Vec3 posF = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2])
					+ impact.nodes[3]->x * abs(impact.w[3]);
		double d[3];
		ImpactPoint *p1 = new ImpactPoint;
		p1->impact = nodalImpact;
		p1->v = Vec3(0);    // p1 is the cloth node, no need to store its velocity.
		ImpactPoint *p2 = new ImpactPoint;
		p2->impact = nodalImpact;
		if (impact.inverted) {
			p1->x = posF;
			p1->u = impact.matPosB;
			for (int v = 1; v <= 3; v++) {
				p1->verts.push_back(impact.verts[v]);
			}
			p2->x = posV;
			p2->v = impact.nodes[0]->v;
			p2->u = impact.matPosA;
			p2->verts.push_back(impact.verts[0]);
		} else {
			p1->x = posV;
			p1->u = impact.matPosA;
			p1->verts.push_back(impact.verts[0]);
			p2->x = posF;
			p2->v = impact.self ? Vec3(0) : impact.nodes[1]->v * abs(impact.w[1]) +
											impact.nodes[2]->v * abs(impact.w[2]) + impact.nodes[3]->v *
																					abs(impact.w[3]);  // If self-impact, p2 is also cloth node, no need to store velocity
			p2->u = impact.matPosB;
			//if (impact.self) {
			for (int v = 1; v <= 3; v++) {
				p2->verts.push_back(impact.verts[v]);
			}
			//}
		}
		p1->type = ImpactPoint::CLOTH;
		p2->type = impact.self ? ImpactPoint::CLOTH : ImpactPoint::OBSTACLE;
		nodalImpact->p1 = p1;
		nodalImpact->p2 = p2;

	}

	void get_EE_impactPos(const Impact &impact, NodalImpact *nodalImpact) {
		assert(impact.type == Impact::EE);
		ImpactPoint *p1 = new ImpactPoint;
		p1->impact = nodalImpact;
		p1->v = Vec3(0);
		p1->x = impact.nodes[0]->x * abs(impact.w[0]) + impact.nodes[1]->x * abs(impact.w[1]);
		p1->u = impact.matPosA;
		p1->type = ImpactPoint::CLOTH;
		p1->verts.push_back(impact.verts[0]);
		p1->verts.push_back(impact.verts[1]);
		ImpactPoint *p2 = new ImpactPoint;
		p2->impact = nodalImpact;
		p2->v = impact.self ? Vec3(0) : impact.nodes[2]->v * abs(impact.w[2]) + impact.nodes[3]->v * abs(impact.w[3]);
		p2->x = impact.nodes[2]->x * abs(impact.w[2]) + impact.nodes[3]->x * abs(impact.w[3]);
		p2->u = impact.matPosB;
		p2->type = impact.self ? ImpactPoint::CLOTH : ImpactPoint::OBSTACLE;
		//if (impact.self) {
		p2->verts.push_back(impact.verts[2]);
		p2->verts.push_back(impact.verts[3]);
		//}
		nodalImpact->p1 = p1;
		nodalImpact->p2 = p2;
	}

	void get_VE_impactPos(const Impact &impact, NodalImpact *nodalImpact) {
		assert(impact.type == Impact::VE);
		Vec3 posV = impact.nodes[0]->x;
		Vec3 posE = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2]);
		ImpactPoint *p1 = new ImpactPoint;
		p1->impact = nodalImpact;
		p1->v = Vec3(0);
		ImpactPoint *p2 = new ImpactPoint;
		p2->impact = nodalImpact;
		if (impact.inverted) {
			p1->x = posE;
			p1->u = impact.matPosB;
			p1->verts.push_back(impact.verts[1]);
			p1->verts.push_back(impact.verts[2]);
			p2->x = posV;
			p2->v = impact.nodes[0]->v;
			p2->u = impact.matPosA;
			p2->verts.push_back(impact.verts[0]);
		} else {
			p1->x = posV;
			p1->u = impact.matPosA;
			p1->verts.push_back(impact.verts[0]);
			p2->x = posE;
			p2->v = impact.self ? Vec3(0) : impact.nodes[1]->v * abs(impact.w[1]) +
											impact.nodes[2]->v * abs(impact.w[2]);
			p2->u = impact.matPosB;
			//if (impact.self) {
			p2->verts.push_back(impact.verts[1]);
			p2->verts.push_back(impact.verts[2]);
			//}
		}
		p1->type = ImpactPoint::CLOTH;
		p2->type = impact.self ? ImpactPoint::CLOTH : ImpactPoint::OBSTACLE;
		nodalImpact->p1 = p1;
		nodalImpact->p2 = p2;
	}

	void get_VV_impactPos(const Impact &impact, NodalImpact *nodalImpact) {
		assert(impact.type == Impact::VV);
		ImpactPoint *p1 = new ImpactPoint;
		p1->impact = nodalImpact;
		p1->x = impact.nodes[0]->x;
		p1->v = Vec3(0);
		p1->u = impact.matPosA;
		p1->type = ImpactPoint::CLOTH;
		p1->verts.push_back(impact.verts[0]);
		ImpactPoint *p2 = new ImpactPoint;
		p2->impact = nodalImpact;
		p2->x = impact.nodes[1]->x;
		p2->v = impact.self ? Vec3(0) : impact.nodes[1]->v;
		p2->u = impact.matPosB;
		p2->type = impact.self ? ImpactPoint::CLOTH : ImpactPoint::OBSTACLE;
		//if (impact.self) {
		p2->verts.push_back(impact.verts[1]);
		//}
		nodalImpact->p1 = p1;
		nodalImpact->p2 = p2;
	}

	void getImpactPos(const Impact &impact, NodalImpact *nodalImpact) {
		switch (impact.type) {
			case Impact::VF : {
				get_VF_impactPos(impact, nodalImpact);
				break;
			}
			case Impact::EE : {
				get_EE_impactPos(impact, nodalImpact);
				break;
			}
			case Impact::VE : {
				get_VE_impactPos(impact, nodalImpact);
				break;
			}
			case Impact::VV : {
				get_VV_impactPos(impact, nodalImpact);
				break;
			}
		}
		nodalImpact->distance = norm(nodalImpact->p1->x - nodalImpact->p2->x);
		nodalImpact->normal = impact.inverted ? -impact.n : impact.n;
	}

	void MergeHelper::mergeImpacts(Cloth &cloth) {
		initNodalImpacts();
		for (int p = 0; p < impactPoints.size(); p++) {
			ImpactPoint *point = impactPoints[p];
			if (moveToCloseNode(point)) {    // impact point moved to a close node
				continue;
			}
			if (moveToCloseEdge(point, cloth)) {    // impact point moved to a close edge and then moved to a close node
				continue;
			}
			insert_nodes(point, cloth);
		}
		// createClusters();
		pruneImpacts();
	}

	void MergeHelper::insertNodes(Cloth &cloth) {
		initNodalImpacts();
		VisualDebugger *vd = VisualDebugger::getInstance();
		for (int p = 0; p < impactPoints.size(); p++) {
			ImpactPoint *point = impactPoints[p];
			if (point->verts.size() > 1) {
				insert_nodes(point, cloth);
				vd->addVisualPoint3(point->x, Vec3(1, 0, 0), 'i');
			}

		}
	}

	void MergeHelper::recycleMemory() {
		// for (int p = 0; p < impactPoints.size(); p++) {
		// 	delete impactPoints[p];
		// }
		for (int i = 0; i < nodalImpacts.size(); i++) {
			delete nodalImpacts[i]->p1;
			delete nodalImpacts[i]->p2;
			delete nodalImpacts[i];
		}
		for (int c = 0; c < clusters.size(); c++) {
			delete clusters[c];
		}
	}

	bool compareDistance(ImpactPoint *p1, ImpactPoint *p2) {
		return (p1->impact->distance) < (p2->impact->distance);
	}

	void MergeHelper::determineInsertOrNot() {
		for (int i = 0; i < nodalImpacts.size(); i++) {
			NodalImpact *nodalImpact = nodalImpacts[i];
			// First impact point
			ImpactPoint *p1 = nodalImpact->p1;
			ImpactPoint *p2 = nodalImpact->p2;
			if (p1->verts.size() > 1) {
				double angle = 0;
				if (p2->verts.size() == 1) {
					angle = getSharpAngle(p2->verts[0]->node);
				} else if (p2->verts.size() == 2) {
					Node *n1 = p2->verts[0]->node;
					Node *n2 = p2->verts[1]->node;
					Edge *e = getCommonEdge(n1, n2);
					angle = getSharpAngle(e);
				} else {
					assert(false);
				}
				if (angle >= arcsim::magic.sharp_point_angle) {
					p1->toInsert = true;
				} else {
					p1->toInsert = false;
				}
			}

			// Second impact point

		}
	}

	void MergeHelper::initNodalImpacts() {
		nodalImpacts.clear();
		impactPoints.clear();
		for (int i = 0; i < impacts.size(); i++) {
			Impact &impact = impacts[i];
			NodalImpact *nodalImpact = new NodalImpact;
			getImpactPos(impact, nodalImpact);
			if (nodalImpact->p1->type == ImpactPoint::CLOTH) {
				impactPoints.push_back(nodalImpact->p1);
				nodalImpact->p1->visited = false;
			}
			if (nodalImpact->p2->type == ImpactPoint::CLOTH) {
				impactPoints.push_back(nodalImpact->p2);
				nodalImpact->p2->visited = false;
			}
			nodalImpacts.push_back(nodalImpact);
		}
		determineInsertOrNot();
		// Sort in distance ascending order.
		sort(impactPoints.begin(), impactPoints.end(), compareDistance);
	}

	Node **get_VF_CloseNodes(Impact &impact) {
		assert(impact.type == Impact::VF);
		Node **nodes = new Node *[2];
		nodes[0] = NULL;
		nodes[1] = NULL;
		if (!impact.inverted && !impact.self) {
			nodes[0] = impact.nodes[0];
			return nodes;
		}
		int index = abs(impact.w[1]) < abs(impact.w[2]) ?
					(abs(impact.w[2]) < abs(impact.w[3]) ? 3 : 2) :
					(abs(impact.w[1]) < abs(impact.w[3]) ? 3 : 1);
		if (impact.inverted) {
			if (norm(impact.matPosB - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
				nodes[0] = impact.nodes[index];
			}
		} else {    // self-contact
			nodes[0] = impact.nodes[0];
			if (norm(impact.matPosB - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
				nodes[1] = impact.nodes[index];
			}
		}
		return nodes;
	}

	Node **get_EE_CloseNodes(Impact &impact) {
		assert(impact.type == Impact::EE);
		Node **nodes = new Node *[2];
		nodes[0] = NULL;
		nodes[1] = NULL;
		int index = abs(impact.w[0]) < abs(impact.w[1]) ? 1 : 0;
		assert(impact.verts[index]->u == impact.nodes[index]->verts[0]->u);
		if (norm(impact.matPosA - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
			nodes[0] = impact.nodes[index];
		}
		if (impact.self) {
			index = abs(impact.w[2]) < abs(impact.w[3]) ? 3 : 2;
			if (norm(impact.matPosB - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
				nodes[1] = impact.nodes[index];
			}
		}
		return nodes;
	}

	Node **get_VE_CloseNodes(Impact &impact) {
		assert(impact.type == Impact::VE);
		Node **nodes = new Node *[2];
		nodes[0] = NULL;
		nodes[1] = NULL;
		if (!impact.inverted && !!impact.self) {
			nodes[0] = impact.nodes[0];
			return nodes;
		}
		int index = abs(impact.w[1]) < abs(impact.w[2]) ? 2 : 1;
		if (impact.inverted) {
			if (norm(impact.matPosB - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
				nodes[0] = impact.nodes[index];
			}
		} else {    // self-contact
			nodes[0] = impact.nodes[0];
			if (norm(impact.matPosB - impact.verts[index]->u) <= arcsim::magic.merge_radius) {
				nodes[1] = impact.nodes[index];
			}
		}
		return nodes;
	}

	Node **get_VV_CloseNodes(Impact &impact) {
		assert(impact.type == Impact::VV);
		Node **nodes = new Node *[2];
		nodes[0] = impact.nodes[0];
		if (impact.self) {
			nodes[1] = impact.nodes[1];
		} else {
			nodes[1] = NULL;
		}
		return nodes;
	}

	Node **getCloseNodes(Impact &impact) {
		switch (impact.type) {
			case Impact::VF : {
				return get_VF_CloseNodes(impact);
			}
			case Impact::EE : {
				return get_EE_CloseNodes(impact);
			}
			case Impact::VE : {
				return get_VE_CloseNodes(impact);
			}
			case Impact::VV : {
				return get_VV_CloseNodes(impact);
			}
            default:
                assert(false);
                return NULL;
		}
	}

	bool MergeHelper::moveToCloseNode(ImpactPoint *point) {
		assert(point->type == ImpactPoint::CLOTH);
		if (point->verts.size() == 1) {    // node impact point
			point->node = point->verts[0]->node;
			return true;
		}
		if (point->verts.size() == 2 || point->verts.size() == 3) {    // Edge/Face impact point
			double min_d = _MAX;
			int index = -1;
			for (int v = 0; v < point->verts.size(); v++) {
				double d = norm(point->u - point->verts[v]->u);
				if (d < min_d) {
					min_d = d;
					index = v;
				}
			}
			// if (min_d < arcsim::magic.merge_radius) {
			// if (!point->toInsert) {
			if (true) {
				point->node = point->verts[index]->node;
				return true;
			}
			Vert *vert = kdTree->getVertWithinRadius(point->u, arcsim::magic.merge_radius, point->verts[index]->kdNode, true);
			if (vert) {
				point->node = vert->node;
				return true;
			}
			return false;
		}
		assert(false);    // Impossible to get here
	}

	bool MergeHelper::moveToCloseEdge(ImpactPoint *point, Cloth &cloth) {
		// if (point->verts.size() == 2) {
		// }
		// if (point->verts.size() == 3) {
		// }
		Face *face = get_enclosing_face(cloth.mesh,
										point->u);    // Temporarily use get_enclosing_face. This is too costly as it loops through all verts.
		assert(face);
		double min_d = _MAX;
		Vec2 min_proj;
		int index = -1;
		for (int v = 0; v < 3; v++) {
			Vec2 proj;
			double d = material_ve_projection(point->u, face->v[v]->u, face->v[(v + 1) % 3]->u, proj);
			if (d < min_d && d >= 0) {
				min_d = d;
				index = v;
				min_proj = proj;
			}
		}
		if (min_d > arcsim::magic.merge_radius) {    // If nearest edge is within radius
			return false;
		}
		point->u = min_proj;
		point->verts.clear();
		point->verts.push_back(face->v[index]);
		point->verts.push_back(face->v[(index + 1) % 3]);
		if (moveToCloseNode(point)) {
			return true;
		}
		return false;
	}

// void MergeHelper::insert_nodes(ImpactPoint* point, vector<Mesh*> &meshes, const vector<Mesh*> &obs_meshes) {

// }

	void MergeHelper::gatherToCloseNodes() {

		vector<Node *> nodes;
		clusters.clear();
		for (int i = 0; i < impacts.size(); i++) {
			Impact &impact = impacts[i];
			Node **cNodes = getCloseNodes(impact);
			if (cNodes[0] != NULL) {
				assert(cNodes[0] != cNodes[1]);
			}
			for (int v = 0; v < 2; v++) {
				ImpactPoint *p = v == 0 ? nodalImpacts[i]->p1 : nodalImpacts[i]->p2;
				if (p->visited) {
					continue;
				}
				Node *n = cNodes[v];
				if (n == NULL) {
					continue;
				}
				int index = find(n, nodes);
				if (index < 0) {
					nodes.push_back(n);
					Cluster *cluster = new Cluster;
					clusters.push_back(cluster);
					cluster->x = n->x;
					cluster->u = n->verts[0]->u;    // TODO: This should be changed when there are seams.
					index = nodes.size() - 1;
					cluster->node = n;
				}
				assert(clusters.size() == nodes.size());
				p->visited = true;
				clusters[index]->points.push_back(p);
			}
		}
	}

// void MergeHelper::createClusters() {

// 	gatherToCloseNodes();
// 	for (int i = 0; i < impactPoints.size(); i++) {
// 		ImpactPoint* pI = impactPoints[i];
// 		if (pI->visited || pI->type != ImpactPoint::CLOTH) {
// 			continue;
// 		}
// 		pI->visited = true;
// 		Cluster *cluster = new Cluster;
// 		clusters.push_back(cluster);
// 		cluster->x = pI->x;
// 		cluster->u = pI->u;
// 		cluster->points.push_back(pI);
// 		for (int j = i + 1; j < impactPoints.size(); j++) {
// 			ImpactPoint* pJ = impactPoints[j];
// 			if (pJ->visited || pJ->type != ImpactPoint::CLOTH) {
// 				continue;
// 			}
// 			if (norm(pI->u - pJ->u) <= ::magic.merge_radius) {
// 				pJ->visited = true;
// 				cluster->points.push_back(pJ);
// 			}
// 		}
// 	}
// 	for (int c = 0; c < clusters.size(); c++) {
// 		Cluster* cluster = clusters[c];
// 		for (int p = 0; p < cluster->points.size(); p++) {
// 			ImpactPoint* point = cluster->points[p];
// 			NodalImpact* impact = point->impact;
// 			if (impact->p1 == point) {
// 				impact->c1 = cluster;
// 			} else {
// 				impact->c2 = cluster;
// 			}
// 		}
// 	}
// }

	void MergeHelper::pruneImpacts() {

		bool *toDelete = new bool[nodalImpacts.size()];
		for (int d = 0; d < nodalImpacts.size(); d++) {
			toDelete[d] = false;
		}

		// Prune those self impacts whose tow impacts points are too close (in the same cluster)
		for (int i = 0; i < nodalImpacts.size(); i++) {
			NodalImpact *imp = nodalImpacts[i];
			if (imp->p1->node == imp->p2->node && imp->p2->type == ImpactPoint::CLOTH) {
				toDelete[i] = true;
			}
		}

		for (int i = 0; i < nodalImpacts.size(); i++) {
			if (toDelete[i]) {
				continue;
			}
			NodalImpact *impI = nodalImpacts[i];
			for (int j = i + 1; j < nodalImpacts.size(); j++) {
				if (toDelete[j]) {
					continue;
				}
				NodalImpact *impJ = nodalImpacts[j];
				if (impI->p2->type != impJ->p2->type) {
					continue;
				}
				double threshold = arcsim::magic.duplication_threshold;
				if (impI->p2->type == ImpactPoint::OBSTACLE) {
					if (impI->p1->node == impJ->p1->node && norm(impI->normal - impJ->normal) < threshold) {
						// toDelete[j] = true;
						double disI = dot((impI->p1->x - (impI->p2->v * dt + impI->p2->x)), impI->normal);
						double disJ = dot((impJ->p1->x - (impJ->p2->v * dt + impJ->p2->x)), impJ->normal);
						if (disI < disJ) {
							toDelete[j] = true;
						} else {
							toDelete[i] = true;
							break;
						}
					}
				} else {
					if ((impI->p1->node == impJ->p1->node) && (impI->p2->node == impJ->p2->node)) {
						if (norm(impI->normal - impJ->normal) < threshold) {
							double disI = dot((impI->p1->x - impI->p2->x), impI->normal);
							double disJ = dot((impJ->p1->x - impJ->p2->x), impJ->normal);
							if (disI < disJ) {
								toDelete[j] = true;
							} else {
								toDelete[i] = true;
								break;
							}
						}
					} else if ((impI->p1->node == impJ->p2->node) && (impI->p2->node == impJ->p1->node)) {
						if (norm(impI->normal + impJ->normal) < threshold) {
							double disI = dot((impI->p1->x - impI->p2->x), impI->normal);
							double disJ = dot((impJ->p1->x - impJ->p2->x), impJ->normal);
							if (disI < disJ) {
								toDelete[j] = true;
							} else {
								toDelete[i] = true;
								break;
							}
						}
					}
				}
			}
		}

		for (int i = nodalImpacts.size() - 1; i >= 0; i--) {
			if (toDelete[i]) {
				delete nodalImpacts[i]->p1;
				delete nodalImpacts[i]->p2;
				delete nodalImpacts[i];
				nodalImpacts.erase(nodalImpacts.begin() + i);
			}
		}
		delete[] toDelete;

		// for (int i = 0; i < nodalImpacts.size(); i++) {
		// 	NodalImpact* imp = nodalImpacts[i];
		// 	cout << "-----" << endl;
		// 	cout << "impact:" << endl;
		// 	cout << "normal:" << imp->normal << endl;
		// 	cout << "cluster1:" << imp->c1 << endl;
		// 	cout << "cluster1:" << (imp->c2 ? imp->c2 : 0) << endl;
		// }

		// cout << "cluster size:" << clusters.size() << endl;
		//    double min = 100;
		//    int ci, cj;
		//    for (int i = 0; i < clusters.size(); i++) {
		//        for (int j = i + 1; j < clusters.size(); j++) {
		//            double dis = norm(clusters[i]->u - clusters[j]->u);
		//            if (min > dis) {
		//                min = dis;
		//                ci = i;
		//                cj = j;
		//            }
		//        }
		//    }
		// cout << "ci:" << clusters[ci] << "," << clusters[ci]->u << endl;
		// for (int i = 0; i < clusters[ci]->points.size(); i++) {
		// 	cout << clusters[ci]->points[i]->impact->normal << endl;
		// }
		// cout << "cj:" << clusters[cj] << ","  << clusters[cj]->u << endl;
		// for (int j = 0; j < clusters[cj]->points.size(); j++) {
		// 	cout << clusters[cj]->points[j]->impact->normal << endl;
		// }
		// cout << min << endl;
		// assert(min > ::magic.merge_radius);
	}

	bool samePrimitive(Node **n1, Node **n2, int num) {
		for (int i = 0; i < num; i++) {
			bool matched = false;
			for (int j = 0; j < num; j++) {
				if (n1[i] == n2[j]) {
					matched = true;
					break;
				}
			}
			if (!matched) {
				return false;
			}
		}
		return true;
	}

	bool isDuplicatedImpacts(Impact &impI, Impact &impJ) {
		if (impI.type != impJ.type) {
			return false;
		}
		if (impI.self != impJ.self) {
			return false;
		}
		switch (impI.type) {
			case Impact::VF: {
				if (impI.nodes[0] == impJ.nodes[0] && samePrimitive(impI.nodes + 1, impJ.nodes + 1, 3)) {
					return true;
				}
				break;
			}
			case Impact::EE: {
				if (samePrimitive(impI.nodes, impJ.nodes, 2) && samePrimitive(impI.nodes + 2, impJ.nodes + 2, 2)) {
					return true;
				}
				if (impI.self && samePrimitive(impI.nodes, impJ.nodes + 2, 2) &&
					samePrimitive(impI.nodes + 2, impJ.nodes, 2)) {
					return true;
				}
				break;
			}
		}
		return false;
	}

	void pruneImpacts(vector<Impact> &impacts) {
		bool *toDelete = new bool[impacts.size()];
		for (int d = 0; d < impacts.size(); d++) {
			toDelete[d] = false;
		}
		for (int i = 0; i < impacts.size(); i++) {
			if (toDelete[i]) {
				continue;
			}
			Impact &impI = impacts[i];
			for (int j = i + 1; j < impacts.size(); j++) {
				if (toDelete[j]) {
					continue;
				}
				Impact &impJ = impacts[j];
				if (isDuplicatedImpacts(impI, impJ)) {
					cout << "duplicated!" << endl;
					print_impacts(impI);
					cout << "----------" << endl;
					print_impacts(impJ);
					exit(-1);
					toDelete[j] = true;
				}
			}
		}
		for (int i = impacts.size() - 1; i >= 0; i--) {
			if (toDelete[i]) {
				impacts.erase(impacts.begin() + i);
			}
		}
		delete[] toDelete;
	}
}
