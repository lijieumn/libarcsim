#include "mergeimpacts.hpp"
#include <algorithm>
#include "magic.hpp"
#include "geometry.hpp"
#include "dynamicremesh.hpp"
#include "util.hpp"


namespace arcsim {

	bool comparePriority(TopologicalImpact *ti1, TopologicalImpact *ti2) {
		return ti1->priority > ti2->priority;
	}

	bool compareTImpactTime(TopologicalImpact *ti1, TopologicalImpact *ti2) {
		return ti1->time < ti2->time;
	}

	bool compareImpactTime(Impact &impact1, Impact &impact2) {
		return impact1.t < impact2.t;
	}

	vector<TopologicalImpact *> createImpactsTopology(vector<Impact> &impacts) {

		vector<TopologicalImpact *> tImpacts;
		for (int i = 0; i < impacts.size(); i++) {
			Impact &impact = impacts[i];
			TopologicalImpact *tImpact = new TopologicalImpact;
			// tImpact->impact = &impact;
			tImpact->time = impact.t;

			if (impact.type == Impact::VF) {
				if (impact.inverted) {

					// Inverted impact so the matpos is going to be point in face, since its obs vert colliding with cloth face
					tImpact->u = impact.matPosB;
					//tImpact->u = impact.nodes[1]->verts[0]->u*abs(impact.w[1]) + impact.nodes[2]->verts[0]->u*abs(impact.w[2])
					//		+ impact.nodes[3]->verts[0]->u*abs(impact.w[3]);

					tImpact->x = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2])
								 + impact.nodes[3]->x * abs(impact.w[3]);

					tImpact->obsX = impact.nodes[0]->x;
					tImpact->obsV = impact.nodes[0]->v;
					tImpact->distance = norm(tImpact->x - tImpact->obsX);
					tImpact->priority = LOW;

				} else {

					// Standard, so matpos is going to be vert in cloth
					tImpact->u = impact.matPosA;
					//tImpact->u = impact.nodes[0]->verts[0]->u;

					tImpact->x = impact.nodes[0]->x;

					tImpact->obsX = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2])
									+ impact.nodes[3]->x * abs(impact.w[3]);
					tImpact->obsV = impact.nodes[1]->v * abs(impact.w[1]) + impact.nodes[2]->v * abs(impact.w[2])
									+ impact.nodes[3]->v * abs(impact.w[3]);
					tImpact->distance = norm(tImpact->x - tImpact->obsX);
					tImpact->priority = HIGH;

				}
			} else if (impact.type == Impact::EE) {
				// Self intersections are not considered for now.
				tImpact->u = impact.matPosA;
				//tImpact->u = impact.nodes[0]->verts[0]->u*abs(impact.w[0]) + impact.nodes[1]->verts[0]->u*abs(impact.w[1]);
				tImpact->x = impact.nodes[0]->x * abs(impact.w[0]) + impact.nodes[1]->x * abs(impact.w[1]);
				tImpact->obsX = impact.nodes[2]->x * abs(impact.w[2]) + impact.nodes[3]->x * abs(impact.w[3]);
				tImpact->obsV = impact.nodes[2]->v * abs(impact.w[2]) + impact.nodes[3]->v * abs(impact.w[3]);
				tImpact->distance = norm(tImpact->x - tImpact->obsX);
				tImpact->priority = NORMAL;
			} else if (impact.type == Impact::VE) {
				if (impact.inverted) {
					tImpact->u = impact.matPosB;
					//tImpact->u = impact.nodes[1]->verts[0]->u*abs(impact.w[1]) + impact.nodes[2]->verts[0]->u*abs(impact.w[2]);
					tImpact->x = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2]);
					tImpact->obsX = impact.nodes[0]->x;
					tImpact->obsV = impact.nodes[0]->v;
					tImpact->distance = norm(tImpact->x - tImpact->obsX);
					tImpact->priority = NORMAL;
				} else {
					tImpact->u = impact.matPosA;
					//tImpact->u = impact.nodes[0]->verts[0]->u;
					tImpact->x = impact.nodes[0]->x;
					tImpact->obsX = impact.nodes[1]->x;
					tImpact->obsV = impact.nodes[1]->v;
					tImpact->distance = norm(tImpact->x - tImpact->obsX);
					tImpact->priority = HIGH;
				}

			} else if (impact.type == Impact::VV) {
				tImpact->u = impact.matPosA;
				//tImpact->u = impact.nodes[0]->verts[0]->u;
				tImpact->x = impact.nodes[0]->x;
				tImpact->obsX = impact.nodes[1]->x;
				tImpact->obsV = impact.nodes[1]->v;
				tImpact->distance = norm(tImpact->x - tImpact->obsX);
				tImpact->priority = HIGH;
			}
			impact.obsX = tImpact->obsX;
			impact.obsV = tImpact->obsV;
			tImpacts.push_back(tImpact);
		}


		double radius = 1e-2;

		for (int i = 0; i < tImpacts.size(); i++) {
			TopologicalImpact *tImpactA = tImpacts[i];
			for (int j = 0; j < tImpacts.size(); j++) {
				if (i == j) {
					continue;
				}
				TopologicalImpact *tImpactB = tImpacts[j];
				// Use the world space to determine if impacts are close
				if (norm(tImpactA->u - tImpactB->u) <= radius) {
					tImpactA->connectedIndex.push_back(j);
					tImpactB->connectedIndex.push_back(i);
				}
			}
		}
		// Sort according to the priority.
		// sort(tImpacts.begin(), tImpacts.end(), compareTImpactTime);
		// sort(impacts.begin(), impacts.end(), compareImpactTime);

		return tImpacts;
	}

	void getMaxIndependentSet(vector<TopologicalImpact *> &tImpacts, vector<Impact> &impacts) {
		bool *toDelete = new bool[tImpacts.size()];
		bool *visited = new bool[tImpacts.size()];
		for (int i = 0; i < tImpacts.size(); i++) {
			toDelete[i] = false;
			visited[i] = false;
		}
		int visitedNum = 0;
		int removed = 0;

		while (visitedNum + removed <= tImpacts.size()) {
			int minIndex = -1;
			double min = 1e10;
			double max = 0;
			for (int i = 0; i < tImpacts.size(); i++) {
				if (toDelete[i] || visited[i]) {
					continue;
				}
				TopologicalImpact *tImpact = tImpacts[i];
				// Vec3 vertical(0, 0, 1);
				// double vertPrj = dot(impacts[i].n, vertical);
				if (tImpact->distance < min) {
					minIndex = i;
					min = tImpact->distance;
				}
				// if (abs(vertPrj) > max) {
				// 	minIndex = i;
				// 	max = abs(vertPrj);
				// }
			}
			if (minIndex == -1) {
				break;
			}
			visited[minIndex] = true;
			visitedNum++;
			TopologicalImpact *tImpact = tImpacts[minIndex];
			Impact &impact = impacts[minIndex];
			impact.adjacentImpacts.clear();
			for (int c = 0; c < tImpact->connectedIndex.size(); c++) {
				if (!toDelete[tImpact->connectedIndex[c]]) {
					removed++;
					impact.adjacentImpacts.push_back(impacts[tImpact->connectedIndex[c]]);
				}
				toDelete[tImpact->connectedIndex[c]] = true;
			}
		}
		// for (int i = 0; i < tImpacts.size(); i++) {
		// 	if (toDelete[i]) {
		// 		continue;
		// 	}
		// 	TopologicalImpact* tImpact = tImpacts[i];
		// 	for (int c = 0; c < tImpact->connectedIndex.size(); c++) {
		// 		toDelete[tImpact->connectedIndex[c]] = true;
		// 	}
		// }
		for (int d = tImpacts.size() - 1; d >= 0; d--) {
			if (toDelete[d]) {
				impacts.erase(impacts.begin() + d);
				tImpacts.erase(tImpacts.begin() + d);
			}
		}
	}

// Get the cloth node that is proximal to an impact.
	Node *get_proximal_node(const Impact &impact) {
		if (impact.type == Impact::VF) {
			if (impact.inverted) {

				// If two of the barycentric coordinates are less than 0.1,
				// the other node is a proximal node.
				int badIndex = 6;    // 3 indices are 1, 2, 3 whose sum is 6.
				int badCount = 0;

				for (int n = 1; n <= 3; n++) {
					if (abs(impact.w[n]) < 0.1) {
						badCount++;
						badIndex -= n;
					}
				}
				if (badCount == 2) {
					return impact.nodes[badIndex];
				}

				// Find the smallest distance between impact point and a node.
				// The node with less than 5e-3 distance is a proximal node.
				Vec3 pos = impact.nodes[1]->x * abs(impact.w[1]) + impact.nodes[2]->x * abs(impact.w[2])
						   + impact.nodes[3]->x * abs(impact.w[3]);
				double smallest = 1;
				badIndex = 0;
				for (int n = 1; n <= 3; n++) {
					double d = norm(impact.nodes[n]->x - pos);
					if (smallest > d) {
						smallest = d;
						badIndex = n;
					}
				}
				if (smallest < 5e-3) {
					return impact.nodes[badIndex];
				}
			} else {
				// If it's a standard VF impact, the vertex is simply the proximal node.
				return impact.nodes[0];
			}
		} else if (impact.type == Impact::EE) {
			// For EE impacts, the node with a bary larger than 0.9 is the proximal node.
			int badIndex = abs(impact.w[0]) > 0.9 ? 0 : (abs(impact.w[1]) > 0.9 ? 1 : -1);
			if (badIndex != -1) {
				return impact.nodes[badIndex];
			}
			// Besides bary coordinates, absolute distance smaller than 5e-3 is also a criterion for proximal node.
			Vec3 pos = impact.nodes[0]->x * abs(impact.w[0]) + impact.nodes[1]->x * abs(impact.w[1]);
			double d1 = norm(impact.nodes[0]->x - pos);
			double d2 = norm(impact.nodes[1]->x - pos);
			badIndex = d1 < d2 ? (d1 < 5e-3 ? 0 : -1) : (d2 < 5e-3 ? 1 : -1);
			if (badIndex != -1) {
				return impact.nodes[badIndex];
			}
		}
		return NULL;
	}

	int get_proximal_node(const TopologicalImpact *tImpact, const Impact &impact) {
		double radius = 1e-2;
		if (impact.type == Impact::VF) {

			// If obstacle vert colliding with cloth face
			if (impact.inverted) {

				Vec3 pos = tImpact->x;
				Vec2 u = tImpact->u;

				double smallest = 1;
				int badIndex = 0;
				for (int n = 1; n <= 3; n++) {

					double d = norm(impact.verts[n]->u - u);
					//double d = norm(impact.nodes[n]->verts[0]->u - u);

					if (smallest > d) {
						smallest = d;
						badIndex = n;
					}
				}
				if (smallest <= radius) {
					return badIndex;
				}
			} else {
				return -1;
			}
		} else if (impact.type == Impact::EE) {
			Vec3 pos = tImpact->x;
			Vec2 u = tImpact->u;

			double d1 = norm(impact.verts[0]->u - u);
			double d2 = norm(impact.verts[1]->u - u);
			//double d1 = norm(impact.nodes[0]->verts[0]->u - u);
			//double d2 = norm(impact.nodes[1]->verts[0]->u - u);

			int badIndex = d1 < d2 ? (d1 < radius ? 0 : -1) : (d2 < radius ? 1 : -1);
			if (badIndex != -1) {
				return badIndex;
			}
		} else if (impact.type == Impact::VE) {
			if (impact.inverted) {
				Vec2 u = tImpact->u;

				double d1 = norm(impact.verts[1]->u - u);
				double d2 = norm(impact.verts[2]->u - u);

				//double d1 = norm(impact.nodes[1]->verts[0]->u - u);
				//double d2 = norm(impact.nodes[2]->verts[0]->u - u);

				int badIndex = d1 < d2 ? (d1 < radius ? 1 : -1) : (d2 < radius ? 2 : -1);
				if (badIndex != -1) {
					return badIndex;
				}
			} else {
				return -1;
			}
		}
		return -1;
	}

	void move_proximal_nodes(vector<TopologicalImpact *> tImpacts, vector<Impact> &impacts, Mesh *mesh) {

		for (int i = 0; i < tImpacts.size(); i++) {
			TopologicalImpact *tImpact = tImpacts[i];
			Impact &impact = impacts[i];
			int proxNode = get_proximal_node(tImpact, impact);

			if (proxNode == -1) {
				continue;
			}

			Node *node = impact.nodes[proxNode];

			if (node->verts.size() > 1) {
				continue;
			}

			bool boundary = false;
			for (int e = 0; e < node->adje.size(); e++) {
				Edge *edge = node->adje[e];
				if (!edge->adjf[0] || !edge->adjf[1]) {
					boundary = true;
					break;
				}
			}
			boundary = true;
			if (!boundary) {
				// Move the node to the impact position and update node's data
				node->verts[0]->u = tImpact->u;
				node->x = tImpact->x;
			} else {
				tImpact->u = node->verts[0]->u;
				tImpact->x = node->x;
			}
			if (impact.type == Impact::VF) {
				if (impact.inverted) {
					// Update the node's physical info
					double b0 = abs(impact.w[1]);
					double b1 = abs(impact.w[2]);
					double b2 = abs(impact.w[3]);
					Node *n0 = impact.nodes[1];
					Node *n1 = impact.nodes[2];
					Node *n2 = impact.nodes[3];
					node->y = b0 * n0->y + b1 * n1->y + b2 * n2->y;
					node->x0 = b0 * n0->x0 + b1 * n1->x0 + b2 * n2->x0;
					node->v = b0 * n0->v + b1 * n1->v + b2 * n2->v;
					node->acceleration = b0 * n0->acceleration + b1 * n1->acceleration + b2 * n2->acceleration;
					int l0 = n0->label, l1 = n1->label, l2 = n2->label;
					node->label = (l0 == l1 && l0 == l2) ? l0 : 0;

					// Update the impact
					impact.type = Impact::VF;
					impact.inverted = false;
					Node *obsNode = impact.nodes[0];
					Face *obsFace = obsNode->verts[0]->adjf[0];
					impact.nodes[1] = obsFace->v[0]->node;
					impact.nodes[2] = obsFace->v[1]->node;
					impact.nodes[3] = obsFace->v[2]->node;
					impact.nodes[0] = node;
					impact.w[0] = 1;
					for (int j = 1; j <= 3; j++) {
						if (impact.nodes[j] == obsNode) {
							impact.w[j] = -1;
						} else {
							impact.w[j] = 0;
						}
					}

					impact.n = -impact.n;
				} else {
					cerr << "Standard impact(VF uninverted) has an proximal node. This couldn't happen!" << endl;
					exit(-1);
				}
			} else if (impact.type == Impact::EE) {
				// Update the node's physical info
				double b0 = abs(impact.w[0]);
				double b1 = abs(impact.w[1]);
				Node *n0 = impact.nodes[0];
				Node *n1 = impact.nodes[1];
				node->y = b0 * n0->y + b1 * n1->y;
				node->x0 = b0 * n0->x0 + b1 * n1->x0;
				node->v = b0 * n0->v + b1 * n1->v;
				node->acceleration = b0 * n0->acceleration + b1 * n1->acceleration;
				node->label = n0->label == n1->label ? n0->label : 0;

				// Update the impact
				impact.type = Impact::VF;
				impact.inverted = false;
				Node *obsNode1 = impact.nodes[2];
				Node *obsNode2 = impact.nodes[3];
				Face *obsFace = NULL;
				for (int e1 = 0; e1 < obsNode1->adje.size(); e1++) {
					for (int e2 = 0; e2 < obsNode2->adje.size(); e2++) {
						if (obsNode1->adje[e1] == obsNode2->adje[e2]) {
							obsFace = obsNode1->adje[e1]->adjf[0];
							break;
						}
					}
				}
				impact.nodes[1] = obsFace->v[0]->node;
				impact.nodes[2] = obsFace->v[1]->node;
				impact.nodes[3] = obsFace->v[2]->node;
				impact.nodes[0] = node;

				impact.w[0] = 1;
				double w1 = impact.w[2], w2 = impact.w[3];
				for (int j = 1; j <= 3; j++) {
					if (impact.nodes[j] == obsNode1) {
						impact.w[j] = w1;
					} else if (impact.nodes[j] == obsNode2) {
						impact.w[j] = w2;
					} else {
						impact.w[j] = 0;
					}
				}

			} else if (impact.type == Impact::VE) {
				if (impact.inverted) {
					double b0 = abs(impact.w[1]);
					double b1 = abs(impact.w[2]);
					Node *n0 = impact.nodes[1];
					Node *n1 = impact.nodes[2];
					node->y = b0 * n0->y + b1 * n1->y;
					node->x0 = b0 * n0->x0 + b1 * n1->x0;
					node->v = b0 * n0->v + b1 * n1->v;
					node->acceleration = b0 * n0->acceleration + b1 * n1->acceleration;
					node->label = n0->label == n1->label ? n0->label : 0;

					// Update the impact
					impact.type = Impact::VF;
					impact.inverted = false;
					Node *obsNode = impact.nodes[0];
					Face *obsFace = obsNode->verts[0]->adjf[0];
					impact.nodes[1] = obsFace->v[0]->node;
					impact.nodes[2] = obsFace->v[1]->node;
					impact.nodes[3] = obsFace->v[2]->node;
					impact.nodes[0] = node;
					impact.w[0] = 1;
					for (int j = 1; j <= 3; j++) {
						if (impact.nodes[j] == obsNode) {
							impact.w[j] = -1;
						} else {
							impact.w[j] = 0;
						}
					}
					impact.n = -impact.n;
				} else {
					cerr << "Standard impact(VE uninverted) has an proximal node. This couldn't happen!" << endl;
					exit(-1);
				}
			}

			// Update the impacts' barycentric cooridnates in adjacent faces.
			// for (int j = 0; j < tImpacts.size(); j++) {
			// 	for (int f = 0; f < node->verts[0]->adjf.size(); f++) {
			// 		Face* face = node->verts[0]->adjf[f];
			// 		TopologicalImpact* adjTImpact = tImpacts[j];
			// 		if (adjTImpact == tImpact) {
			// 			continue;
			// 		}
			// 		Vec3 bary = get_barycentric_coords(adjTImpact->u, face);
			// 		double eps = 1e-6;
			// 		bool inside = (bary[0] >= -eps) && (bary[1] >= -eps) && (bary[2] >= -eps);
			// 		if (inside) {
			// 			Impact* adjImpact = adjTImpact->impact;
			// 			if (adjImpact->type == Impact::VF) {
			// 				if (adjImpact->inverted) {
			// 					adjImpact->w[1] = bary[0];
			// 					adjImpact->w[2] = bary[1];
			// 					adjImpact->w[3] = bary[2];
			// 				}
			// 			} else {
			// 				adjImpact->w[0] = bary
			// 			}

			// 			break;
			// 		}
			// 	}

			// }

		}
	}

	Edge *get_proximal_edge(const TopologicalImpact *tImpact, const Impact &impact) {
		double radius = 1e-2;
		if (impact.type == Impact::VF) {
			if (impact.inverted) {
				double smallest = 1;
				int badIndex = 0;
				for (int n = 0; n < 3; n++) {
					Vec2 u1 = impact.nodes[n + 1]->verts[0]->u;
					Vec2 u2 = impact.nodes[(n + 1) % 3 + 1]->verts[0]->u;
					Vec2 e = u1 - u2;
					Vec2 u = tImpact->u - u1;
					double proj = dot(e, u) / norm(e);
					double d = sqrt(pow(norm(u), 2) - pow(proj, 2));
					if (smallest > d) {
						smallest = d;
						badIndex = n;
					}
				}
				if (smallest <= radius) {
					Node *n1 = impact.nodes[badIndex + 1];
					Node *n2 = impact.nodes[(badIndex + 1) % 3 + 1];
					for (int e = 0; e < n1->adje.size(); e++) {
						Edge *edge = n1->adje[e];
						if (edge->n[0] == n2 || edge->n[1] == n2) {
							return edge;
						}
					}
				}
			}
		}
		return NULL;
	}

	void move_to_proximal_edge(vector<TopologicalImpact *> tImpacts, vector<Impact> &impacts, Mesh *mesh) {
		for (int i = 0; i < tImpacts.size(); i++) {
			Impact &impact = impacts[i];
			if (impact.type != Impact::VF) {
				continue;
			}
			TopologicalImpact *tImpact = tImpacts[i];
			Edge *edge = get_proximal_edge(tImpact, impact);
			if (edge == NULL) {
				continue;
			}

			if (edge_vert(edge, 0, 0) != edge_vert(edge, 1, 0)) {
				continue;
			}

			if (edge->n[0]->verts.size() > 0 || edge->n[1]->verts.size() > 0) {
				continue;
			}

			// bool boundary = false;
			// if (!edge->adjf[0] || !edge->adjf[1]) {
			// 	boundary = true;
			// }
			if (true) {
				impact.type = Impact::EE;
				impact.inverted = false;
				Vec3 x = tImpact->x;
				Vec3 x1 = edge->n[0]->x;
				Vec3 x2 = edge->n[1]->x;
				double edgeLength = norm(x2 - x1);
				double proj1 = dot(x - x1, x2 - x1) / edgeLength;
				double proj2 = edgeLength - proj1;
				double bary1 = proj1 / edgeLength;
				double bary2 = proj2 / edgeLength;

				Node *obsNode = impact.nodes[0];
				Edge *obsEdge = obsNode->adje[0];
				int index = obsEdge->n[0] == obsNode ? 1 : 0;
				Node *obsNode2 = obsEdge->n[index];
				impact.nodes[0] = edge->n[0];
				impact.w[0] = bary2;
				impact.nodes[1] = edge->n[1];
				impact.w[1] = bary1;
				impact.nodes[2] = obsNode;
				impact.w[2] = -1;
				impact.nodes[3] = obsNode2;
				impact.w[3] = 0;
				impact.n = -impact.n;
				impact.matPosA = tImpact->u;

				tImpact->x = x1 * bary1 + x2 * bary2;
				tImpact->u = edge->n[0]->verts[0]->u * bary1 + edge->n[1]->verts[0]->u * bary2;


			}
		}
	}

	void merge_proximal_impacts(vector<Impact> &impacts, Mesh *mesh) {
		vector<TopologicalImpact *> tImpacts = createImpactsTopology(impacts);
		getMaxIndependentSet(tImpacts, impacts);
		move_proximal_nodes(tImpacts, impacts, mesh);
		move_to_proximal_edge(tImpacts, impacts, mesh);

	}
}