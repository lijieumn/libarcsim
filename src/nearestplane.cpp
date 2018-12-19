#include "nearestplane.hpp"
#include "collisionutil.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "visualdebug.hpp"
#include "timer.hpp"
#include "omp.h"
#include <algorithm>

namespace arcsim {
	static vector<TangentialPlane> planes;
	static vector<double> distances;
	static vector<bool> isStandard;
	static vector<int> pieceIndices;    // Indicate which cloth piece the node is in
	static double dmin;

	static vector<Mesh *> clothMeshes;
	static vector<Mesh *> obsMeshes;

	bool isClothNode(const Node *node) {
		for (int m = 0; m < clothMeshes.size(); m++) {
			Mesh *mesh = clothMeshes[m];
			for (int i = 0; i < mesh->nodes.size(); i++) {
				if (node == mesh->nodes[i]) {
					return true;
				}
			}
		}
		return false;
	}

	bool isClothFace(const Face *face) {
		for (int m = 0; m < clothMeshes.size(); m++) {
			Mesh *mesh = clothMeshes[m];
			for (int i = 0; i < mesh->faces.size(); i++) {
				if (face == mesh->faces[i]) {
					return true;
				}
			}
		}
		return false;
	}

	void updatePieceIndices(const std::vector<Mesh *> &meshes) {
		uint32_t pieceCount = 0;
		for (int m = 0; m < meshes.size(); m++) {
			Mesh *mesh = meshes[m];
			vector<Node *> stack;
			if (mesh->nodes.empty()) {
				continue;
			}
			for (int n = 0; n < mesh->nodes.size(); n++) {
				Node *node = mesh->nodes[n];
				if (pieceIndices[n] >= 0) {
					continue;
				}
				stack.push_back(node);
				while (!stack.empty()) {
					Node *node = stack.back();
					stack.pop_back();
					pieceIndices[node->index] = pieceCount;
					for (int i = 0; i < node->adje.size(); i++) {
						Edge *edge = node->adje[i];
						Node *other = edge->n[0] == node ? edge->n[1] : edge->n[0];
						if (pieceIndices[other->index] == -1) {
							stack.push_back(other);
						}
					}
				}
				pieceCount++;
			}
		}
	}

	double getMaterialDistance(const Node *node, const Face *face) {
		double INF = 0xFFFFFFFF;
		int i0 = node->index;
		int i1 = face->v[0]->node->index;
		if (pieceIndices[i0] != pieceIndices[i1]) {
			return INF;
		}
		double materialDistance = INF;
		for (int i = 0; i < 3; i++) {
			Vert *v0 = face->v[i];
			Vert *v1 = face->v[(i + 1) % 3];
			Vec2 proj;
			for (int v = 0; v < node->verts.size(); v++) {
				double md = material_ve_projection(node->verts[v]->u, v0->u, v1->u, proj);
				if (materialDistance > md) {
					materialDistance = md;
				}
			}
		}
		return materialDistance;
	}

	double getMaterialDistance(const Edge *edge0, const Edge *edge1) {
		double INF = 0xFFFFFFFF;
		int i0 = edge0->n[0]->index;
		int i1 = edge1->n[1]->index;
		if (pieceIndices[i0] != pieceIndices[i1]) {
			return INF;
		}
		double materialDistance = INF;
		Node *n0 = edge0->n[0];
		Node *n1 = edge0->n[1];
		Node *n2 = edge1->n[0];
		Node *n3 = edge1->n[1];
		for (int i = 0; i < n0->verts.size(); i++) {
			for (int j = 0; j < n2->verts.size(); j++) {
				double d = norm(n0->verts[i]->u - n2->verts[j]->u);
				if (materialDistance > d) {
					materialDistance = d;
				}
			}
		}
		for (int i = 0; i < n0->verts.size(); i++) {
			for (int j = 0; j < n3->verts.size(); j++) {
				double d = norm(n0->verts[i]->u - n3->verts[j]->u);
				if (materialDistance > d) {
					materialDistance = d;
				}
			}
		}
		for (int i = 0; i < n1->verts.size(); i++) {
			for (int j = 0; j < n2->verts.size(); j++) {
				double d = norm(n1->verts[i]->u - n2->verts[j]->u);
				if (materialDistance > d) {
					materialDistance = d;
				}
			}
		}
		for (int i = 0; i < n1->verts.size(); i++) {
			for (int j = 0; j < n3->verts.size(); j++) {
				double d = norm(n1->verts[i]->u - n3->verts[j]->u);
				if (materialDistance > d) {
					materialDistance = d;
				}
			}
		}
		// Vert* v0 = edge0->n[0]->verts[0];
		// Vert* v1 = edge0->n[1]->verts[0];
		// Vert* v2 = edge1->n[0]->verts[0];
		// Vert* v3 = edge1->n[1]->verts[0];
		// double d1 = norm(v0->u - v2->u);
		// double d2 = norm(v0->u - v3->u);
		// double d3 = norm(v1->u - v2->u);
		// double d4 = norm(v1->u - v3->u);
		return materialDistance;
	}

	static vector<CuttingPlanes> *nodePlanes;
	static vector<CuttingPlanes> *edgePlanes;
	static vector<CuttingPlanes> *facePlanes;

	void updatePlaneNew(const Node *node, const Face *face, bool first, bool second) {
		int t = omp_get_thread_num();
		Node *tNodes[3];
		for (int i = 0; i < 3; i++) {
			tNodes[i] = face->v[i]->node;
			if (node == tNodes[i]) {
				return;
			}
		}
		if (!overlap(node_box(node, false), face_box(face, false), arcsim::dmin)) {
			return;
		}
		Vec3 n;
		double w[4];
		double d = signed_vf_distance(node->x, tNodes[0]->x, tNodes[1]->x, tNodes[2]->x, &n, w);
		bool inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);


		if (fabs(d) < arcsim::dmin && inside) {
			if (first) {
				int index = node->index;
				CuttingPlane plane;
				plane.pos = tNodes[0]->x * (-w[1]) + tNodes[1]->x * (-w[2]) + tNodes[2]->x * (-w[3]);
				Vec3 dir = normalize(face->n);
				plane.directions.push_back(dir);
				arcsim::nodePlanes[t][index].push_back(plane);
			}
			if (second) {
				int index = face->index;
				CuttingPlane plane;
				plane.pos = node->x;
				for (int v = 0; v < node->verts.size(); v++) {
					Vert *vert = node->verts[v];
					for (int f = 0; f < vert->adjf.size(); f++) {
						Face *af = vert->adjf[f];
						plane.directions.push_back(af->n);
					}
				}
				arcsim::facePlanes[t][index].push_back(plane);
			}
		}
	}

	void updatePlaneNew(const Edge *edge0, const Edge *edge1, bool first, bool second) {
		if (edge0 == edge1) {
			return;
		}
		int t = omp_get_thread_num();
		Node *nodes[4];
		nodes[0] = edge0->n[0];
		nodes[1] = edge0->n[1];
		nodes[2] = edge1->n[0];
		nodes[3] = edge1->n[1];
		if (nodes[0] == nodes[2] || nodes[0] == nodes[3] || nodes[1] == nodes[2] || nodes[1] == nodes[3]) {
			return;
		}
		if (getCommonEdge(nodes[0], nodes[2]) || getCommonEdge(nodes[0], nodes[3])
			|| getCommonEdge(nodes[1], nodes[2]) || getCommonEdge(nodes[1], nodes[3])) {
			return;
		}
		if (!overlap(edge_box(edge0, false), edge_box(edge1, false), arcsim::dmin)) {
			return;
		}
		Vec3 n;
		double w[4];
		double d = signed_ee_distance(nodes[0]->x, nodes[1]->x, nodes[2]->x, nodes[3]->x, &n, w);
		bool inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);

		if (fabs(d) < arcsim::dmin && inside) {
			if (first) {
				int index = edge0->index;
				CuttingPlane plane;
				plane.pos = (-w[2]) * nodes[2]->x + (-w[3]) * nodes[3]->x;
				for (int f = 0; f < 2; f++) {
					Face *af = edge1->adjf[f];
					if (af) {
						plane.directions.push_back(af->n);
					}
				}
				arcsim::edgePlanes[t][index].push_back(plane);
			}
			if (second) {
				int index = edge1->index;
				CuttingPlane plane;
				plane.pos = w[0] * nodes[0]->x + w[1] * nodes[1]->x;
				for (int f = 0; f < 2; f++) {
					Face *af = edge0->adjf[f];
					if (af) {
						plane.directions.push_back(af->n);
					}
				}
				arcsim::edgePlanes[t][index].push_back(plane);
			}
		}
	}

	void NearestSearch::findProximities(const Face *face0, const Face *face1) {
		if (face0 == face1) {
			return;
		}
		if (getCommonEdge(face0, face1)) {
			return;
		}
		// bool first = isClothFace(face0);
		bool first = face0->v[0]->node->inMesh;
		// bool second = isClothFace(face1);
		bool second = face1->v[0]->node->inMesh;
		if (arcsim::magic.use_representative_triangles) {
			for (int n = 0; n < face0->rNode.size(); n++) {
				Node *node = face0->rNode[n];
				updatePlaneNew(node, face1, first, second);
			}
		} else {
			for (int i = 0; i < 3; i++) {    // cloth node, obstacle face
				Node *node = face0->v[i]->node;
				updatePlaneNew(node, face1, first, second);
			}
		}

		if (arcsim::magic.use_representative_triangles) {
			for (int n = 0; n < face1->rNode.size(); n++) {
				Node *node = face1->rNode[n];
				updatePlaneNew(node, face0, second, first);
			}
		} else {
			for (int i = 0; i < 3; i++) {    // cloth face, obstacle node
				Node *node = face1->v[i]->node;
				updatePlaneNew(node, face0, second, first);
			}
		}

		if (arcsim::magic.use_representative_triangles) {
			for (int e0 = 0; e0 < face0->rEdge.size(); e0++) {
				Edge *edge0 = face0->rEdge[e0];
				for (int e1 = 0; e1 < face1->rEdge.size(); e1++) {
					Edge *edge1 = face1->rEdge[e1];
					updatePlaneNew(edge0, edge1, first, second);
				}
			}
		} else {
			for (int e0 = 0; e0 < 3; e0++) {    // Edge edge
				Edge *edge0 = face0->adje[e0];
				for (int e1 = 0; e1 < 3; e1++) {
					Edge *edge1 = face1->adje[e1];
					updatePlaneNew(edge0, edge1, first, second);
				}
			}
		}
	}

	CuttingPlaneSet NearestSearch::getCuttingPlanes(const vector<Mesh *> &meshes, const vector<Mesh *> &obs_meshes) {
		int nthreads = omp_get_max_threads();
		if (!arcsim::nodePlanes) {
			arcsim::nodePlanes = new vector<CuttingPlanes>[nthreads];
		}
		if (!arcsim::edgePlanes) {
			arcsim::edgePlanes = new vector<CuttingPlanes>[nthreads];
		}
		if (!arcsim::facePlanes) {
			arcsim::facePlanes = new vector<CuttingPlanes>[nthreads];
		}
		for (int t = 0; t < nthreads; t++) {
			arcsim::nodePlanes[t].clear();
			arcsim::edgePlanes[t].clear();
			arcsim::facePlanes[t].clear();
			arcsim::nodePlanes[t].resize(meshes[0]->nodes.size());
			arcsim::edgePlanes[t].resize(meshes[0]->edges.size());
			arcsim::facePlanes[t].resize(meshes[0]->faces.size());
		}

		arcsim::dmin = magic.proximity_min;
		vector<AccelStruct *> accs = create_accel_structs(meshes, false);
		vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, false);
		// updatePieceIndices(meshes);
		if (arcsim::magic.use_stack_overloop_face_check) {
			find_overlapping_faces(accs, obs_accs, arcsim::magic.projection_thickness, findProximities);
		} else {
			for_overlapping_faces(accs, obs_accs, arcsim::magic.projection_thickness, findProximities);
		}
		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);

		CuttingPlaneSet planeSet;
		planeSet.nodePlanes.resize(meshes[0]->nodes.size());
		planeSet.edgePlanes.resize(meshes[0]->edges.size());
		planeSet.facePlanes.resize(meshes[0]->faces.size());
		for (int t = 0; t < nthreads; t++) {
			for (int n = 0; n < meshes[0]->nodes.size(); n++) {
				append(planeSet.nodePlanes[n], arcsim::nodePlanes[t][n]);
			}
			for (int e = 0; e < meshes[0]->edges.size(); e++) {
				append(planeSet.edgePlanes[e], arcsim::edgePlanes[t][e]);
			}
			for (int f = 0; f < meshes[0]->faces.size(); f++) {
				append(planeSet.facePlanes[f], arcsim::facePlanes[t][f]);
			}
		}
		if (arcsim::nodePlanes) {
			delete[] arcsim::nodePlanes;
			arcsim::nodePlanes = NULL;
		}
		if (arcsim::edgePlanes) {
			delete[] arcsim::edgePlanes;
			arcsim::edgePlanes = NULL;
		}
		if (arcsim::facePlanes) {
			delete[] arcsim::facePlanes;
			arcsim::facePlanes = NULL;
		}

		return planeSet;
	}

// Currently only considered one mesh
	vector<TangentialPlane>
	NearestSearch::getNearestPlanes(const std::vector<Mesh *> &meshes, const vector<Mesh *> &obs_meshes) {
		clothMeshes = meshes;
		obsMeshes = obs_meshes;
		// ::dmin = 10*::magic.repulsion_thickness;
		arcsim::dmin = magic.proximity_min;
		vector<AccelStruct *> accs = create_accel_structs(meshes, false);
		vector<AccelStruct *> obs_accs = create_accel_structs(obs_meshes, false);
		arcsim::planes.clear();
		arcsim::distances.clear();
		arcsim::isStandard.clear();
		arcsim::pieceIndices.clear();
		arcsim::planes.resize(meshes[0]->nodes.size(), make_pair(Vec3(0), Vec3(0)));
		arcsim::distances.resize(arcsim::planes.size(), arcsim::dmin);
		arcsim::isStandard.resize(arcsim::planes.size(), false);
		arcsim::pieceIndices.resize(arcsim::planes.size(), -1);
		updatePieceIndices(meshes);
		if (arcsim::magic.use_stack_overloop_face_check) {
			find_overlapping_faces(accs, obs_accs, arcsim::magic.projection_thickness, findCloseFaces);
		} else {
			for_overlapping_faces(accs, obs_accs, arcsim::magic.projection_thickness, findCloseFaces);
		}
		destroy_accel_structs(accs);
		destroy_accel_structs(obs_accs);
		return arcsim::planes;
	}

	void updatePlane(const Node *node, const Face *face, bool first, bool second) {
		Node *tNodes[3];
		for (int i = 0; i < 3; i++) {
			tNodes[i] = face->v[i]->node;
			if (node == tNodes[i]) {
				return;
			}
		}
		// for (int i = 0; i < 3; i++) {
		// 	Edge* edge = face->adje[i];
		// 	Vert* v0 = edge_opp_vert(edge, 0);
		// 	Vert* v1 = edge_opp_vert(edge, 1);
		// 	if (v0 && v0->node == node) {
		// 		return;
		// 	}
		// 	if (v1 && v1->node == node) {
		// 		return;
		// 	}
		// }
		Vec3 n;
		double w[4];
		double d = unsigned_vf_distance(node->x, tNodes[0]->x, tNodes[1]->x, tNodes[2]->x, &n, w);
		if ((first && second) &&
			d >= arcsim::magic.false_proximity_threshold * getMaterialDistance(node, face)) {    // Avoid self false proximity
			return;
		}
		if (first) {    // node is in cloth mesh
			uint32_t index = node->index;
			if (d < arcsim::distances[index] || ((!arcsim::isStandard[index]) && d < arcsim::dmin)) {
				arcsim::distances[index] = d;
				Vec3 x = abs(w[1]) * tNodes[0]->x + abs(w[2]) * tNodes[1]->x + abs(w[3]) * tNodes[2]->x;
				arcsim::planes[index].first = x;
				arcsim::planes[index].second = normalize(node->x - x);
				arcsim::isStandard[index] = true;
			}
		}
		if (second) {    // face is in cloth mesh
			for (int i = 0; i < 3; i++) {
				uint32_t index = tNodes[i]->index;
				if (isStandard[index]) {    // captured by standard case
					continue;
				}
				if (d < arcsim::distances[index]) {
					arcsim::distances[index] = d;
					arcsim::planes[index].first = node->x;
					arcsim::planes[index].second = normalize(tNodes[i]->x - node->x);
				}
			}
		}
	}

	void updatePlane(const Edge *edge0, const Edge *edge1, bool first, bool second) {
		if (edge0 == edge1) {
			return;
		}
		Node *nodes[4];
		nodes[0] = edge0->n[0];
		nodes[1] = edge0->n[1];
		nodes[2] = edge1->n[0];
		nodes[3] = edge1->n[1];
		if (nodes[0] == nodes[2] || nodes[0] == nodes[3] || nodes[1] == nodes[2] || nodes[1] == nodes[3]) {
			return;
		}
		if (getCommonEdge(nodes[0], nodes[2]) || getCommonEdge(nodes[0], nodes[3])
			|| getCommonEdge(nodes[1], nodes[2]) || getCommonEdge(nodes[1], nodes[3])) {
			return;
		}
		Vec3 n;
		double w[4];
		double d = unsigned_ee_distance(nodes[0]->x, nodes[1]->x, nodes[2]->x, nodes[3]->x, &n, w);
		if ((first && second) && d >= arcsim::magic.false_proximity_threshold * getMaterialDistance(edge0, edge1)) {
			return;
		}
		Vec3 x1 = w[0] * nodes[0]->x + w[1] * nodes[1]->x;
		Vec3 x2 = abs(w[2]) * nodes[2]->x + abs(w[3]) * nodes[3]->x;
		if (first) {    // The first edge is in cloth
			for (int i = 0; i < 2; i++) {
				uint32_t index = nodes[i]->index;
				if (isStandard[index]) {    // captured by standard case
					continue;
				}
				if (d < arcsim::distances[index]) {
					arcsim::distances[index] = d;
					arcsim::planes[index].first = x2;
					arcsim::planes[index].second = normalize(nodes[i]->x - x2);
				}
			}
		}
		if (second) {    // The second edge is in cloth
			for (int i = 2; i < 4; i++) {
				uint32_t index = nodes[i]->index;
				if (isStandard[index]) {    // captured by standard case
					continue;
				}
				if (d < arcsim::distances[index]) {
					arcsim::distances[index] = d;
					arcsim::planes[index].first = x1;
					arcsim::planes[index].second = normalize(nodes[i]->x - x1);
				}
			}
		}
	}

	void NearestSearch::findCloseFaces(const Face *face0, const Face *face1) {
		if (face0 == face1) {
			return;
		}
		if (getCommonEdge(face0, face1)) {
			return;
		}
		// bool first = isClothFace(face0);
		bool first = face0->v[0]->node->inMesh;
		// bool second = isClothFace(face1);
		bool second = face1->v[0]->node->inMesh;
		if (arcsim::magic.use_representative_triangles) {
			for (int n = 0; n < face0->rNode.size(); n++) {
				Node *node = face0->rNode[n];
				updatePlane(node, face1, first, second);
			}
		} else {
			for (int i = 0; i < 3; i++) {    // cloth node, obstacle face
				Node *node = face0->v[i]->node;
				updatePlane(node, face1, first, second);
			}
		}

		if (arcsim::magic.use_representative_triangles) {
			for (int n = 0; n < face1->rNode.size(); n++) {
				Node *node = face1->rNode[n];
				updatePlane(node, face0, second, first);
			}
		} else {
			for (int i = 0; i < 3; i++) {    // cloth face, obstacle node
				Node *node = face1->v[i]->node;
				updatePlane(node, face0, second, first);
			}
		}

		if (arcsim::magic.use_representative_triangles) {
			for (int e0 = 0; e0 < face0->rEdge.size(); e0++) {
				Edge *edge0 = face0->rEdge[e0];
				for (int e1 = 0; e1 < face1->rEdge.size(); e1++) {
					Edge *edge1 = face1->rEdge[e1];
					updatePlane(edge0, edge1, first, second);
				}
			}
		} else {
			for (int e0 = 0; e0 < 3; e0++) {    // Edge edge
				Edge *edge0 = face0->adje[e0];
				for (int e1 = 0; e1 < 3; e1++) {
					Edge *edge1 = face1->adje[e1];
					updatePlane(edge0, edge1, first, second);
				}
			}
		}
	}
}