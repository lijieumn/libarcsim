#ifndef MERGEHELPER_HPP
#define MERGEHELPER_HPP

#include "collision.hpp"
#include "kdtree.hpp"

namespace arcsim {
	struct ImpactPoint;
	struct Cluster {
		Vec3 x;
		Vec2 u;
		Node *node = NULL;
		vector<ImpactPoint *> points;    // contact points contained.
	};


	struct NodalImpact {
		double distance;
		Vec3 normal;

		ImpactPoint *p1 = NULL;
		ImpactPoint *p2 = NULL;

		// Cluster* c1 = NULL;
		// Cluster* c2 = NULL;
	};

	struct ImpactPoint {
		enum Type {    // 3 different types of impact point:
			// 1. impact point on cloth; 2. impact point on obstacle. 3. node that close to impact (not a impact point)
					CLOTH, OBSTACLE
		} type;
		vector<const Vert *> verts;    // Enclosing verts. Face impact has 3 verts, edge impact has 2 and node impact has 1.
		Node *node = NULL;    //	impact node.
		Vec3 x;
		Vec2 u;
		Vec3 v;
		NodalImpact *impact = NULL;
		bool visited;
		bool toInsert;
	};

	struct MergeHelper {
		vector<Impact> impacts;
		vector<ImpactPoint *> impactPoints;
		vector<NodalImpact *> nodalImpacts;
		vector<Cluster *> clusters;
		KDTree *kdTree;
		double dt;

		// Below are to be modified ***


		// vector<NodalImpact> nodalImpacts;
		vector<pair<double, int> > distanceAndIndex;

		MergeHelper(vector<Impact> &_impacts, KDTree *_tree, double _dt) : impacts(_impacts), kdTree(_tree), dt(_dt) {}

		void mergeImpacts(Cloth &cloth);

		void initNodalImpacts();

		void determineInsertOrNot();

		void gatherToCloseNodes();

		bool moveToCloseNode(ImpactPoint *point);

		bool moveToCloseEdge(ImpactPoint *point, Cloth &cloth);

		// void createClusters();
		void pruneImpacts();

		void insertNodes(Cloth &cloth);

		void recycleMemory();
		// void insert_nodes(ImpactPoint* point, vector<Mesh*> &meshes, const vector<Mesh*> &obs_meshes);
	};

	void pruneImpacts(vector<Impact> &impacts);
}
#endif