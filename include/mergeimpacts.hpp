#ifndef MERGEIMPACTS_HPP
#define MERGEIMPACTS_HPP

#include "collision.hpp"

namespace arcsim {
	enum Priority {
		LOW, NORMAL, HIGH
	};

	struct TopologicalImpact {

		Vec2 u;    // material space
		Vec3 x;    // world space
		Vec3 obsX;
		Vec3 obsV;
		double distance;
		Face *enclosingFace;
		int index;
		vector<int> connectedIndex;
		Priority priority;
		double time;
	};

// vector<TopologicalImpact*> createImpactsTopology(vector<Impact> impacts);

// void getMaxIndependentSet(vector<TopologicalImpact*>& tImpacts);

	void merge_proximal_impacts(vector<Impact> &impacts, Mesh *mesh);
}

#endif