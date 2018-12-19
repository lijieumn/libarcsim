#include "kdtree.hpp"
#include "util.hpp"
#include <stdlib.h>
#include <time.h>
#include <algorithm>

namespace arcsim {

	void KDBox::reset() {
		for (int d = 0; d < dimension; d++) {
			min[d] = _MAX;
			max[d] = _MIN;
		}
	}

	KDBox KDBox::operator+(const KDBox &box) {
		KDBox sum;
		for (int d = 0; d < dimension; d++) {
			sum.min[d] = std::min(min[d], box.min[d]);
			sum.max[d] = std::max(max[d], box.max[d]);
		}
		return sum;
	}

	KDBox &KDBox::operator+=(const KDBox &box) {
		for (int d = 0; d < dimension; d++) {
			min[d] = std::min(min[d], box.min[d]);
			max[d] = std::max(max[d], box.max[d]);
		}
		return *this;
	}

	KDBox &KDBox::operator+=(const Vec2 &u) {
		assert(dimension == 2);
		for (int d = 0; d < 2; d++) {
			min[d] = std::min(min[d], u[d]);
			max[d] = std::max(max[d], u[d]);
		}
		return *this;
	}

	bool KDBox::operator==(const KDBox &box) {
		for (int d = 0; d < dimension; d++) {
			if (max[d] != box.max[d] || min[d] != box.min[d]) {
				return false;
			}
		}
		return true;
	}

	double KDBox::mid(unsigned int axis) {
		return (max[axis] + min[axis]) / 2;
	}

	bool KDBox::inside(const Vec2 &u) {
		assert(dimension == 2);
		double eps = 1e-6;
		if (u[0] - min[0] >= -eps && u[0] - max[0] <= eps &&
			u[1] - min[1] >= -eps && u[1] - max[1] <= eps) {
			return true;
		}
		return false;
	}

	double KDBox::distance(const Vec2 &u) {
		if (u[0] < min[0] && u[1] < min[1]) {
			return norm(u - Vec2(min[0], min[1]));
		}
		if (u[0] < min[0] && u[1] > max[1]) {
			return norm(u - Vec2(min[0], max[1]));
		}
		if (u[0] > max[0] && u[1] > max[1]) {
			return norm(u - Vec2(max[0], max[1]));
		}
		if (u[0] > max[0] && u[1] < min[1]) {
			return norm(u - Vec2(max[0], min[1]));
		}
		if (u[0] < min[0]) {
			return abs(u[0] - min[0]);
		}
		if (u[0] > max[0]) {
			return abs(u[0] - max[0]);
		}
		if (u[1] < min[1]) {
			return abs(u[1] - min[1]);
		}
		if (u[1] > max[1]) {
			return abs(u[1] - max[1]);
		}
        // inside
        return 0;
	}

	double KDBox::insideDistance(const Vec2 &u) {
		if (!inside(u)) {
			return 0;
		}
		double dx = u[0] < mid(0) ? abs(u[0] - min[0]) : abs(max[0] - u[0]);
		double dy = u[1] < mid(1) ? abs(u[1] - min[1]) : abs(max[1] - u[1]);
		return dx < dy ? dx : dy;
	}

	void KDTree::createKDTree(const vector<Vert *> &verts) {
		box = new KDBox;
		for (int v = 0; v < verts.size(); v++) {
			*box += verts[v]->u;
		}

		vector<Vert *> verts_buffer = verts;
		root = new KDNode(verts_buffer, 0, box, NULL);
	}

	void KDTree::updateKDTree(const Vert *&vert) {

	}

// In argus logic, material positions of vertices won't change. Only adding new vert and removing
// old vert can happen. Also new vert is not going to be outside the entire bounding box.
	bool KDTree::addVert(Vert *vert, KDNode *node) {
		if (node == NULL) {
			return false;
		}
		if (!node->box->inside(vert->u)) {
			return false;
		}
		if (node->isLeaf()) {    // leaf
			Vert *vert2 = node->vert;
			if (vert == vert2) {
				return true;
			}
			// If two nodes have the same value along one axis, then divide by another axis.
			// Loop won't be infinite unless two nodes are on exactly the same position.
			assert(vert->u != vert2->u);
			while (vert->u[node->axis] == vert2->u[node->axis]) {
				node->axis = (node->axis + 1) % dimension;
			}
			node->vert = NULL;
			double division = (vert->u[node->axis] + vert2->u[node->axis]) / 2;
			KDBox *left_box = new KDBox(node->box);
			left_box->max[node->axis] = division;
			KDBox *right_box = new KDBox(node->box);
			right_box->min[node->axis] = division;
			vector<Vert *> left_verts(1), right_verts(1);
			if (vert->u[node->axis] < division) {
				left_verts[0] = vert;
				right_verts[0] = vert2;
			} else {
				left_verts[0] = vert2;
				right_verts[0] = vert;
			}
			node->left = new KDNode(left_verts, (node->axis + 1) % dimension, left_box, node);
			node->right = new KDNode(right_verts, (node->axis + 1) % dimension, right_box, node);
			return true;
		}
		if (node->empty()) {    // empty space
			node->vert = vert;
			vert->kdNode = node;
			return true;
		}
		assert(node->left && node->right);
		// Right child has higher priority, as nodes on division line are to be grouped into right tree.
		if (!addVert(vert, node->right)) {
			return addVert(vert, node->left);
		}
		return true;
	}

	bool KDTree::addVert(Vert *vert) {
		return addVert(vert, root);
	}

	void KDTree::mergeEmptySpace(KDNode *node) {
		KDNode *parent = node->parent;
		KDNode *sibling = node->getSibling();
		if (node->empty() && sibling->empty()) {
			parent->left = parent->right = NULL;
			delete node->box;
			delete node;
			delete sibling->box;
			delete sibling;
			mergeEmptySpace(parent);
		}
	}

	bool KDTree::removeVert(Vert *vert, KDNode *node) {
		if (node == NULL) {
			return false;
		}
		if (!node->box->inside(vert->u)) {
			return false;
		}
		if (node->vert) {
			Vert *vert2 = node->vert;
			if (vert == vert2) {
				KDNode *parent = node->parent;
				if (parent == NULL) {    // root
					cerr << "#Error: the root can't be removed." << endl;
					exit(-1);
				}
				KDNode *sibling = node->getSibling();
				KDNode *grandparent = parent->parent;
				if (grandparent == NULL) {    // parent is root
					root = sibling;
					root->parent = NULL;
					box = root->box;
					delete node->box;
					delete node;
					delete parent->box;
					delete parent;
					return true;
				}
				vert->kdNode = NULL;
				node->vert = NULL;
				// Merge the space if both of them are empty.
				mergeEmptySpace(node);
				return true;
			}
			return false;
		}
		if (node->left == NULL && node->right == NULL) {
			return false;
		}
		assert(node->left && node->right);
		if (!removeVert(vert, node->right)) {
			return removeVert(vert, node->left);
		}
		return true;
	}

	bool KDTree::removeVert(Vert *vert) {
		return removeVert(vert, root);
	}

	void KDBox::drawBox() {
		VisualDebugger *vd = VisualDebugger::getInstance();
		double min0 = min[0];
		double min1 = min[1];
		double max0 = max[0];
		double max1 = max[1];
		Vec3 color = Vec3(1, 0, 0);
		vd->addVisualLine2(Vec2(min0, min1), Vec2(min0, max1), color, 'k');
		vd->addVisualLine2(Vec2(max0, min1), Vec2(max0, max1), color, 'k');
		vd->addVisualLine2(Vec2(min0, min1), Vec2(max0, min1), color, 'k');
		vd->addVisualLine2(Vec2(min0, max1), Vec2(max0, max1), color, 'k');
	}

	void drawDivision(KDNode *node) {
		node->box->drawBox();
	}

	void KDTree::traverse(KDNode *node, TraverseCallback callback) {
		if (node == NULL) {
			return;
		}
		if (node->isLeaf()) {
			callback(node);
			return;
		}
		if (node->empty()) {
			return;
		}
		traverse(node->left, callback);
		traverse(node->right, callback);
	}

	void KDTree::traverseAll() {
		traverse(root, drawDivision);
	}

	unsigned int KDTree::leavesNumber(KDNode *node) {
		if (node == NULL) {
			node = root;
		}
		if (node->isLeaf()) {
			return 1;
		}
		if (node->empty()) {
			return 0;
		}
		return leavesNumber(node->left) + leavesNumber(node->right);
	}

	Vert *KDTree::getVertWithinRadius(const Vec2 &u, const double r, KDNode *node, bool findGlobally) {
		if (node->isLeaf()) {
			Vert *vert = node->vert;
			if (norm(vert->u - u) <= r) {
				return vert;
			}
		} else if (!node->empty()) {
			Vert *vert = NULL;
			if (node->right->box->distance(u) <= r) {
				vert = getVertWithinRadius(u, r, node->right);
				if (vert) {
					return vert;
				}
			}
			if (node->left->box->distance(u) <= r) {
				vert = getVertWithinRadius(u, r, node->left);
				if (vert) {
					return vert;
				}
			}
		}
		if (!findGlobally) {
			return NULL;
		}
		KDNode *curNode = node;
		while (curNode && curNode->parent && curNode->box->insideDistance(u) <= r) {
			KDNode *sibling = curNode->getSibling();
			if (sibling->box->distance(u) <= r) {
				Vert *vert = getVertWithinRadius(u, r, sibling);
				if (vert) {
					return vert;
				}
			}
			curNode = curNode->parent;
		}
		return NULL;
	}

	KDTree::KDTree(const vector<Vert *> &verts) {
		createKDTree(verts);
		traverseAll();
	}

	bool KDNode::empty() {
		return (vert == NULL && left == NULL && right == NULL);
	}

	bool KDNode::isLeaf() {
		return vert;
	}

	KDNode *KDNode::getSibling() {
		if (parent == NULL) {
			return NULL;
		}
		return parent->left == this ? parent->right : parent->left;
	}


	bool compare_u(pair<int, double> i, pair<int, double> j) {
		return i.second < j.second;
	}

	KDNode::KDNode(vector<Vert *> &verts, unsigned int _axis, KDBox *_box, KDNode *_parent) {
		parent = _parent;
		box = _box;
		axis = _axis;
		if (verts.size() == 0) {    // empty space
			return;
		}
		if (verts.size() == 1) {    // leaf node.
			vert = verts[0];
			left = right = NULL;
			vert->kdNode = this;
			return;
		}
		vector<Vert *> left_verts;
		vector<Vert *> right_verts;
		// double division = _box->mid(axis);	// Evenly divide the space by the axis

		// Divide the vertices evenly by the number.
		vector<pair<int, double> > toSort(verts.size());
		for (int v = 0; v < verts.size(); v++) {
			toSort[v] = make_pair(v, verts[v]->u[axis]);
		}
		sort(toSort.begin(), toSort.end(), compare_u);

		int m = (verts.size() + 1) / 2;
		int median = toSort[m].first;
		double division = verts[median]->u[axis];
		KDBox *left_box = new KDBox(box);
		left_box->max[axis] = division;
		KDBox *right_box = new KDBox(box);
		right_box->min[axis] = division;

		left_verts.resize(m);
		right_verts.resize(verts.size() - m);
		for (int v = 0; v < verts.size(); v++) {
			if (v < m) {
				left_verts[v] = verts[toSort[v].first];
			} else {
				right_verts[v - m] = verts[toSort[v].first];
			}
		}

		// for (int v = 0; v < verts.size(); v++) {
		// 	if (verts[v]->u[axis] < division) {
		// 		left_verts.push_back(verts[v]);
		// 	} else {
		// 		right_verts.push_back(verts[v]);
		// 	}
		// }

		// If one side has much more nodes than the other side, then move some to make the tree balanced.
		// As the mesh at the beginning is uniformly distributed, the difference between two sides won't be large.
		// Thus this is kind of more efficient than sorting and obtaining the median.
		// int diff = (left_verts.size() - right_verts.size())/2;
		// if (diff < 0) {
		// 	vector<Vert*> toMove = popSomeMin(right_verts, -diff, axis);
		// 	for (int m = 0; m < toMove.size(); m++) {
		// 		left_verts.push_back(toMove[m]);
		// 	}
		// } else if (diff > 0) {
		// 	vector<Vert*> toMove = popSomeMax(left_verts, diff, axis);
		// 	for (int m = 0; m < toMove.size(); m++) {
		// 		right_verts.push_back(toMove[m]);
		// 	}
		// }
		verts.clear();
		left = new KDNode(left_verts, (axis + 1) % dimension, left_box, this);
		right = new KDNode(right_verts, (axis + 1) % dimension, right_box, this);
	}

	void KDNode::deleteKDNode() {
		if (box) delete box;
		if (left) {
			left->deleteKDNode();
			delete left;
		}
		if (right) {
			right->deleteKDNode();
			delete right;
		}
	}

	void KDTree::deleteKDTree() {
		root->deleteKDNode();
		delete root;
		// delete box;
	}
}
