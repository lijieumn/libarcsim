#ifndef KDTREE_HPP
#define KDTREE_HPP

#include "mesh.hpp"

#include <vector>
using namespace std;

namespace arcsim {
// Currently this k-d tree structure only supports 2-d space.
    static const unsigned int dimension = 2;    // dimension of the space.

#define _MAX 0xffffffff
#define _MIN -0xfffffff

    struct KDNode;
    struct Vert;

    typedef void (*TraverseCallback)(KDNode *node);

    struct KDBox {
        double min[dimension];
        double max[dimension];

        void reset();

        KDBox operator+(const KDBox &box);

        KDBox &operator+=(const KDBox &box);

        KDBox &operator+=(const Vec2 &u);

        bool operator==(const KDBox &box);

        double mid(unsigned int axis);

        bool inside(const Vec2 &u);

        KDBox() { reset(); }

        KDBox(KDBox *box) {
            reset();
            (*this) += (*box);
        }

        // Get the distance of a vert from this box. Return 0 if it's inside.
        double distance(const Vec2 &u);

        // Get the distance of a vert from the boundary of the box if the vert is inside.
        double insideDistance(const Vec2 &u);

        // Display the box in material space
        void drawBox();
    };

    struct KDNode {
        unsigned int axis;    // axis used to select the splitting plane
        KDBox *box;    // division point along the corresponding axis.

        Vert *vert = NULL;    // vertex contained. NULL if this kd-node isn't a leaf
        KDNode *left = NULL;    // left child
        KDNode *right = NULL;    // right child
        KDNode *parent = NULL;    // parent

        KDNode(vector<Vert *> &verts, unsigned int _axis, KDBox *_box, KDNode *_parent);

        // Empty if it doesn't have any node inside it
        bool empty();

        // Is it leaf
        bool isLeaf();

        // Get sibling
        KDNode *getSibling();

        void deleteKDNode();
    };

    struct KDTree {
        KDNode *root = NULL;    // root node.

        KDBox *box = NULL;    // The box that contain all vertices.

        // Construct the k-d tree
        void createKDTree(const vector<Vert *> &verts);

        // Update vertex in k-d tree
        void updateKDTree(const Vert *&vert);

        // Add new vertex into k-d tree at certain node
        bool addVert(Vert *vert, KDNode *node);

        // Add new vertex into k-d tree at root
        bool addVert(Vert *vert);

        // Remove vertex from k-d tree starting from a certain node
        bool removeVert(Vert *vert, KDNode *node);

        // Remove vertex from k-d tree starting from the root
        bool removeVert(Vert *vert);

        // Recursively traverse nodes
        void traverse(KDNode *node, TraverseCallback callback);

        // Traverse from the root
        void traverseAll();

        // Merge space if both two (the node and its sibling) are empty.
        void mergeEmptySpace(KDNode *node);

        // Get number of leaves
        unsigned int leavesNumber(KDNode *node = NULL);

        // Find a vert that is withing a certain radius from all descendants of this node.
        // If findGlobally is true, then if it fails to find from descendants, it keeps looking at its sibling and parent.
        Vert *getVertWithinRadius(const Vec2 &u, const double r, KDNode *node, bool findGlobally = false);

        void deleteKDTree();

        KDTree(const vector<Vert *> &verts);
    };
}

#endif
