#ifndef VISUALDEBUG_HPP
#define VISUALDEBUG_HPP

#include "vectors.hpp"
#include <vector>
using namespace std;

namespace arcsim {
	const Vec3 DEFAULT_COLOR(0, 0, 0);
	const unsigned char DEFAULT_KEY = '0';
	const Vec3 COLORS[] = {
			Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(0, 1, 1), Vec3(1, 0, 1),
			Vec3(1, 1, 0), Vec3(1, .5, 0), Vec3(.5, 1, 0), Vec3(1, 0, .5), Vec3(1, .5, 0),
			Vec3(0, 1, .5), Vec3(0, .5, 1), Vec3(1, .5, .5), Vec3(.5, 1, .5), Vec3(.5, .5, 1)
	};


	struct VisualInfo {
		unsigned char key;    // Keyboard control
		vector<Vec3> pos; // World positions
		vector<Vec2> u_pos; // Material positions
		vector<pair<Vec3, Vec3> > lines;    // World space lines;
		vector<pair<Vec2, Vec2> > u_lines;    // Material space lines;
		Vec3 color;

		VisualInfo(unsigned char _key, Vec3 _color) : key(_key), color(_color) {}
	};

	struct VisualDebugger {
		static VisualDebugger *getInstance();

		vector<VisualInfo> infoList;
		vector<unsigned char> keys;
		vector<bool> display;

		void clearData();    // clear all data


		// Add a display node in world space
		void addVisualPoint3(Vec3 p);

		void addVisualPoint3(Vec3 p, Vec3 color);

		void addVisualPoint3(Vec3 p, unsigned char key);

		void addVisualPoint3(Vec3 p, Vec3 color, unsigned char key);

		// Add a display node in material space
		void addVisualPoint2(Vec2 p);

		void addVisualPoint2(Vec2 p, Vec3 color);

		void addVisualPoint2(Vec2 p, unsigned char key);

		void addVisualPoint2(Vec2 p, Vec3 color, unsigned char key);

		// Add a display line in world space
		void addVisualLine3(Vec3 p1, Vec3 p2);

		void addVisualLine3(Vec3 p1, Vec3 p2, Vec3 color);

		void addVisualLine3(Vec3 p1, Vec3 p2, unsigned char key);

		void addVisualLine3(Vec3 p1, Vec3 p2, Vec3 color, unsigned char key);

		// Add a display line in material space
		void addVisualLine2(Vec2 p1, Vec2 p2);

		void addVisualLine2(Vec2 p1, Vec2 p2, Vec3 color);

		void addVisualLine2(Vec2 p1, Vec2 p2, unsigned char key);

		void addVisualLine2(Vec2 p1, Vec2 p2, Vec3 color, unsigned char key);

	private:
		static VisualDebugger *instance;

	};

}
#endif