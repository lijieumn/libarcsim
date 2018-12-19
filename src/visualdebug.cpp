#include "visualdebug.hpp"
#include "util.hpp"

namespace arcsim {
	VisualDebugger *VisualDebugger::instance;

	VisualDebugger *VisualDebugger::getInstance() {
		if (instance == NULL) {
			instance = new VisualDebugger;
		}
		return instance;
	}

	void VisualDebugger::clearData() {
		// keys.clear();
		// display.clear();
		// infoList.clear();
		for (int i = 0; i < infoList.size(); i++) {
			infoList[i].pos.clear();
			infoList[i].u_pos.clear();
			infoList[i].lines.clear();
			infoList[i].u_lines.clear();
		}
	}

	void VisualDebugger::addVisualPoint3(Vec3 p) {
		addVisualPoint3(p, DEFAULT_COLOR, DEFAULT_KEY);
	}

	void VisualDebugger::addVisualPoint3(Vec3 p, Vec3 color) {
		addVisualPoint3(p, color, DEFAULT_KEY);
	}


	void VisualDebugger::addVisualPoint3(Vec3 p, unsigned char key) {
		addVisualPoint3(p, DEFAULT_COLOR, key);
	}

	void VisualDebugger::addVisualPoint3(Vec3 p, Vec3 color, unsigned char key) {
		int index = find(key, keys);
		if (index < 0) {
			keys.push_back(key);
			display.push_back(false);
			VisualInfo info(key, color);
			infoList.push_back(info);
			index = infoList.size() - 1;
		}
		infoList[index].pos.push_back(p);
	}


	void VisualDebugger::addVisualPoint2(Vec2 p) {
		addVisualPoint2(p, DEFAULT_COLOR, DEFAULT_KEY);
	}

	void VisualDebugger::addVisualPoint2(Vec2 p, Vec3 color) {
		addVisualPoint2(p, color, DEFAULT_KEY);
	}


	void VisualDebugger::addVisualPoint2(Vec2 p, unsigned char key) {
		addVisualPoint2(p, DEFAULT_COLOR, key);
	}

	void VisualDebugger::addVisualPoint2(Vec2 p, Vec3 color, unsigned char key) {
		int index = find(key, keys);
		if (index < 0) {
			keys.push_back(key);
			display.push_back(false);
			VisualInfo info(key, color);
			infoList.push_back(info);
			index = infoList.size() - 1;
		}
		infoList[index].u_pos.push_back(p);
	}

	void VisualDebugger::addVisualLine3(Vec3 p1, Vec3 p2) {
		addVisualLine3(p1, p2, DEFAULT_COLOR, DEFAULT_KEY);
	}

	void VisualDebugger::addVisualLine3(Vec3 p1, Vec3 p2, Vec3 color) {
		addVisualLine3(p1, p2, color, DEFAULT_KEY);
	}


	void VisualDebugger::addVisualLine3(Vec3 p1, Vec3 p2, unsigned char key) {
		addVisualLine3(p1, p2, DEFAULT_COLOR, key);
	}

	void VisualDebugger::addVisualLine3(Vec3 p1, Vec3 p2, Vec3 color, unsigned char key) {
		int index = find(key, keys);
		if (index < 0) {
			keys.push_back(key);
			display.push_back(false);
			VisualInfo info(key, color);
			infoList.push_back(info);
			index = infoList.size() - 1;
		}
		infoList[index].lines.push_back(make_pair(p1, p2));
	}

	void VisualDebugger::addVisualLine2(Vec2 p1, Vec2 p2) {
		addVisualLine2(p1, p2, DEFAULT_COLOR, DEFAULT_KEY);
	}

	void VisualDebugger::addVisualLine2(Vec2 p1, Vec2 p2, Vec3 color) {
		addVisualLine2(p1, p2, color, DEFAULT_KEY);
	}


	void VisualDebugger::addVisualLine2(Vec2 p1, Vec2 p2, unsigned char key) {
		addVisualLine2(p1, p2, DEFAULT_COLOR, key);
	}

	void VisualDebugger::addVisualLine2(Vec2 p1, Vec2 p2, Vec3 color, unsigned char key) {
		int index = find(key, keys);
		if (index < 0) {
			keys.push_back(key);
			display.push_back(false);
			VisualInfo info(key, color);
			infoList.push_back(info);
			index = infoList.size() - 1;
		}
		infoList[index].u_lines.push_back(make_pair(p1, p2));
	}
}