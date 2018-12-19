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

#include "log.hpp"
#include <fstream>
#include "util.hpp"

namespace arcsim {
    Log *Log::instance;

    Log *Log::getInstance() {
      if (instance == NULL) {
        instance = new Log;
      }
      return instance;
    }

    static string keys[] = {
            "frame", "nodeNumber", "contactNumber", "intersectionNumber", "intersectionNumberAfter", "physicsStepTime",
            "remeshingTime", "collisionTime", "extraTime", "mergingTime", "solverTime",
            "icmTime", "error", "iterations"
    };

    void Log::init(string logFile) {
      reset();
      this->logFile = logFile;
      // ofstream out(logFile, ios::app);
      // int n = sizeof(keys)/sizeof(string);
      // out << "%";
      // for (int k = 0; k < n; k++) {
      //   out << keys[k] << "\t";
      // }
      // out << endl;
      // out.close();
    }

    void Log::reset() {
      frame = 0;
      nodeNumber = 0;
      contactNumber = 0;
      contactNumberBefore = 0;
      intersectionNumber = 0;
      intersectionNumberAfter = 0;
      physicsStepTime = 0;
      remeshingTime = 0;
      collisionTime = 0;
      extraTime = 0;
      mergingTime = 0;
      solverTime = 0;
      contactSolveTime = 0;
      icmTime = 0;
      error = 0;
      iterations = 0;
    }

    void Log::output() {
      ofstream out(logFile, ios::app);
      // int n = sizeof(keys)/sizeof(string);
      out << frame << "\t";
      out << nodeNumber << "\t";
      out << contactNumber << "\t";
      out << contactNumberBefore << "\t";
      out << intersectionNumber << "\t";
      out << intersectionNumberAfter << "\t";
      out << physicsStepTime << "\t";
      out << remeshingTime << "\t";
      out << collisionTime << "\t";
      out << extraTime << "\t";
      out << mergingTime << "\t";
      out << solverTime << "\t";
      out << contactSolveTime << "\t";
      out << icmTime << "\t";
      out << error << "\t";
      out << iterations << "\t";
      out << endl;
      out.close();
    }



    void Log::saveTotalTime(double tt) {
      string totalTimeFile = stringf("%s_%s", logFile.c_str(), "total_time.txt").c_str();
      ofstream out(totalTimeFile, ios::app);
      out << tt << endl;
    }

    void Log::setFrame(int _frame) {
      frame = _frame;
    }

    void Log::setNodeNumber(int _nodeNumber) {
      nodeNumber = _nodeNumber;
    }

    void Log::setContactNumber(int _contactNumber) {
      contactNumber = _contactNumber;
    }

    void Log::setContactNumberBefore(int _contactNumberBefore) {
      contactNumberBefore = _contactNumberBefore;
    }

    void Log::setIntersectionNumber(int _intersectionNumber) {
      intersectionNumber = _intersectionNumber;
    }

    void Log::setIntersectionNumberAfter(int _intersectionNumberAfter) {
      intersectionNumberAfter = _intersectionNumberAfter;
    }

    void Log::setPhysicsStepTime(double _physicsStepTime) {
      physicsStepTime = _physicsStepTime;
    }

    void Log::setRemeshingTime(double _remeshingTime) {
      remeshingTime = _remeshingTime;
    }

    void Log::setCollisionTime(double _collisionTime) {
      collisionTime = _collisionTime;
    }

    void Log::setExtraTime(double _extraTime) {
      extraTime = _extraTime;
    }

    void Log::setMergingTime(double _mergingTime) {
      mergingTime = _mergingTime;
    }

    void Log::setSolverTime(double _solverTime) {
      solverTime = _solverTime;
    }

    void Log::setContactSolveTime(double _contactSolveTime) {
      contactSolveTime = _contactSolveTime;
    }

    void Log::setICMTime(double _icmTime) {
      icmTime = _icmTime;
    }

    void Log::setError(double _error) {
      error = _error;
    }

    void Log::setIterations(int _iterations) {
      iterations = _iterations;
    }
}