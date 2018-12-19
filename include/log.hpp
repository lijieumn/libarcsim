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

#ifndef LOG_HPP
#define LOG_HPP

#include <string>
using namespace std;

namespace arcsim {
    class Log {
    private:
        int frame = 0;
        int nodeNumber = 0;
        int contactNumber = 0;
        int contactNumberBefore = 0;
        int intersectionNumber = 0;
        int intersectionNumberAfter = 0;
        double physicsStepTime = 0;
        double remeshingTime = 0;
        double collisionTime = 0;
        double extraTime = 0;
        double mergingTime = 0;
        double solverTime = 0;
        double contactSolveTime = 0;
        double icmTime = 0;
        double error = 0;
        int iterations;

        string logFile;

        static Log *instance;

    public:

        static Log *getInstance();

        void init(string logFile);

        void reset();

        void output();

        void saveTotalTime(double tt);

        void setFrame(int _frame);

        void setNodeNumber(int _nodeNumber);

        void setContactNumber(int _contactNumber);

        void setContactNumberBefore(int _contactNumberBefore);

        void setIntersectionNumber(int _intersectionNumber);

        void setIntersectionNumberAfter(int _intersectionNumberAfter);

        void setPhysicsStepTime(double _physicsStepTime);

        void setRemeshingTime(double _remeshingTime);

        void setCollisionTime(double _collisionTime);

        void setExtraTime(double _extraTime);

        void setMergingTime(double _mergingTime);

        void setSolverTime(double _solverTime);

        void setContactSolveTime(double _contactSolveTime);

        void setICMTime(double _icmTime);

        void setError(double _error);

        void setIterations(int _iterations);


        Log() {};

        ~Log() {};

    };
}
#endif
