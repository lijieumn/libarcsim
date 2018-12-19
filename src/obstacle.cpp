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

#include "io.hpp"
#include "obstacle.hpp"
#include "util.hpp"
#include <cstdio>

using namespace std;

namespace arcsim {

    Mesh &Obstacle::get_mesh() {
        return curr_state_mesh;
    }

    const Mesh &Obstacle::get_mesh() const {
        return curr_state_mesh;
    }

    Mesh &Obstacle::get_mesh(double time) {
        if (time > end_time)
            delete_mesh(curr_state_mesh);
        if (time < start_time || time > end_time)
            return curr_state_mesh;
        if (!activated) {
            delete_mesh(curr_state_mesh);
            curr_state_mesh = deep_copy(base_mesh);
        }

        //    if( frame_sequence.size() > 0 ){
        // 	interpolate_frames(time);
        // }
        if (file_names.size() > 0) {
            interpolate_frames(time);
        }

        if (transform_spline) {
            DTransformation dtrans = get_dtrans(*transform_spline, time);
            Mesh &mesh = curr_state_mesh;
            for (int n = 0; n < curr_state_mesh.nodes.size(); n++) {
                mesh.nodes[n]->x = base_mesh.nodes[n]->x - transform_spline->center;
                mesh.nodes[n]->x = apply_dtrans(dtrans, mesh.nodes[n]->x);
                mesh.nodes[n]->x += transform_spline->center;
            }
            compute_ws_data(mesh);
        }

        if (!activated)
            update_x0(curr_state_mesh);
        activated = true;
        return curr_state_mesh;
    }

    void Obstacle::interpolate_frames(double t) {

        // int currFrame = int(t / frame_time);
        double r = t / (frame_time * per_duration);
        int curIndex = int(r);

        if (curIndex >= file_names.size() - 1) {
            if (prevprev) {
                delete_mesh(*prevprev);
                delete prevprev;
                prevprev = NULL;
            }
            if (prev) {
                delete_mesh(*prev);
                delete prev;
                prev = NULL;
            }
            if (next) {
                delete_mesh(*next);
                delete next;
                next = NULL;
            }
            if (nextnext) {
                delete_mesh(*nextnext);
                delete nextnext;
                nextnext = NULL;
            }
            return;
        }

        // if( currFrame >= (file_names.size() - 1)*per_duration ){
        // 	return;
        // }

        double weight = (r - curIndex);

        // double tol = 1e-8;
        // if (weight < tol) {
        if (curIndex > prevIndex) {
            if (prevprev) {
                delete_mesh(*prevprev);
                delete prevprev;
                prevprev = NULL;
            }
            prevprev = prev;
            prevprevIndex = prevIndex;
            prev = next;
            prevIndex = nextIndex;
            next = nextnext;
            nextIndex = nextnextIndex;
            nextnext = NULL;
            nextnextIndex = -1;
        }

        // double weight = (t / frame_time) - double(currFrame);
        // double weight = double(currFrame%per_duration)/double(per_duration);

        // int curIndex = currFrame/per_duration;
        Mesh &mesh = curr_state_mesh;
        // if (!prev) {
        if (curIndex != prevIndex) {
            if (prev) {
                delete_mesh(*prev);
                delete prev;
                prev = NULL;
            }
            prevIndex = curIndex;
            prev = new Mesh;
            string &prev_file = file_names[prevIndex];
            load_obj(*prev, prev_file);
            apply_transformation(*prev, transform);
        }
        if (curIndex + 1 != nextIndex) {
            if (next) {
                delete_mesh(*next);
                delete next;
                next = NULL;
            }
            nextIndex = curIndex + 1;
            next = new Mesh;
            string &next_file = file_names[nextIndex];
            load_obj(*next, next_file);
            apply_transformation(*next, transform);
        }
        // Mesh &prev = frame_sequence[ currFrame ];
        // Mesh &next = frame_sequence[ currFrame+1 ];
        if (curIndex > 0 && (curIndex - 1 != prevprevIndex)) {
            if (prevprev) {
                delete_mesh(*prevprev);
                delete prevprev;
                prevprev = NULL;
            }
            prevprevIndex = curIndex - 1;
            prevprev = new Mesh;
            load_obj(*prevprev, file_names[prevprevIndex]);
            apply_transformation(*prevprev, transform);
        }
        if (curIndex < file_names.size() - 2 && curIndex + 2 != nextnextIndex) {
            if (nextnext) {
                delete_mesh(*nextnext);
                delete nextnext;
                nextnext = NULL;
            }
            nextnextIndex = curIndex + 2;
            nextnext = new Mesh;
            load_obj(*nextnext, file_names[nextnextIndex]);
            apply_transformation(*nextnext, transform);
        }
        for (int n = 0; n < mesh.nodes.size(); n++) {
            Node *node = mesh.nodes[n];
            Node *prevNode = prev->nodes[n];
            Node *nextNode = next->nodes[n];
            double t0 = curIndex * per_duration * frame_time;
            double t1 = (curIndex + 1) * per_duration * frame_time;
            double s = (t - t0) / (t1 - t0), s2 = s * s, s3 = s2 * s;

            if (curIndex == 0) {
                prevNode->v = Vec3(0);
            } else {
                prevNode->v = (nextNode->x - prevprev->nodes[n]->x) * .5 / (t1 - t0);
            }
            if (curIndex == file_names.size() - 2) {
                nextNode->v = Vec3(0);
            } else {
                nextNode->v = (nextnext->nodes[n]->x - prevNode->x) * .5 / (t1 - t0);
            }
            node->x = prevNode->x * (2 * s3 - 3 * s2 + 1) + nextNode->x * (-2 * s3 + 3 * s2)
                      + (prevNode->v * (s3 - 2 * s2 + s) + nextNode->v * (s3 - s2)) * (t1 - t0);
        }

        // Linear interpolation
        // for(int n = 0; n < mesh.nodes.size(); n++){
        // 	Node *node = mesh.nodes[n];
        // 	Node *prevNode = prev.nodes[n];
        // 	Node *nextNode = next.nodes[n];
        // 	node->x = (1.0-weight)*prevNode->x + (weight)*nextNode->x;
        // }
        compute_ws_data(mesh);


    }

    void Obstacle::blend_with_previous(double t, double dt, double blend) {
        const Motion *spline = transform_spline;
        Transformation trans = (spline)
                               ? get_trans(*spline, t)
                                 * inverse(get_trans(*spline, t - dt))
                               : identity();
        Mesh &mesh = curr_state_mesh;
        for (int n = 0; n < mesh.nodes.size(); n++) {
            Node *node = mesh.nodes[n];
            Vec3 x0 = trans.apply(node->x0);
            node->x = x0 + blend * (node->x - x0);
        }
        compute_ws_data(mesh);
    }
}