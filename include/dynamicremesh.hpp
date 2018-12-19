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

#ifndef DYNAMICREMESH_HPP
#define DYNAMICREMESH_HPP

#include <map>
#include "cloth.hpp"
#include "nearobs.hpp"
#include "nearestplane.hpp"
#include "collision.hpp"
#include "mergehelper.hpp"

namespace arcsim {
// void static_remesh (Cloth &cloth);
    void static_remesh(Cloth &cloth, const vector<AccelStruct *> &obs_accs, bool collapse = true);

    std::vector<ArgusImpact>
    collision_fixed(Cloth &cloth, std::vector<Mesh *> obs_meshes, const std::vector<Plane> &planes,
                    const std::vector<Impact> &impacts, bool plasticity);

    std::vector<ArgusImpact>
    collision_refine(Cloth &cloth, std::vector<Mesh *> obs_meshes, const std::vector<Plane> &planes,
                     const std::vector<Impact> &impacts, bool plasticity);

    void collision_coarsen(Cloth &cloth, std::vector<Mesh *> obs_meshes, const vector<Plane> &planes, bool plasticity);

// void dynamic_remesh (Cloth &cloth, const std::vector<Plane> &planes,
//                      bool plasticity);
    void dynamic_remesh(Cloth &cloth,
                        bool plasticity, const vector<AccelStruct *> &obs_accs,
                        bool collapse = true);

    void collapse_edges(Cloth &cloth, const std::vector<Plane> &planes, const vector<AccelStruct *> &obs_accs);

    void print_impacts(const Impact &impact);

    void print_argus_impacts(ArgusImpact &aImpact);

    void insert_nodes(vector<Cluster *> &clusters, Cloth &cloth, const vector<Mesh *> obs_meshes);

    void insert_nodes(ImpactPoint *point, Cloth &cloth);

    RemeshOp argus_based_flip_edges(vector<Face *> &active, Mesh &mesh, const vector<AccelStruct *> &obs_accs);

    vector<ArgusImpact> convert_to_argus(const vector<NodalImpact *> nImpacts);

    void setPlaneSet(CuttingPlaneSet &planeSet);

    void setPlanes(vector<Plane> &planes);
}
#endif
