#include "collision.hpp"
#include "magic.hpp"

namespace arcsim {
// Bridson-Harmon Collision Response
    bool collision_response_bh(vector<Mesh *> &meshes, const vector<Constraint *> &cons,
                               const vector<Mesh *> &obs_meshes, double dt);

    void collision_response_arcsim(std::vector<Mesh *> &meshes,
                                   const std::vector<Constraint *> &cons,
                                   const std::vector<Mesh *> &obs_meshes);
}