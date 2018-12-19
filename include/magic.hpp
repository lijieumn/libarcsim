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

#ifndef MAGIC_HPP
#define MAGIC_HPP

#include <string>
// Magic numbers and other hacks

namespace arcsim {
	struct Magic {
		bool fixed_high_res_mesh;
		bool proximity_collisions;
		double handle_stiffness, collision_stiffness;
		double repulsion_thickness, projection_thickness, detection_thickness;
		double proximity_min;
		double friction_coeff;
		double edge_flip_threshold;
		double rib_stiffening;
		double merge_radius;
		double sharp_point_angle;
		double false_proximity_threshold;
		double duplication_threshold;
		double init_frames;
		double argus_tolerance;
		double dt;
		bool combine_tensors;
		bool preserve_creases;
		bool self_contact;
		bool conserve_momentum;
		bool relax_initial_state;
		bool merge_proximal_impacts;
		bool face_edge_constraints;
		bool consider_contact_force;
		bool inelastic_projection;
		bool harmon_inelastic_projection;
		bool use_representative_triangles;
		bool use_stack_overloop_face_check;
		bool use_eigen_solver;
		bool use_proximity_metric;
		bool ccd_forward;
		bool forward_proximity;
		bool apply_adhesion;
		bool facet_solver;
		int node_lifetime;
		bool vfImpacts;
		bool fvImpacts;
		bool eeImpacts;
		bool veImpacts;
		bool evImpacts;
		bool vvImpacts;
		bool enable_localopt;
		std::string sim_type;

		Magic() :
				fixed_high_res_mesh(false),
				proximity_collisions(false),
				handle_stiffness(1e3),
				collision_stiffness(1e9),
				repulsion_thickness(1e-3),
				projection_thickness(2e-3),
				detection_thickness(3e-3),
				proximity_min(4e-3),
				friction_coeff(0.3),
				edge_flip_threshold(0),
				rib_stiffening(1),
				merge_radius(1e-2),
				false_proximity_threshold(1e-1),
				sharp_point_angle(30.0),
				duplication_threshold(0.17),
				init_frames(0),
				argus_tolerance(1e-8),
				combine_tensors(true),
				preserve_creases(false),
				self_contact(true),
				conserve_momentum(false),
				relax_initial_state(true),
				merge_proximal_impacts(true),
				face_edge_constraints(false),
				consider_contact_force(false),
				node_lifetime(-1),
				vfImpacts(true),
				fvImpacts(true),
				eeImpacts(true),
				veImpacts(false),
				evImpacts(false),
				vvImpacts(false),
				inelastic_projection(false),
				harmon_inelastic_projection(false),
				use_representative_triangles(true),
				use_stack_overloop_face_check(false),
				use_eigen_solver(false),
				use_proximity_metric(true),
				ccd_forward(true),
				forward_proximity(false),
				apply_adhesion(false),
				facet_solver(false),
				enable_localopt(false),
				sim_type("argus") {}
	};

	extern Magic magic;
}
#endif
