/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/particle/teleporter/most_compressive_stress.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>
#include <aspect/parameters.h>
#include <fstream>
#include <string>
#include <deal.II/base/mpi.h>

// Maybe need something from this?
#include <aspect/geometry_model/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/prescribed_dilation/interface.h>
#include <aspect/material_model/interface.h>


namespace aspect
{
  namespace Particle
  {
    namespace Teleporter
    {

      template <int dim>
      void
      MostCompressiveStress<dim>::initialize()
      {
        this->get_signals().modify_reaction_outputs.connect(
          [&](const SimulatorAccess<dim> &simulator_access,
              const MaterialModel::MaterialModelInputs<dim> &in,
              MaterialModel::MaterialModelOutputs<dim> &out)
        {
          this->create_particles(simulator_access, in, out);
        }
        );

        this->get_signals().post_compute_reactions.connect(
          [&](const SimulatorAccess<dim> &simulator_access)
        {
          this->advect_particles(simulator_access);
        }
        );

      }



      template <int dim>
      bool MostCompressiveStress<dim>::new_teleportation_step()
      {
        // Here, I need to add a flag that checks to see how close the particle is to the boundary of the domain.
        return false;
      }



      template <int dim>
      void MostCompressiveStress<dim>::local_teleport_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                                           const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                                           const std::vector<Tensor<1,dim>> &vector_field)
      {
        // Implement the local teleportation step here. This is likely going to just replace `do_particle_stuff`.
      }



      template <int dim>
      std::array<std::pair<double, Tensor<1,dim>>, dim>
      MostCompressiveStress<dim>::principal_stresses_and_directions(const MaterialModel::MaterialModelInputs<dim> &in,
                                                                    MaterialModel::MaterialModelOutputs<dim> &out,
                                                                    const SimulatorAccess<dim> &simulator_access,
                                                                    const unsigned int q)
      {
        SymmetricTensor<2,dim> stress;

        const double eta = out.viscosities[q];

        const SymmetricTensor<2, dim> deviatoric_strain_rate = (simulator_access.get_material_model().is_compressible()
                                                                ? in.strain_rate[q] - 1. / 3. * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>()
                                                                : in.strain_rate[q]);

        // Add elastic stresses if existent
        if (simulator_access.get_parameters().enable_elasticity == false)
          {
            // Compressive stress is positive in geoscience applications
            stress = -2. * eta * deviatoric_strain_rate;
          }
        else
          {
            // Get the total deviatoric stress from the material model.
            const std::shared_ptr<MaterialModel::ElasticAdditionalOutputs<dim>> elastic_additional_out
              = out.template get_additional_output_object<MaterialModel::ElasticAdditionalOutputs<dim>>();

            Assert(elastic_additional_out != nullptr, ExcMessage("Elastic Additional Outputs are needed for the 'principal stress' postprocessor, but they have not been created."));

            stress = -(elastic_additional_out->deviatoric_stress[q]);
          }

        const std::array<std::pair<double, Tensor<1,dim>>, dim> principal_stresses_and_directions = eigenvectors(stress);

        return principal_stresses_and_directions;
      }



      template <int dim>
      void MostCompressiveStress<dim>::create_particles (const SimulatorAccess<dim> &simulator_access,
                                                         const MaterialModel::MaterialModelInputs<dim> &in,
                                                         MaterialModel::MaterialModelOutputs<dim> &out)
      {
        // This function is called during the loop over all cells in the compute_reactions step just
        // after the material model has been evaluated. Therefore, we do not need to loop over all cells,
        // but instead can just check if the current cell is valid and locally owned. If it is, then we
        // compute the melt volume in that cell, and if it satisfies the criteria for generating a
        // particle, then we assign all that volume to a single particle, and modify the reaction rates
        // to remove all of the melt from the cell.

        // I THINK I HAVE AN IDEA FOR HOW TO ACHIEVE EVERYTHING THAT I NEED TO WITHOUT DOING TOO MUCH MODIFICATION
        // WITHIN THIS LOOP, SINCE THIS IS THE ONLY WAY THAT WE CAN DIRECTLY MODIFY THE MATERIAL MODEL INPUTS REACTIONS
        // WE NEED TO BOTH REMOVE POROSITY FROM DEPTH AND ADD POROSITY TO THE SURFACE. THE TRICK I THINK IS TO NOT DELETE
        // PARTICLES IN THE ADVECTION FUNCTION BELOW. WHAT WE NEED TO DO INSTEAD IS ITERATE OVER EACH PARTICLE IN THE
        // ADVECTION FUNCTION AND EITHER:
        //    1. STOP THE ADVECTION AS SOON AS THE PARTICLE REACHES A CELL THAT IS ON THE MODEL BOUNDARY CONTAINING
        //       A FACE THAT LIES ON THE TOP BOUNDARY, OR
        //    2. CONTINUE ADVECTION UNTIL THE NUMBER OF JUMPS IS REACHED. IF THE PARTICLES DOES NOT END UP IN A CELL
        //       ON THE TOP BOUNDARY, THEN IT SHOULD BE REMOVED.
        // BECAUSE THIS WILL PREVENT PARTICLES FROM BEING REMOVED, THEN IN THIS FUNCTION, WE CAN CHECK EACH CELL THAT
        // CONTAINS A PARTICLE. IF IT DOES, WE KNOW THAT THIS CELL NEEDS TO HAVE THE MELT ON THAT PARTICLE ASSIGNED
        // TO THE REACTION RATES IN THE CURRENT CELL ON THE SURFACE NODES. THEN, WE JUST NEED TO MAKE SURE WE LOOP OVER
        // THESE PARTICLES FIRST, AND DELETE THEM BEFORE WE MAKE ANY NEW PARTICLES.

        // If the current cell being processed in compute_reactions is not valid or locally owned, immediate skip.
        if (in.current_cell.state() != IteratorState::valid || !in.current_cell->is_locally_owned())
          return;

        // Mesh deformation replaces the simulator mapping in setup_dofs(). Initialize
        // this handler on first use so it observes the final mapping.
        if (!teleport_particle_handler_initialized)
          {
            teleport_particle_handler.initialize(simulator_access.get_triangulation(),
                                                 simulator_access.get_mapping(),
                                                 1);
            teleport_particle_handler_initialized = true;
          }

        // The current cell is valid and locally owned, so we can proceed with computing the melt volume.
        const Quadrature<dim> &quadrature_formula = simulator_access.introspection().quadratures.compositional_field_max;
        const unsigned int n_q_points = quadrature_formula.size();
        const Mapping<dim> &mapping = simulator_access.get_mapping();
        FEValues<dim> fe_values ( mapping,
                                  simulator_access.get_fe(),
                                  quadrature_formula,
                                  update_values |
                                  update_gradients |
                                  update_quadrature_points |
                                  update_JxW_values);
        fe_values.reinit (in.current_cell);

        // Get the reaction rate outputs so that we can modify it's entries if we need to.
        const std::shared_ptr<MaterialModel::ReactionRateOutputs<dim>> reaction_rate_outputs
          = out.template get_additional_output_object<MaterialModel::ReactionRateOutputs<dim>>();
        const unsigned int porosity_idx = simulator_access.introspection().compositional_index_for_name("porosity");

        std::vector<double> porosity_values(n_q_points);
        const unsigned int surf_lava_idx = simulator_access.introspection().compositional_index_for_name("surface_lava_volume");

        // If the local melt volume exceeds the user-specified threshold, add a single particle at the center of the cell
        // that contains all of the melt within the cell.
        const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                               std::max(this->get_parameters().reaction_steps_per_advection_step,1U));
        const double reaction_time_step_size = this->get_timestep() / static_cast<double>(number_of_reaction_steps);

        // Zero out all of the surface porosity from the previous advection step. I think that this still
        // needs to be done explicitly before applying the melt from any particles in the currente cell as
        // a reaction rate.
        for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
          reaction_rate_outputs->reaction_rates[q][surf_lava_idx] = -in.composition[q][surf_lava_idx] / reaction_time_step_size;

        // Now check to see if any particles ended up in the current cell from the advection of the particles.
        if (teleport_particle_handler.n_particles_in_cell(in.current_cell) > 0)
          {
            for (auto &particle: teleport_particle_handler.particles_in_cell(in.current_cell))
              {
                const double particle_melt = particle.get_properties()[0];
                std::vector<double> distance_weights(in.n_evaluation_points());
                for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
                  {
                    const Point<dim> &quadrature_point = in.position[q];
                    distance_weights[q] = quadrature_point.distance(particle.get_location());
                  }

                // Normalize the distance weights so that they sum to 1.
                const double total_distance = std::accumulate(distance_weights.begin(), distance_weights.end(), 0.0);
                if (total_distance > 0)
                  for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
                    distance_weights[q] /= total_distance;

                // Apply the particle melt to the reaction rates based on the distance weights.
                for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
                  reaction_rate_outputs->reaction_rates[q][surf_lava_idx] += particle_melt * distance_weights[q] / reaction_time_step_size;
              }

            while (teleport_particle_handler.n_particles_in_cell(in.current_cell) > 0)
              teleport_particle_handler.remove_particle(
                teleport_particle_handler.particles_in_cell(in.current_cell).begin());
          }


        // Integrate the total melt volume within the cell.
        double local_melt_volume = 0;
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            const Point<dim> &quadrature_point = fe_values.quadrature_point(q);
            fe_values[this->introspection().extractors.compositional_fields[porosity_idx]]
            .get_function_values(this->get_solution(),
                                 porosity_values);
            if (this->get_geometry_model().depth(quadrature_point) < dike_extraction_depth)
              local_melt_volume += porosity_values[q] * fe_values.JxW(q);
          }

        if (local_melt_volume > melt_volume_threshold)
          {
            dealii::Particles::Particle<dim> particle;
            particle.set_location(in.current_cell->center());
            particle.set_id(next_particle_id++);

            typename Particle::ParticleHandler<dim>::particle_iterator it =
              teleport_particle_handler.insert_particle(particle, in.current_cell);

            const std::vector<double> stored_melt = {local_melt_volume};
            it->set_properties(stored_melt);

            // Create a fresh ASCII file for this particle and record its initial location.
            // DELETE THIS FILE STUFF ONCE DONE WITH TESTING
            std::ofstream location_file("particle_paths/particle_id_" +
                                        std::to_string(it->get_id()) + ".txt");

            location_file << "# step";
            for (unsigned int d = 0; d < dim; ++d)
              location_file << " x" << d;
            location_file << '\n' << 0;
            for (unsigned int d = 0; d < dim; ++d)
              location_file << ' ' << it->get_location()[d];
            location_file << '\n';

            // Remove the melt from the cell by setting the porosity reaction term
            // Also remove any melt stored as a surface composition.
            for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
              {
                if (this->get_geometry_model().depth(in.position[q]) < dike_extraction_depth)
                  reaction_rate_outputs->reaction_rates[q][porosity_idx] = -in.composition[q][porosity_idx] / reaction_time_step_size;
              }
          }
      }



      template <int dim>
      void MostCompressiveStress<dim>::advect_particles (const SimulatorAccess<dim> &simulator_access)
      {
        // This function takes all of the particles created in the create_particles function and then
        // advects them along the most compressive stress direction. It does this instantaneously, until
        // the particle either exits the model domain, or reaches a maximum of number of allowed particle
        // "jumps".
        std::cout << "Advecting particles..." << std::endl;
        if (!teleport_particle_handler_initialized)
          {
            teleport_particle_handler.initialize(simulator_access.get_triangulation(),
                                                 simulator_access.get_mapping(),
                                                 1);
            teleport_particle_handler_initialized = true;
          }

        const Mapping<dim> &mapping = simulator_access.get_mapping();
        const Quadrature<dim> &quadrature_formula = simulator_access.introspection().quadratures.temperature;
        const unsigned int n_q_points = quadrature_formula.size();
        FEValues<dim> fe_values ( mapping,
                                  simulator_access.get_fe(),
                                  quadrature_formula,
                                  update_values |
                                  update_gradients |
                                  update_quadrature_points |
                                  update_JxW_values);
        MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, simulator_access.n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, simulator_access.n_compositional_fields());

        in.requested_properties = MaterialModel::MaterialProperties::viscosity | MaterialModel::MaterialProperties::additional_outputs;
        const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

        // Loop over the allowed number of jumps for particle advection,
        // then loop over all cells, then loop over all particles, then
        // loop over all quadrature points. This is obviously expensive,
        // but there are several checks that will reduce the actual number
        // of iterations that happen.
        for (unsigned int jump=0; jump<allowed_number_of_jumps; ++jump)
          {
            for (auto &cell: simulator_access.get_dof_handler().active_cell_iterators())
              {
                bool is_top_boundary_cell = false;
                if (cell->at_boundary())
                  {
                    for (const unsigned int f : cell->face_indices())
                      if (cell->face(f)->at_boundary() &&
                          cell->face(f)->boundary_id() == top_boundary_id)
                        {
                          is_top_boundary_cell = true;
                          break;
                        }
                  }

                if (is_top_boundary_cell ||
                    !cell->is_locally_owned() ||
                    teleport_particle_handler.n_particles_in_cell(cell) == 0)
                  continue;

                fe_values.reinit(cell);
                in.reinit(fe_values, cell, simulator_access.introspection(), simulator_access.get_solution());
                simulator_access.get_material_model().create_additional_named_outputs(out);
                simulator_access.get_material_model().evaluate(in, out);

                std::vector<std::array<std::pair<double, Tensor<1,dim>>, dim>> principal_stresses_and_directions(n_q_points);
                for (unsigned int q = 0; q < n_q_points; ++q)
                  principal_stresses_and_directions[q] = this->principal_stresses_and_directions(in,
                                                                                                 out,
                                                                                                 simulator_access,
                                                                                                 q);

                for (auto &particle : teleport_particle_handler.particles_in_cell(cell))
                  {
                    std::ofstream location_file("particle_paths/particle_id_" +
                                                std::to_string(particle.get_id()) + ".txt",
                                                std::ios::app);

                    unsigned int closest_q = numbers::invalid_unsigned_int;
                    double min_distance_squared = std::numeric_limits<double>::max();
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      {
                        const double distance_squared =
                          (particle.get_location() - fe_values.quadrature_point(q)).norm_square();
                        if (distance_squared < min_distance_squared)
                          {
                            min_distance_squared = distance_squared;
                            closest_q = q;
                          }
                      }

                    Tensor<1,dim> move_vector = principal_stresses_and_directions[closest_q][0].second;

                    // Choose the direction opposite gravity, then move half a cell diameter.
                    const Tensor<1,dim> gravity_vector = this->get_gravity_model().gravity_vector(in.position[closest_q]);
                    const double scale = move_vector * gravity_vector < 0
                                         ? cell->diameter() / 2.0
                                         : -cell->diameter() / 2.0;
                    move_vector *= scale;

                    const Point<dim> new_location = particle.get_location() + move_vector;
                    particle.set_location(new_location);

                    location_file << jump;
                    for (unsigned int d = 0; d < dim; ++d)
                      location_file << ' ' << new_location[d];
                    location_file << '\n';
                  }
              }
            teleport_particle_handler.sort_particles_into_subdomains_and_cells();
          }

        // I still need to add the condition that checks to see if any particles
        // remain AFTER the jump loop and are NOT within a boundary cell containing
        // a face located on the top boundary. If any of these particles exist, then
        // they need to be removed.

        // Clean up, destroy particles
        // teleport_particle_handler.clear_particles();
        return;
      }



      template <int dim>
      void MostCompressiveStress<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection ("Particle teleportation");
        {
          prm.declare_entry ("Number of jumps", "100",
                             Patterns::Integer (),
                             "Defines how many jumps a particle should do");
          prm.enter_subsection ("Most compressive stress");
          {
            prm.declare_entry ("Melt volume threshold", "1e50",
                               Patterns::Double (),
                               "Defines the melt volume threshold for teleportation");
            prm.declare_entry ("Dike extraction depth", "0.0",
                               Patterns::Double (),
                               "Defines the dike extraction depth for teleportation");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void MostCompressiveStress<dim>::parse_parameters(ParameterHandler &prm)
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("The MostCompressiveStress teleporter only works if "
                               "there is a compositional field called porosity."));

        AssertThrow(this->introspection().compositional_name_exists("surface_lava_volume"),
                    ExcMessage("The MostCompressiveStress teleporter only works if "
                               "there is a compositional field called surface_lava_volume."));

        const unsigned int surface_lava_volume_idx = this->introspection().compositional_index_for_name("surface_lava_volume");
        AssertThrow(this->introspection().get_composition_descriptions()[surface_lava_volume_idx].type == aspect::CompositionalFieldDescription::generic,
                    ExcMessage("The MostCompressiveStress teleporter requires the surface_lava_volume field to use the generic advection method."));

        AssertThrow(this->get_parameters().use_operator_splitting,
                    ExcMessage("The MostCompressiveStress teleporter must be used with operator splitting."));

        prm.enter_subsection ("Particle teleportation");
        {
          allowed_number_of_jumps = prm.get_integer ("Number of jumps");
          prm.enter_subsection ("Most compressive stress");
          {
            melt_volume_threshold = prm.get_double ("Melt volume threshold");
            dike_extraction_depth = prm.get_double ("Dike extraction depth");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Teleporter
    {
      ASPECT_REGISTER_PARTICLE_TELEPORTER(MostCompressiveStress,
                                          "most compressive stress",
                                          "A plugin in which the particle property is "
                                          "defined as the evolving chemical composition "
                                          "that results from the formation of oceanic crust "
                                          "and lithosphere as mantle material approaches the "
                                          "surface and melts. Note that this does not necessarily "
                                          "conserves the bulk chemical composition of the mantle, "
                                          "since the conversion only depends on the mantle flow. "
                                          "See the crust and lithosphere formation reaction model "
                                          "documentation for more detailed information.")

    }
  }
}
