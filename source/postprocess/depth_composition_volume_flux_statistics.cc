// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/depth_composition_volume_flux_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DepthCompositionVolumeFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      // Determine the units for the compositional volume fluxes based on the output settings.
      const std::string length = dim == 2 ? "m^2" : "m^3";
      const std::string time = (this->convert_output_to_years()) ? "yr" : "s";
      const std::string unit = length + "/" + time;

      const double in_years = (this->convert_output_to_years())
                              ?
                              year_in_seconds
                              :
                              1.0;

      // Determine if the model is spherical or Cartesian. Ellipsoidal geometries are not supported.
      bool is_spherical;
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::spherical)
        is_spherical = true;
      else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian)
        is_spherical = false;
      else
        AssertThrow(false, ExcMessage("The DepthCompositionVolumeFluxStatistics postprocessor only supports Cartesian and spherical geometries."));

      // We need a volume quadrature for the net flux within the cell, and a face quadrature to determine
      // the length/area of the face of the cell that is normal to gravity.
      const Quadrature<dim> &quadrature_formula_volume = this->introspection().quadratures.velocities;
      const Quadrature<dim-1> &quadrature_formula_face = this->introspection().face_quadratures.velocities;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula_volume,
                               update_values            | update_gradients |
                               update_quadrature_points | update_JxW_values);

      FEFaceValues<dim> face_values(this->get_mapping(),
                                    this->get_fe(),
                                    quadrature_formula_face,
                                    update_normal_vectors    |
                                    update_quadrature_points | update_JxW_values);

      // Allocate some storage for volume flux and material inputs.
      std::vector<std::vector<double>> local_volume_fluxes (n_depth_zones, std::vector<double> (this->n_compositional_fields()));
      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());

      // Loop over all cells. If locally owned, first compute the length of cell face that has a normal
      // that is parallel to gravity. Then, construct a bounding box using the current cell. Loop over
      // the depth zones specified in the input file and determine if the depth lies in the bounding box.
      // If it does, compute the average compositional flux within the entire cell, and divide by the cell
      // volume and multiply by the length/area of the face that is normal to gravity.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // The cell is locally owned. Proceed with the rest of the calculation.

            // Determine the length/area of the face of the cell that is normal
            // to the gravity vector. Loop over all faces and if the dot product
            // of the cell normal with gravity is not 0, then that means they are
            // parallel.
            double gravity_normal_face_length = 0.0;
            bool gravity_normal_face_found = false;
            for (const unsigned int f : cell->face_indices())
              {
                face_values.reinit(cell, f);
                for (unsigned int q = 0; q < face_values.n_quadrature_points; ++q)
                  {
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(face_values.quadrature_point(q));
                    if (abs(face_values.normal_vector(q) * gravity) < std::numeric_limits<double>::epsilon())
                      continue;
                    else
                      {
                        gravity_normal_face_length += face_values.JxW(q);
                        gravity_normal_face_found = true;
                      }
                  }
                if (gravity_normal_face_found)
                  break;
              }

            // Now initialize the volume fe values.
            fe_values.reinit (cell);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());

            // Compute cell volume.
            double cell_volume = 0.0;
            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              cell_volume += fe_values.JxW(q);

            // Set the bounding box. We will create a check point that is initially the center
            // of the cell, and then we will shift this point to the desired depth. This depends
            // on the geometry. If spherical, then we need to adjust the radius. If Cartesian, we
            // need to adjust the y/z coordinate depending on the dimension.
            const BoundingBox<dim> cell_bounding_box = cell->bounding_box();
            Point<dim> check_point = cell->center();
            const Point<dim> surface_point = this->get_geometry_model().representative_point(0);

            // Loop over all depth zones, and check if they lie within the current cell.
            for (unsigned int d_idx = 0; d_idx < n_depth_zones; ++d_idx)
              {
                if (is_spherical)
                  {
                    // The model is spherical. First, convert the check point to spherical coordiantes
                    // so we can easily adjust the point to the desired depth. Then, convert back to
                    // Cartesian coordinates to check whether the point lies in the bounding box.
                    const double max_radius = surface_point[0];
                    std::array<double,dim> cell_center_in_spherical = Utilities::Coordinates::cartesian_to_spherical_coordinates(cell->center());

                    cell_center_in_spherical[0] = max_radius - flux_depths[d_idx];
                    Point<dim> shifted_cell_center_cartesian = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(cell_center_in_spherical);
                    check_point                       = shifted_cell_center_cartesian;
                  }
                else
                  check_point[dim-1] = surface_point[dim-1] - flux_depths[d_idx];

                // If the cell does not lie in the bounding box, skip to the next depth zone and repeat.
                if (!cell_bounding_box.point_inside(check_point))
                  continue;

                // Additionally, we need to check whether the point lies ON the bounding box. If it does,
                // then the flux will be double counted because the point will be "within" two cells simultaneously.
                // Divide by two later if this is the case.
                bool on_bounding_box = false;
                if (std::abs(cell_bounding_box.signed_distance(check_point)) < std::numeric_limits<double>::epsilon())
                  on_bounding_box = true;

                // Determine the local contribution of the flux.
                std::vector<double> local_volume_flux_in_gravity(this->n_compositional_fields(), 0.0);
                for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
                  {
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));
                    for (unsigned int c=0; c < this->n_compositional_fields(); ++c)
                      {
                        local_volume_flux_in_gravity[c]
                        += (in.velocity[q] * gravity / gravity.norm())
                           * fe_values.JxW(q)
                           * in.composition[q][c]
                           * in_years
                           * (on_bounding_box ? 0.5 : 1.0);
                      }

                  }

                // Scale the local volume flux by the cell volume and the gravity-normal face length.
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  local_volume_fluxes[d_idx][c] += local_volume_flux_in_gravity[c] / cell_volume * gravity_normal_face_length;
              }
          }

      std::vector<std::vector<double>> global_volume_fluxes (n_depth_zones, std::vector<double>(this->n_compositional_fields(), 0.0));
      {
        // Now collect contributions from all processors and sum them to a global volume flux.
        for (unsigned int d_idx=0; d_idx<n_depth_zones; ++d_idx)
          Utilities::MPI::sum (local_volume_fluxes[d_idx], this->get_mpi_communicator(), global_volume_fluxes[d_idx]);
      }

      // Write to the statistics file.
      for (unsigned int d_idx=0; d_idx<n_depth_zones; ++d_idx)
        {
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              const std::string comp_name = this->introspection().name_for_compositional_index(c);
              const std::string name = comp_name + " vertical volume flux at "
                                       + Utilities::to_string(flux_depths[d_idx])
                                       + " m depth (" + unit + ")";
              statistics.add_value (name, global_volume_fluxes[d_idx][c]);

              // also make sure that the other columns filled by this object
              // all show up with sufficient accuracy and in scientific notation
              statistics.set_precision (name, 8);
              statistics.set_scientific (name, true);
            }
        }

      // Write to the log file.
      std::ostringstream screen_text;
      for (unsigned int d_idx=0; d_idx<n_depth_zones; ++d_idx)
        {
          screen_text << flux_depths[d_idx] << " m: ";
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              screen_text.precision(4);
              screen_text << global_volume_fluxes[d_idx][c]
                          << (c == this->n_compositional_fields()-1 ? "" : " // ");
            }
          screen_text << " " << unit << (d_idx == n_depth_zones-1 ? "" : ",") << std::endl;
        }

      return std::pair<std::string, std::string> ("Volume fluxes through depth: ",
                                                  screen_text.str());
    }



    template <int dim>
    void
    DepthCompositionVolumeFluxStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Depth composition volume flux statistics");
        {
          prm.declare_entry("Depths for flux evaluation",
                            "0.0",
                            Patterns::List(Patterns::Double()),
                            "The depths at which the fluxes are evaluated.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DepthCompositionVolumeFluxStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Depth composition volume flux statistics");
        {
          flux_depths   = Utilities::string_to_double(Utilities::split_string_list(
                                                        prm.get("Depths for flux evaluation")));
          n_depth_zones = flux_depths.size();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DepthCompositionVolumeFluxStatistics,
                                  "depth composition volume flux statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the depth composition volume flux across boundaries. For each provided "
                                  "depth, the depth composition volume flux is computed in the gravity  "
                                  "parallel direction. Therefore, positive flux is in the direction of gravity, "
                                  "and negative flux is opposite gravity. The volume flux is determined with the "
                                  "following expression: "
                                  "1 / V * \\int_{\\Gamma_v} \\mathbf v \\cdot \\mathbf g * C_i * dA"
                                  "Where $V$ is the volume of the cell, $\\Gamma_v$ is an integral over the volume of "
                                  "the cell, $\\mathbf v$ is the velocity, $\\mathbf g$ is the gravity vector, $C_i$ "
                                  "is the i-th compositional field, and dA is the area of the cell face with a face normal "
                                  "perpendicular to the gravity direction."
                                  "\n\n")
  }
}
