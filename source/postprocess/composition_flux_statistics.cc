/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/composition_flux_statistics.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    CompositionFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      if (this->n_compositional_fields() == 0)
        return {"", ""};

      const double time_factor = this->convert_output_to_years() ? year_in_seconds : 1.0; 

      // create a quadrature formula based on the compositional element alone.

      for (unsigned int depth_index=0; depth_index<number_of_depths; ++depth_index)
        {
          const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.compositional_field_max;
          const unsigned int n_q_points = quadrature_formula.size();

          FEValues<dim> fe_values (this->get_mapping(),
                                  this->get_fe(),
                                  quadrature_formula,
                                  update_values   |
                                  update_quadrature_points |
                                  update_gradients |
                                  update_JxW_values);

          std::vector<double> compositional_flux_values(n_q_points);

          MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
          in.requested_properties = MaterialModel::MaterialProperties::density;

          std::vector<double> local_compositional_flux_integrals (this->n_compositional_fields());

          // compute the integral quantities by quadrature
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                if (std::abs(cell->center()[dim-1] - cutoff_depths[depth_index]) > depth_layer_half_thickness)
                  continue;
                  
                fe_values.reinit (cell);
                in.reinit(fe_values, cell, this->introspection(), this->get_solution());

                this->get_material_model().evaluate(in, out);

                for (unsigned int q=0; q<n_q_points; ++q)
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    {
                      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(in.position[q]);
                      const Tensor<1,dim> velocity = in.velocity[q];
                      const double normal_velocity = velocity * gravity / gravity.norm();
                      local_compositional_flux_integrals[c] += out.densities[q] * in.composition[q][c] * fe_values.JxW(q) / 
                                                              (2.0 * depth_layer_half_thickness) * normal_velocity * time_factor;
                    }
              }

          // compute the sum over all processors
          std::vector<double> global_compositional_flux_integrals (local_compositional_flux_integrals.size());
          Utilities::MPI::sum (local_compositional_flux_integrals,
                              this->get_mpi_communicator(),
                              global_compositional_flux_integrals);


          // finally produce something for the statistics file
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              statistics.add_value ("Global flux for composition " + this->introspection().name_for_compositional_index(c) + " below depth " + std::to_string(cutoff_depths[depth_index]),
                                    global_compositional_flux_integrals[c]);
            }

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              const std::string columns[] = {"Global flux for composition " + this->introspection().name_for_compositional_index(c) + " below depth " + std::to_string(cutoff_depths[depth_index])};
              for (const auto &col : columns)
                {
                  statistics.set_precision (col, 8);
                  statistics.set_scientific (col, true);
                }
            }
        }

      return std::pair<std::string, std::string> ("","");
    }

    template <int dim>
    void
    CompositionFluxStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition mass flux statistics");
        {
          prm.declare_entry ("Number of depth intervals", "1",
                             Patterns::Integer (1),
                             "The number of depth intervals where the compositional flux is computed. ");
          prm.declare_entry ("Cutoff depths", "0.0",
                             Patterns::Anything(),
                             "The depths at which the compositional flux is computed. ");
          prm.declare_entry ("Depth layer half thickness", "0.0",
                             Patterns::Double (0.0),
                             "The half thickness of the depth layer used to compute the compositional flux. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    CompositionFluxStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition mass flux statistics");
        {
          number_of_depths = prm.get_integer("Number of depth intervals");
          cutoff_depths = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cutoff depths"))),
                                                                  number_of_depths,
                                                                  "Cutoff depths");
          AssertThrow(cutoff_depths.size() == number_of_depths,
                      ExcMessage("The number of cutoff depths must match the number of depth intervals."));

          depth_layer_half_thickness = prm.get_double("Depth layer half thickness");
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
    ASPECT_REGISTER_POSTPROCESSOR(CompositionFluxStatistics,
                                  "composition mass flux statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total flux of this "
                                  "field as defined by the integral "
                                  "$F_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}