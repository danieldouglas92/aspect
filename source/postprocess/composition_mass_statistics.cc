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


#include <aspect/postprocess/composition_mass_statistics.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    CompositionMassStatistics<dim>::execute (TableHandler &statistics)
    {
      if (this->n_compositional_fields() == 0)
        return {"", ""};

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

          std::vector<double> compositional_mass_values(n_q_points);

          MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
          in.requested_properties = MaterialModel::MaterialProperties::density;

          std::vector<double> local_compositional_mass_integrals (this->n_compositional_fields());

          // compute the integral quantities by quadrature
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit (cell);
                in.reinit(fe_values, cell, this->introspection(), this->get_solution());

                this->get_material_model().evaluate(in, out);

                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  {
                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                        compositional_mass_values);
                    for (unsigned int q=0; q<n_q_points; ++q)
                      {
                          const double cutoff_depth = cutoff_depths[depth_index];
                          if (this->get_geometry_model().depth(in.position[q]) >= cutoff_depth)
                            local_compositional_mass_integrals[c] += out.densities[q] * compositional_mass_values[q] * fe_values.JxW(q);
                      }
                  }
              }
          // compute the sum over all processors
          std::vector<double> global_compositional_mass_integrals (local_compositional_mass_integrals.size());
          Utilities::MPI::sum (local_compositional_mass_integrals,
                              this->get_mpi_communicator(),
                              global_compositional_mass_integrals);


          // finally produce something for the statistics file
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              statistics.add_value ("Global mass for composition " + this->introspection().name_for_compositional_index(c) + " below depth " + std::to_string(cutoff_depths[depth_index]),
                                    global_compositional_mass_integrals[c]);
            }

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              const std::string columns[] = {"Global mass for composition " + this->introspection().name_for_compositional_index(c) + " below depth " + std::to_string(cutoff_depths[depth_index])};
              for (const auto &col : columns)
                {
                  statistics.set_precision (col, 8);
                  statistics.set_scientific (col, true);
                }
            }
        }

      return std::pair<std::string, std::string> ("Composition mass:",
                                                  "Done");
    }

    template <int dim>
    void
    CompositionMassStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition mass statistics");
        {
          prm.declare_entry ("Number of depth intervals", "1",
                             Patterns::Integer (1),
                             "The number of depth intervals below which the compositional mass is computed. ");
          prm.declare_entry ("Cutoff depths", "0.0",
                             Patterns::Anything(),
                             "The depths below which the compositional mass is computed. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    CompositionMassStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition mass statistics");
        {
          number_of_depths = prm.get_integer("Number of depth intervals");
          cutoff_depths = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cutoff depths"))),
                                                                  number_of_depths,
                                                                  "Cutoff depths");
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
    ASPECT_REGISTER_POSTPROCESSOR(CompositionMassStatistics,
                                  "composition mass statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
