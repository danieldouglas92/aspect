/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      DislocationCreepParameters::DislocationCreepParameters()
        :  prefactor (numbers::signaling_nan<double>()),
           activation_energy (numbers::signaling_nan<double>()),
           activation_volume (numbers::signaling_nan<double>()),
           stress_exponent (numbers::signaling_nan<double>())
      {}



      template <int dim>
      DislocationCreep<dim>::DislocationCreep ()
        = default;



      template <int dim>
      const DislocationCreepParameters
      DislocationCreep<dim>::compute_creep_parameters (const unsigned int composition,
                                                       const std::vector<double> &phase_function_values,
                                                       const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        DislocationCreepParameters creep_parameters;
        if (phase_function_values == std::vector<double>())
          {
            // no phases
            creep_parameters.prefactor = prefactors[composition];
            creep_parameters.activation_energy = activation_energies[composition];
            creep_parameters.activation_volume = activation_volumes[composition];
            creep_parameters.stress_exponent = stress_exponents[composition];
          }
        else
          {
            // Average among phases
            creep_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                         prefactors, composition,
                                         MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            creep_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_energies, composition);
            creep_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_volumes , composition);
            creep_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                               stress_exponents, composition);
          }

        return creep_parameters;
      }



      template <int dim>
      double
      DislocationCreep<dim>::compute_viscosity (const double strain_rate,
                                                const double pressure,
                                                const double temperature,
                                                const unsigned int composition,
                                                const std::vector<double> &phase_function_values,
                                                const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        const DislocationCreepParameters p = compute_creep_parameters(composition,
                                                                      phase_function_values,
                                                                      n_phase_transitions_per_composition);

        // Power law creep equation:
        //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * std::exp((E + P*V)/(nRT))
        // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
        // E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        double viscosity = 0.5 * std::pow(p.prefactor,-1/p.stress_exponent) *
                           std::exp((p.activation_energy + pressure*p.activation_volume)/
                                    (constants::gas_constant*temperature*p.stress_exponent)) *
                           std::pow(strain_rate,((1. - p.stress_exponent)/p.stress_exponent));

        Assert (viscosity > 0.0,
                ExcMessage ("Negative dislocation viscosity detected. This is unphysical and should not happen. "
                            "Check for negative parameters. Temperature and pressure are "
                            + Utilities::to_string(temperature) + " K, " + Utilities::to_string(pressure) + " Pa. "));

        // Creep viscosities become extremely large at low
        // temperatures and can therefore provoke floating-point overflow errors. In
        // real rocks, other deformation mechanisms become dominant at low temperatures,
        // so these high viscosities are never achieved. It is therefore both reasonable
        // and desirable to require the single-mechanism viscosity to be smaller than
        // std::sqrt(max_double).
        viscosity = std::min(viscosity, std::sqrt(std::numeric_limits<double>::max()));

        return viscosity;
      }



      template <int dim>
      std::pair<double, double>
      DislocationCreep<dim>::compute_strain_rate_and_derivative (const double stress,
                                                                 const double pressure,
                                                                 const double temperature,
                                                                 const DislocationCreepParameters creep_parameters) const
      {
        // Power law creep equation:
        //   edot_ii_partial = A * stress^n * exp(-(E + P*V)/(RT))
        //   d(edot_ii_partial)/d(stress) = A * n * stress^(n-1) * exp(-(E + P*V)/(RT))
        // A: prefactor, edot_ii_partial: square root of second invariant of deviatoric strain rate tensor attributable to the creep mechanism,
        // stress: deviatoric stress, E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double strain_rate = creep_parameters.prefactor *
                                   std::pow(stress,creep_parameters.stress_exponent) *
                                   std::exp(-(creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                            (constants::gas_constant*temperature));

        const double dstrain_rate_dstress = creep_parameters.prefactor *
                                            creep_parameters.stress_exponent *
                                            std::pow(stress,creep_parameters.stress_exponent-1.) *
                                            std::exp(-(creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                                     (constants::gas_constant*temperature));

        return std::make_pair(strain_rate, dstrain_rate_dstress);
      }


      template <int dim>
      std::pair<double, double>
      DislocationCreep<dim>::compute_log_strain_rate_and_derivative (const double log_stress,
                                                                     const double pressure,
                                                                     const double temperature,
                                                                     const DislocationCreepParameters creep_parameters) const
      {
        // Power law creep equation
        // log(edot_ii_partial) = std::log(A) + n*std::log(stress) - m*std::log(d) - (E + P*V)/(RT)
        //   d(log_edot_ii_partial)/d(log_stress) = n
        // A: prefactor, edot_ii_partial: square root of second invariant of deviatoric strain rate tensor attributable to the creep mechanism,
        // stress: deviatoric stress, E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double log_strain_rate = std::log(creep_parameters.prefactor) +
                                       creep_parameters.stress_exponent * log_stress -
                                       (creep_parameters.activation_energy + pressure*creep_parameters.activation_volume)/
                                       (constants::gas_constant*temperature);

        const double dlog_strain_rate_dlog_stress = creep_parameters.stress_exponent;

        return std::make_pair(log_strain_rate, dlog_strain_rate_dlog_stress);
      }



      template <int dim>
      void
      DislocationCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}$^{-n_{\\text{dislocation}}}$ \\si{\\per\\second}.");
        prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                           Patterns::Anything(),
                           "List of stress exponents, $n_{\\text{dislocation}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                           Patterns::Anything(),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                           Patterns::Anything(),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
      }



      template <int dim>
      void
      DislocationCreep<dim>::parse_parameters (ParameterHandler &prm,
                                               const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        // Make options file for parsing maps to double arrays
        Utilities::MapParsing::Options options(chemical_field_names, "Prefactors for dislocation creep");
        options.list_of_allowed_keys = compositional_field_names;
        options.allow_multiple_values_per_key = true;
        if (expected_n_phases_per_composition)
          {
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }

        // Read parameters, each of size of number of composition + number of phases + 1
        prefactors = Utilities::MapParsing::parse_map_to_double_array(prm.get("Prefactors for dislocation creep"),
                                                                      options);

        options.property_name = "Stress exponents for dislocation creep";
        stress_exponents = Utilities::MapParsing::parse_map_to_double_array(prm.get("Stress exponents for dislocation creep"),
                                                                            options);

        options.property_name = "Activation energies for dislocation creep";
        activation_energies = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation energies for dislocation creep"),
                                                                               options);

        options.property_name = "Activation volumes for dislocation creep";
        activation_volumes = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation volumes for dislocation creep"),
                                                                              options);

        // Check that there are no prefactor entries set to zero,
        // for example because the entry is for a field
        // that is masked anyway, like strain. Despite
        // these compositions being masked, their viscosities
        // are computed anyway and this will lead to division by zero.
        for (const double prefactor : prefactors)
          AssertThrow(prefactor > 0., ExcMessage("The dislocation prefactor should be larger than zero."));
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
#define INSTANTIATE(dim) \
  template class DislocationCreep<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
