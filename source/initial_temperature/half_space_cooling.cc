/*
  Copyright (C) 2016 - 2021 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/half_space_cooling.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    HalfSpaceCooling<dim>::HalfSpaceCooling ()
      :
      surface_boundary_id(numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void
    HalfSpaceCooling<dim>::initialize ()
    {
      // Find the boundary indicator that represents the surface
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);

      // The input ascii table contains one data column (LAB depths(m)) in addition to the coordinate columns.
      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    1);
    }

    template <int dim>
    double
    HalfSpaceCooling<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double seafloor_age = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 0);
      const double T_lithosphere = T_surface + (T_mantle - T_surface) * std::erf(depth / (2 * std::sqrt(kappa * seafloor_age)));
      return T_lithosphere;
    }

    template <int dim>
    void
    HalfSpaceCooling<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/",
                                                          "adiabatic_boundary.txt",
                                                          "Half space cooling");
        prm.enter_subsection("Half space cooling");
        {
          prm.declare_entry ("Mantle temperature", "1673",
                             Patterns::Double (0.),
                             "The value of the mantle temperature. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Surface temperature", "273.15",
                             Patterns::Double (0.),
                             "The value of the surface temperature. Units: \\si{\\kelvin}.");
	  prm.declare_entry ("Lithosphere thickness", "125e3",
			     Patterns::Double (0.),
			     "The thickness of the lithosphere. Units: \\si{\\m}.");
	  prm.declare_entry ("Thermal diffusivity", "0.8e-6",
			     Patterns::Double (0.),
			     "The thermal diffusivity.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    HalfSpaceCooling<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm,"Half space cooling");

        prm.enter_subsection("Half space cooling");
        {
          T_mantle = prm.get_double("Mantle temperature");
          T_surface  = prm.get_double("Surface temperature");
	  lithosphere_thickness = prm.get_double("Lithosphere thickness");
	  kappa = prm.get_double("Thermal diffusivity");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(HalfSpaceCooling,
                                              "half space cooling",
                                              "An initial temperature condition that allows for discretizing "
                                              "the model domain into two layers separated by a user-defined "
                                              "isothermal boundary. The user includes an input ascii data file "
                                              "that is formatted as 3 columns of `longitude(radians)', "
                                              "`colatitude(radians)', and `isotherm depth(meters)', where `isotherm depth' represents the depth "
                                              "of an initial temperature of 1673.15 K (by default). "
                                              "The first lines in the data file may contain any number of comments if they begin "
                                              "with `#', but one of these lines needs to contain the number of grid points "
                                              "in each dimension as for example `# POINTS: 69 121'. Note that the coordinates need "
                                              "to be sorted in a specific order: the `longitude' coordinate needs to ascend first, "
                                              "followed by the `colatitude' coordinate in order to assign the correct data (isotherm depth) to the "
                                              "prescribed coordinates. "
                                              "The temperature is defined from the surface (273.15 K) to the isotherm depth (1673.15 K) "
                                              "as a linear gradient. Below the isotherm depth the temperature increases "
                                              "approximately adiabatically (0.0005 K per meter). "
                                              "This plugin should work for all geometry models, but is currently only tested for spherical models.")
  }
}
