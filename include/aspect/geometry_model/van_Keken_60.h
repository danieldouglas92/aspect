/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_geometry_model_vanKeken60_h
#define _aspect_geometry_model_vanKeken60_h

#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <aspect/utilities.h>
#include <deal.II/dofs/dof_tools.h>

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A geometry model based on the 2D Cartesian van Keken 2008 subduction
     * benchmark. A custom mesh that is better suited to deal with
     * the unique constraints of this benchmark is necessary to obtain a
     * continuous pressure field.
     *
     * Because this mesh is highly specific, many of these functions are
     * present only because of the way that we need to interface with the
     * GeometryModel class, but they return nothing.
     */

    template <int dim>
    class vanKeken60 : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        vanKeken60();

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        virtual
        Point<dim> get_extents () const;


        virtual
        Point<dim> get_origin () const;
        

        double length_scale () const override;

        /**
         * Return the depth that corresponds to the given
         * position. The documentation of the base class (see
         * GeometryModel::Interface::depth()) describes in detail how
         * "depth" is interpreted in general.
         */
        double depth(const Point<dim> &position) const override;


        double height_above_reference_surface(const Point<dim> &position) const override;

        /**
         * @copydoc Interface<dim>::representative_point()
         */
        Point<dim> representative_point(const double depth) const override;

        /**
         * @copydoc Interface<dim>::maximal_depth()
         */
        double maximal_depth() const override;

        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         */
        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        /**
         * Return symbolic names for all boundary components. Their names are
         * described in the documentation of this plugin, at the bottom of the
         * .cc file.
         */
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        bool
        point_is_in_domain(const Point<dim> &point) const override;

        /**
         * Returns what the natural coordinate system for this geometry model is,
         * which in this case is cartesian.
         */
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        /**
         * Declare the parameters this class takes through input files. However,
         * for this geometry there are no parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &);

        /**
         * Read the parameters this class declares from the parameter file.
         * Again, there are no parameters to be declared for this class.
         */
        void
        parse_parameters (ParameterHandler &) override;
    };
  }
}

#endif
