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


#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <aspect/geometry_model/van_Keken_45deep.h>
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
    
    template <int dim>
    vanKeken45Deep<dim>::vanKeken45Deep()
    {}

    template <>
    void
    vanKeken45Deep<2>::
    create_coarse_mesh (parallel::distributed::Triangulation<2> &coarse_grid) const
    {
      // };
      const std::vector<Point<2>> vertices =
      {
        //   0            1             2             3             4
        {0.0, 0.0}, {80e3, 0.0}, {330e3, 0.0}, {600e3, 0.0}, {660e3, 0.0},
        //   5              6               7               8               9
        {0.0, 300e3}, {80e3, 300e3}, {233e3, 183e3}, {330e3, 270e3}, {660e3, 270e3},
        //    10             11            12              13
        {80e3, 520e3}, {660e3, 520e3}, {0.0, 600e3}, {660e3, 600e3}
      };
      // Next, we have to define the cells and the vertices they contain.
      const std::vector<std::array<int, GeometryInfo<2>::vertices_per_cell>>
      cell_vertices =
      {
        {{0, 1, 5, 6}},
        {{8, 9, 10, 11}},
        {{10, 11, 12, 13}},
        {{1, 2, 6, 7}},
        {{2, 3, 7, 8}},
        {{3, 4, 8, 9}},
        {{5, 6, 12, 10}},
        {{6, 7, 10, 8}}
      };

      std::vector<CellData<2>> cells(cell_vertices.size(), CellData<2>());
      for (unsigned int i = 0; i < cell_vertices.size(); ++i)
        {
          for (unsigned int j = 0; j < cell_vertices[i].size(); ++j)
            cells[i].vertices[j] = cell_vertices[i][j];
          cells[i].material_id = 0;
        }

      // Here we define the boundaries and the corresponding boundary id
      // explicitly, based on the vertices that lie on them
      SubCellData subcell_data;
      subcell_data.boundary_lines.resize(10);

      subcell_data.boundary_lines[0].vertices = {0,1};
      subcell_data.boundary_lines[0].boundary_id = 0;

      subcell_data.boundary_lines[1].vertices = {1,2};
      subcell_data.boundary_lines[1].boundary_id = 0;

      subcell_data.boundary_lines[2].vertices = {2,3};
      subcell_data.boundary_lines[2].boundary_id = 0;

      subcell_data.boundary_lines[3].vertices = {3,4};
      subcell_data.boundary_lines[3].boundary_id = 0;

      subcell_data.boundary_lines[4].vertices = {4,9};
      subcell_data.boundary_lines[4].boundary_id = 3;

      subcell_data.boundary_lines[5].vertices = {9,11};
      subcell_data.boundary_lines[5].boundary_id = 3;

      subcell_data.boundary_lines[6].vertices = {11,13};
      subcell_data.boundary_lines[6].boundary_id = 3;

      subcell_data.boundary_lines[7].vertices = {13,12};
      subcell_data.boundary_lines[7].boundary_id = 1;

      subcell_data.boundary_lines[8].vertices = {12,5};
      subcell_data.boundary_lines[8].boundary_id = 2;

      subcell_data.boundary_lines[9].vertices = {5,0};
      subcell_data.boundary_lines[9].boundary_id = 2;

      GridTools::consistently_order_cells(cells);

      // Then create the actual mesh:
      coarse_grid.create_triangulation(vertices, cells, subcell_data);
    }


    // Not allow to be 3D
    template <>
    void
    vanKeken45Deep<3>::
    create_coarse_mesh (parallel::distributed::Triangulation<3> &) const
    {
      Assert (false, ExcNotImplemented());
    }


    template <int dim>
    std::set<types::boundary_id>
    vanKeken45Deep<dim>::
    get_used_boundary_indicators () const
    {
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }



    template <>
    std::map<std::string,types::boundary_id>
    vanKeken45Deep<2>::
    get_symbolic_boundary_names_map () const
    {
      static const std::pair<std::string,types::boundary_id> mapping[]
        = { std::pair<std::string,types::boundary_id> ("bottom", 0),
            std::pair<std::string,types::boundary_id> ("top", 1),
            std::pair<std::string,types::boundary_id> ("left",  2),
            std::pair<std::string,types::boundary_id> ("right", 3)
          };
      return std::map<std::string,types::boundary_id> (std::begin(mapping),
                                                       std::end(mapping));
    }

    template <>
    std::map<std::string,types::boundary_id>
    vanKeken45Deep<3>::
    get_symbolic_boundary_names_map () const
    {
      Assert (false, ExcNotImplemented());
      std::map<std::string,types::boundary_id> empty;
      return empty;
    }

    template <int dim>
    Point<dim>
    vanKeken45Deep<dim>::get_extents () const
    {
      Point<dim> extents = {660e3, 600e3};
      return extents;
    }

    template <int dim>
    Point<dim>
    vanKeken45Deep<dim>::get_origin () const
    {
      Point<dim> box_origin = {0, 0};
      return box_origin;
    }

    template <int dim>
    double
    vanKeken45Deep<dim>::
    length_scale () const
    {
      return 1e4;
    }

    template <int dim>
    double
    vanKeken45Deep<dim>::depth(const Point<dim> &position) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = position[d];

      std::vector<double> extents = {660e3, 600e3};
      const double d = extents[dim-1] - (position(dim-1));
      return std::min (std::max (d, 0.), maximal_depth());
    }



    template <int dim>
    double
    vanKeken45Deep<dim>::height_above_reference_surface(const Point<dim> &) const
    {
      return 0.0;
    }



    template <int dim>
    Point<dim>
    vanKeken45Deep<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      p(dim-1) = depth;
      return p;
    }



    template <int dim>
    double
    vanKeken45Deep<dim>::maximal_depth() const
    {
      return 600e3;
    }



    template <int dim>
    bool
    vanKeken45Deep<dim>::point_is_in_domain(const Point<dim> &) const
    {
      return true;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    vanKeken45Deep<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }



    template <int dim>
    void
    vanKeken45Deep<dim>::declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    vanKeken45Deep<dim>::parse_parameters (ParameterHandler &)
    {}
  }
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(vanKeken45Deep,
                                   "vanKeken45Deep box",
                                   "A specialized geometry model used for the van Keken 2008 "
                                   "corner flow style subduction benchmark. While the geometry "
                                   "described by this plugin is just a square, a mesh consisting of non-square cells "
                                   "is required to allow the boundaries of the mesh to lie along "
                                   "the slab-mantle wedge interface, which allows the pressure "
                                   "solution to be continuous."
                                  )
  }
}
