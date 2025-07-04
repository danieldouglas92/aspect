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


#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <aspect/utilities.h>
#include <deal.II/dofs/dof_tools.h>

namespace aspect
{
  namespace GeometryModel
  {
    namespace
    {
      void append_face_to_subcell_data(SubCellData &subcell_data, const CellData<1> &face)
      {
        subcell_data.boundary_lines.push_back(face);
      }



      void append_face_to_subcell_data(SubCellData &subcell_data, const CellData<2> &face)
      {
        subcell_data.boundary_quads.push_back(face);
      }
    }


    namespace internal
    {
      template <int dim>
      SphericalManifoldWithTopography<dim>::
      SphericalManifoldWithTopography(const std::shared_ptr<const InitialTopographyModel::Interface<dim>> topography,
                                      const double inner_radius,
                                      const double outer_radius)
        :
        SphericalManifold<dim>(Point<dim>()),
        topo (topography),
        R0 (inner_radius),
        R1 (outer_radius)
      {}


      template <int dim>
      std::unique_ptr<Manifold<dim, dim>>
      SphericalManifoldWithTopography<dim>::clone() const
      {
        return std::make_unique<SphericalManifoldWithTopography<dim>>(*this);
      }



      template <int dim>
      double
      SphericalManifoldWithTopography<dim>::
      topography_for_point(const Point<dim> &x_y_z) const
      {
        if (dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(topo.get()) != nullptr)
          return 0;
        else
          {
            // The natural coordinate system of the sphere geometry is r/phi/theta.
            // This is what we need to query the topography with. So start by
            // converting into this coordinate system
            const std::array<double, dim> r_phi_theta = Utilities::Coordinates::cartesian_to_spherical_coordinates(x_y_z);

            // Grab lon,lat coordinates
            Point<dim-1> surface_point;
            for (unsigned int d=0; d<dim-1; ++d)
              surface_point[d] = r_phi_theta[d+1];
            return topo->value(surface_point);
          }
      }



      template <int dim>
      Point<dim>
      SphericalManifoldWithTopography<dim>::
      push_forward_from_sphere(const Point<dim> &p) const
      {
        const double topography = topography_for_point(p);

        // We can short-circuit the computations here if the topography
        // is zero.
        if (topography == 0)
          return p;
        else
          {
            // For the case with (non-zero) topography, we have to project
            // the point out/inward.

            // We start from the undeformed sphere. So the current radius must be
            // between R0 and R1. Check this, with a slight tolerance to account
            // for the fact that our mesh is not *completely* spherical.
            const double r = p.norm();
            Assert (r>=R0*0.99, ExcInternalError());
            Assert (r<=R1*1.01, ExcInternalError());

            // Then we need to stretch the point p outward by a factor that is
            // zero at the core-mantle boundary and results in the right topography
            // at the top.
            //
            // The stretching factor is relative to the distance from the core
            // mantle boundary.
            const double cmb_stretching_factor = (R1+topography-R0)/(R1-R0);

            // From this we can compute what the desired radius is going to be, and
            // scale the given point accordingly:
            const double new_radius = (r-R0)*cmb_stretching_factor + R0;
            return p * (new_radius/r);
          }
      }



      template <int dim>
      Point<dim>
      SphericalManifoldWithTopography<dim>::
      pull_back_to_sphere(const Point<dim> &p) const
      {
        const double topography = topography_for_point(p);

        // We can short-circuit the computations here if the topography
        // is zero.
        if (topography == 0)
          return p;
        else
          {
            // For the case with (non-zero) topography, we have to project
            // the point in/outward.

            // We start from the deformed sphere. So the current radius must be
            // between R0 and R1+topography. Check this, with a slight tolerance
            // to account for the fact that our mesh is not *completely* spherical.
            const double r = p.norm();
            Assert (r>=R0*0.99, ExcInternalError());
            Assert (r<=(R1+topography)*1.01, ExcInternalError());

            // Then we need to stretch(=shrink) the point p outward by a factor that is
            // zero at the core-mantle boundary and results in the right topography
            // at the top.
            const double cmb_stretching_factor = (R1-R0)/(R1+topography-R0);

            // From this we can compute what the desired radius is going to be, and
            // scale the given point accordingly:
            const double new_radius = (r-R0)*cmb_stretching_factor + R0;
            return p * (new_radius/r);
          }
      }



      template <int dim>
      Point<dim>
      SphericalManifoldWithTopography<dim>::
      get_intermediate_point(const Point<dim> &p1,
                             const Point<dim> &p2,
                             const double      w) const
      {
        return
          push_forward_from_sphere
          (SphericalManifold<dim>::get_intermediate_point (pull_back_to_sphere(p1),
                                                           pull_back_to_sphere(p2), w));
      }



      template <int dim>
      Tensor<1, dim>
      SphericalManifoldWithTopography<dim>::
      get_tangent_vector(const Point<dim> &x1,
                         const Point<dim> &x2) const
      {
        // TODO: Deal with pull back and push forward
        return SphericalManifold<dim>::get_tangent_vector (x1, x2);
      }



      template <int dim>
      Tensor<1, dim>
      SphericalManifoldWithTopography<dim>::
      normal_vector(const typename Triangulation<dim, dim>::face_iterator &face,
                    const Point<dim> &p) const
      {
        // We calculate radial, rather than *normal*, vectors here if a face is
        // at the boundary. This is as described in the documentation of this
        // function.

        // TODO: Add an input parameter that determines whether we use this
        //   or the "geometrically correct" behavior.
        return SphericalManifold<dim>::normal_vector (face, p);
      }



      template <int dim>
      void
      SphericalManifoldWithTopography<dim>::
      get_normals_at_vertices(
        const typename Triangulation<dim, dim>::face_iterator &face,
        typename Manifold<dim, dim>::FaceVertexNormals &face_vertex_normals) const
      {
        // TODO: Deal with pull back and push forward
        SphericalManifold<dim>::get_normals_at_vertices(face, face_vertex_normals);
      }



      template <int dim>
      void
      SphericalManifoldWithTopography<dim>::
      get_new_points(const ArrayView<const Point<dim>> &surrounding_points,
                     const Table<2, double>            &weights,
                     ArrayView<Point<dim>>              new_points) const
      {
        // First compute the points in the untransformed coordinate system
        // without surface topography:
        std::vector<Point<dim>> untransformed_points;
        untransformed_points.reserve(surrounding_points.size());
        for (const Point<dim> &p : surrounding_points)
          untransformed_points.emplace_back (pull_back_to_sphere(p));

        // Call the function in the base class:
        SphericalManifold<dim>::get_new_points(untransformed_points, weights, new_points);

        // Since the function in the base class returned its information
        // in a writable array, we can in-place push forward to the deformed
        // sphere with topography:
        for (Point<dim> &p : new_points)
          p = push_forward_from_sphere(p);
      }



      template <int dim>
      Point<dim>
      SphericalManifoldWithTopography<dim>::
      get_new_point(const ArrayView<const Point<dim>> &surrounding_points,
                    const ArrayView<const double>          &weights) const
      {
        // First compute the points in the untransformed coordinate system
        // without surface topography:
        std::vector<Point<dim>> untransformed_points;
        untransformed_points.reserve(surrounding_points.size());
        for (const Point<dim> &p : surrounding_points)
          untransformed_points.emplace_back (pull_back_to_sphere(p));

        // Call the function in the base class, and transform the
        // returned value:
        return push_forward_from_sphere(SphericalManifold<dim>::get_new_point(untransformed_points,
                                                                              weights));
      }
    }



    template <int dim>
    void
    SphericalShell<dim>::initialize ()
    {
      manifold = std::make_unique<internal::SphericalManifoldWithTopography<dim>>(this->get_initial_topography_model_pointer(),
                                                                                   R0, R1);
    }



    template <int dim>
    void
    SphericalShell<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      AssertThrow (phi == 360 || phi == 90 || ((phi == 180) && (dim == 2)),
                   ExcMessage ("The only opening angles that are allowed for "
                               "this geometry are 90, 180, and 360 in 2d; "
                               "and 90 and 360 in 3d."));

      // Custom mesh extrusion is only implemented for full spherical shells
      // for now, so check that custom mesh schemes are only being used when
      // the opening angle is 360 degrees.
      AssertThrow (phi == 360 || (custom_mesh == none),
                   ExcMessage ("The only opening angle that is allowed for "
                               "this geometry with a custom mesh is 360."));

      if (phi == 360)
        {
          if (custom_mesh == none)
            {
              // If we are not using a custom mesh scheme, the mesh is generated
              // as per the original code.
              GridGenerator::hyper_shell (coarse_grid,
                                          Point<dim>(),
                                          R0,
                                          R1,
                                          (n_cells_along_circumference == 0
                                           ?
                                           // automatic choice that leads to reasonable
                                           // meshes with the typical aspect ratio of
                                           // the Earth
                                           (dim==3 ? 96 : 12)
                                           :
                                           // user choice
                                           n_cells_along_circumference),
                                          true);
            }
          else
            {
              Assert(custom_mesh == slices || custom_mesh == list, ExcNotImplemented());
              // If we are using a custom mesh scheme, we need to create
              // a new triangulation to extrude (this will be a 1D line in
              // 2d space, or a 2d surface in 3d space).
              Triangulation<dim-1,dim> sphere_mesh;
              GridGenerator::hyper_sphere (sphere_mesh);
              sphere_mesh.refine_global (initial_lateral_refinement);

              // Calculate the number of R_values wrt custom mesh scheme
              unsigned int n_R_values;
              if (custom_mesh == slices)
                n_R_values = n_slices+1;
              else
                n_R_values = R_values_list.size()+2;

              // Allocate R_values wrt the number of slices
              std::vector<double> R_values (n_R_values);
              if (custom_mesh == slices)
                {
                  for (unsigned int s=0; s<n_R_values; ++s)
                    R_values[s] = R0 + (R1-R0)/n_slices * s;
                }
              if (custom_mesh == list)
                {
                  R_values[0] = R0;
                  R_values[n_R_values-1] = R1;
                  for (unsigned int s=1; s<(n_R_values-1); ++s)
                    R_values[s] = R_values_list[s-1];
                }
              std::vector<Point<dim>>    points(R_values.size() * sphere_mesh.n_vertices());


              // Copy the array of points as many times as there will be slices,
              // one slice at a time. The z-axis value are defined in slices_coordinates
              for (unsigned int point_layer = 0; point_layer < R_values.size(); ++point_layer)
                {
                  for (unsigned int vertex_n = 0; vertex_n < sphere_mesh.n_vertices();
                       ++vertex_n)
                    {
                      const Point<dim> vertex = sphere_mesh.get_vertices()[vertex_n];
                      points[point_layer * sphere_mesh.n_vertices() + vertex_n] =
                        vertex * R_values[point_layer];
                    }
                }

              // Then create the cells of each of the slices, one stack at a
              // time
              std::vector<CellData<dim>> cells;
              cells.reserve((R_values.size() - 1) * sphere_mesh.n_active_cells());

              SubCellData               subcell_data;

              for (const auto &cell : sphere_mesh.active_cell_iterators())
                {
                  for (unsigned int cell_layer = 0; cell_layer < R_values.size() - 1; ++cell_layer)
                    {
                      CellData<dim> this_cell;
                      for (const unsigned int vertex_n : cell->vertex_indices())
                        {
                          this_cell.vertices[vertex_n] =
                            cell->vertex_index(vertex_n) + cell_layer * sphere_mesh.n_vertices();
                          this_cell.vertices[vertex_n + cell->n_vertices()] =
                            cell->vertex_index(vertex_n) +
                            (cell_layer + 1) * sphere_mesh.n_vertices();
                        }

                      this_cell.material_id = cell->material_id();


                      Assert(GridTools::cell_measure (points, this_cell.vertices) > 0, ExcInternalError());

                      cells.push_back(this_cell);

                      // Mark the bottom face of the cell as boundary 0 if we are in
                      // the bottom layer of cells
                      if (cell_layer == 0)
                        {
                          CellData<dim-1> face;
                          for (const unsigned int vertex_n : cell->vertex_indices())
                            face.vertices[vertex_n] =
                              cell->vertex_index(vertex_n) + cell_layer * sphere_mesh.n_vertices();
                          face.boundary_id = 0;

                          append_face_to_subcell_data(subcell_data, face);
                        }

                      // Mark the top face of the cell as boundary 1 if we are in
                      // the top layer of cells
                      if (cell_layer == R_values.size()-2)
                        {
                          CellData<dim-1> face;
                          for (const unsigned int vertex_n : cell->vertex_indices())
                            face.vertices[vertex_n] =
                              cell->vertex_index(vertex_n) +
                              (cell_layer + 1) * sphere_mesh.n_vertices();
                          face.boundary_id = 1;

                          append_face_to_subcell_data(subcell_data, face);
                        }

                    }
                }

              // Then create the actual mesh:
              GridTools::consistently_order_cells (cells);
              coarse_grid.create_triangulation(points, cells, subcell_data);
            }
        }
      else if (phi == 90)
        {
          GridGenerator::quarter_hyper_shell (coarse_grid,
                                              Point<dim>(),
                                              R0,
                                              R1,
                                              0,
                                              true);

          // there was a bug with boundary colorization of thin shells
          // before deal.II 9.7. Use a fixed version of that function,
          // for deal.II versions that need it.
#if !DEAL_II_VERSION_GTE(9,7,0)
          if (dim == 3)
            colorize_quarter_hyper_shell(coarse_grid,
                                         Point<dim>(),
                                         R0,
                                         R1);
#endif

          if (periodic)
            {
              // Tell p4est about the periodicity of the mesh.
              std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator>>
              matched_pairs;
              FullMatrix<double> rotation_matrix(dim);
              rotation_matrix[0][1] = 1.;
              rotation_matrix[1][0] = -1.;

              GridTools::collect_periodic_faces(coarse_grid, /*b_id1*/ 2, /*b_id2*/ 3,
                                                /*direction*/ 1, matched_pairs,
                                                Tensor<1, dim>(), rotation_matrix);

              if (matched_pairs.size() > 0)
                coarse_grid.add_periodicity (matched_pairs);
            }
        }
      else if (phi == 180)
        {
          GridGenerator::half_hyper_shell (coarse_grid,
                                           Point<dim>(),
                                           R0,
                                           R1,
                                           0,
                                           true);
        }
      else
        {
          Assert (false, ExcInternalError());
        }

      // Then add the topography to the mesh by moving all vertices from the
      // undeformed sphere to the sphere with topography:
      GridTools::transform (
        [&](const Point<dim> &p) -> Point<dim>
      {
        return manifold->push_forward_from_sphere(p);
      },
      coarse_grid);

      // Use a manifold description for all cells. Use manifold_id 99 in order
      // not to step on the boundary indicators used below.
      coarse_grid.set_manifold (my_manifold_id, *manifold);
      set_manifold_ids(coarse_grid);
    }



    template <int dim>
    void
    SphericalShell<dim>::set_manifold_ids (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        cell->set_all_manifold_ids (my_manifold_id);
    }



    template <int dim>
    std::set<types::boundary_id>
    SphericalShell<dim>::
    get_used_boundary_indicators () const
    {
      // Follow what is described in the documentation of this class.
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      if (phi == 360)
        {
          const types::boundary_id s[] = { 0, 1 };
          return std::set<types::boundary_id>(std::begin(s), std::end(s));
        }
      else if (phi == 90 && dim == 3)
        {
          const types::boundary_id s[] = { 0, 1, 2, 3, 4};
          return std::set<types::boundary_id>(std::begin(s), std::end(s));
        }
      else
        {
          const types::boundary_id s[] = { 0, 1, 2, 3 };
          return std::set<types::boundary_id>(std::begin(s), std::end(s));
        }
    }



    template <int dim>
    std::map<std::string,types::boundary_id>
    SphericalShell<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
            = { {"bottom", 0},
              {"top", 1},
              {"left",  2},
              {"right", 3}
            };

            if (phi == 360)
              return std::map<std::string,types::boundary_id> (std::begin(mapping),
                                                               std::begin(mapping)+2);
            else
              return std::map<std::string,types::boundary_id> (std::begin(mapping),
                                                               std::begin(mapping)+4);
          }

          case 3:
          {
            if (phi == 360)
              {
                return
                {
                  {"bottom", 0},
                  {"top",    1}
                };
              }
            else if (phi == 90)
              {
                return
                {
                  {"bottom", 0},
                  {"top",    1},
                  {"east",   2},
                  {"west",   3},
                  {"south",  4}
                };
              }
            else
              Assert (false, ExcNotImplemented());
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }



    template <int dim>
    std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>
    SphericalShell<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>> periodic_boundaries;
      if (periodic)
        {
          periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2, 3), 1) );
        }
      return periodic_boundaries;
    }



    template <int dim>
    void
    SphericalShell<dim>::adjust_positions_for_periodicity (Point<dim> &position,
                                                           const ArrayView<Point<dim>> &connected_positions,
                                                           const ArrayView<Tensor<1, dim>> &connected_velocities) const
    {
      AssertThrow(dim == 2,
                  ExcMessage("Periodic boundaries currently "
                             "only work with 2d spherical shell."));
      AssertThrow(phi == 90,
                  ExcMessage("Periodic boundaries currently "
                             "only work with 90 degree opening angle in spherical shell."));

      if (periodic)
        {
          // define a rotation matrix for the new position depending on the boundary
          Tensor<2,dim> rotation_matrix;

          if (position[0] < 0.)
            {
              rotation_matrix[0][1] = 1.;
              rotation_matrix[1][0] = -1.;
            }
          else if (position[1] < 0.)
            {
              rotation_matrix[0][1] = -1.;
              rotation_matrix[1][0] = 1.;
            }
          else
            return;

          position = rotation_matrix * position;

          for (auto &connected_position: connected_positions)
            connected_position = rotation_matrix * connected_position;

          for (auto &connected_velocity: connected_velocities)
            connected_velocity = rotation_matrix * connected_velocity;
        }

      return;
    }



    template <int dim>
    double
    SphericalShell<dim>::
    length_scale () const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. So use a formulation that yields this
      // value for the R0,R1 corresponding to earth, and that more
      // generally scales with (R1-R0), which we can get by calling
      // maximal_depth(). In essence, the factor in front of the call
      // to maximal_depth() is just a magic number that has turned out
      // to work well.
      return (1e4 / (6336000.-3481000.)) * maximal_depth();
    }



    template <int dim>
    double
    SphericalShell<dim>::depth(const Point<dim> &position) const
    {
      if (this->simulator_is_past_initialization() &&
          !Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()))
        return std::min(std::max (R1 + manifold->topography_for_point(position) - position.norm(), 0.), maximal_depth());
      else
        return std::min (std::max (R1-position.norm(), 0.), maximal_depth());
    }



    template <int dim>
    double
    SphericalShell<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return position.norm()-R1;
    }



    template <int dim>
    Point<dim>
    SphericalShell<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point along the axes toward the north pole, at the
      // requested depth.
      Point<dim> p;

      p[dim-1] = std::min (std::max(R1 + manifold->topography_for_point(p) - depth, R0), R1);

      return p;
    }



    template <int dim>
    double
    SphericalShell<dim>::maximal_depth() const
    {
      return R1 + this->get_initial_topography_model().max_topography() - R0;
    }



    template <int dim>
    double SphericalShell<dim>::inner_radius () const
    {
      return R0;
    }



    template <int dim>
    double SphericalShell<dim>::outer_radius () const
    {
      return R1;
    }



    template <int dim>
    double SphericalShell<dim>::opening_angle () const
    {
      return phi;
    }



    template <int dim>
    bool
    SphericalShell<dim>::has_curved_elements () const
    {
      return true;
    }



    template <int dim>
    bool
    SphericalShell<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  this->simulator_is_past_initialization() == false,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      const std::array<double, dim> spherical_point = Utilities::Coordinates::cartesian_to_spherical_coordinates(point);

      std::array<double, dim> point1, point2;
      point1[0] = R0;
      if (this->simulator_is_past_initialization() &&
          !Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()))
        point2[0] =  R1 + manifold->topography_for_point(point);
      else
        point2[0] = R1;
      point1[1] = 0.0;
      point2[1] = phi * constants::degree_to_radians;
      if (dim == 3)
        {
          point1[2] = 0.0;

          if (phi == 90.0)
            // Octant
            point2[2] = 0.5 * numbers::PI;
          else
            // Full shell
            point2[2] = numbers::PI;
        }

      for (unsigned int d = 0; d < dim; ++d)
        if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
            spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
          return false;

      return true;
    }



    template <int dim>
    std::array<double,dim>
    SphericalShell<dim>::cartesian_to_natural_coordinates(const Point<dim> &position) const
    {
      return Utilities::Coordinates::cartesian_to_spherical_coordinates<dim>(position);
    }



    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    SphericalShell<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::spherical;
    }



    template <int dim>
    Point<dim>
    SphericalShell<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position) const
    {
      return Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(position);
    }



    template <int dim>
    void
    SphericalShell<dim>::make_periodicity_constraints(const DoFHandler<dim> &dof_handler,
                                                      AffineConstraints<double> &constraints) const
    {
      if (periodic)
        {
          std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
          matched_pairs;
          FullMatrix<double> rotation_matrix(dim);
          rotation_matrix[0][1] = 1.;
          rotation_matrix[1][0] = -1.;

          GridTools::collect_periodic_faces(dof_handler, /*b_id1*/ 2, /*b_id2*/ 3,
                                            /*direction*/ 1, matched_pairs,
                                            Tensor<1, dim>(), rotation_matrix);

          DoFTools::make_periodicity_constraints<dim,dim,double>(matched_pairs,
                                                                 constraints,
                                                                 ComponentMask(),
          {0},
          1.);
        }
    }



    template <int dim>
    void
    SphericalShell<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Spherical shell");
        {
          prm.declare_entry ("Custom mesh subdivision", "none",
                             Patterns::Selection ("none|list of radial values|number of slices"),
                             "Choose how the spherical shell mesh is generated. "
                             "By default, a coarse mesh is generated with respect to the "
                             "inner and outer radius, and an initial number of cells along "
                             "circumference. "
                             "In the other cases, a surface mesh is first generated and "
                             "refined as desired, before it is extruded radially following "
                             "the specified subdivision scheme.");
          prm.declare_entry ("List of radial values", "",
                             Patterns::List(Patterns::Double ()),
                             "List of radial values for the custom mesh scheme. Units: "
                             "$\\si{m}$. "
                             "A list of radial values subdivides the spherical shell at "
                             "specified radii. The list must be strictly ascending, and "
                             "the first value must be greater than the inner radius "
                             "while the last must be less than the outer radius.");
          prm.declare_entry ("Number of slices", "1",
                             Patterns::Integer (1),
                             "Number of slices for the custom mesh subdivision scheme. "
                             "The number of slices subdivides the spherical shell into N "
                             "slices of equal thickness. Must be greater than 0.");
          prm.declare_entry ("Initial lateral refinement", "0",
                             Patterns::Integer (0),
                             "Initial lateral refinement for the custom mesh subdivision "
                             "schemes."
                             "The number of refinement steps performed on the initial "
                             "coarse surface mesh, before the surface is extruded "
                             "radially. This parameter allows the user more control "
                             "over the ratio between radial and lateral refinement of "
                             "the mesh.");
          prm.declare_entry ("Inner radius", "3481000.",  // 6371-2890 in km
                             Patterns::Double (0.),
                             "Inner radius of the spherical shell. Units: \\si{\\meter}."
                             "\n\n"
                             ":::{note}\n"
                             "The default value of 3,481,000 m equals the "
                             "radius of a sphere with equal volume as Earth (i.e., "
                             "6371 km) minus the average depth of the core-mantle "
                             "boundary (i.e., 2890 km).\n"
                             ":::");
          prm.declare_entry ("Outer radius", "6336000.",  // 6371-35 in km
                             Patterns::Double (0.),
                             "Outer radius of the spherical shell. Units: \\si{\\meter}."
                             "\n\n"
                             ":::{note}\n"
                             "The default value of 6,336,000 m equals the "
                             "radius of a sphere with equal volume as Earth (i.e., "
                             "6371 km) minus the average depth of the mantle-crust "
                             "interface (i.e., 35 km).\n"
                             ":::");
          prm.declare_entry ("Opening angle", "360.",
                             Patterns::Double (0., 360.),
                             "Opening angle in degrees of the section of the shell "
                             "that we want to build. "
                             "The only opening angles that are allowed for "
                             "this geometry are 90, 180, and 360 in 2d; "
                             "and 90 and 360 in 3d. "
                             "Units: degrees.");
          prm.declare_entry ("Cells along circumference", "0",
                             Patterns::Integer (0),
                             "The number of cells in circumferential direction that are "
                             "created in the coarse mesh in 2d. If zero, this number "
                             "is chosen automatically in a way that produces meshes "
                             "in which cells have a reasonable aspect ratio for models "
                             "in which the depth of the mantle is roughly that of the "
                             "Earth. For planets with much shallower mantles and larger "
                             "cores, you may want to chose a larger number to avoid "
                             "cells that are elongated in tangential and compressed in "
                             "radial direction."
                             "\n\n"
                             "In 3d, the number of cells is computed differently and does "
                             "not have an easy interpretation. Valid values for this parameter "
                             "in 3d are 0 (let this class choose), 6, 12 and 96. "
                             "Other possible values may be discussed in the documentation "
                             "of the deal.II function GridGenerator::hyper_shell. "
                             "The parameter is best left at its default in 3d."
                             "\n\n"
                             "In either case, this parameter is ignored unless the opening "
                             "angle of the domain is 360 degrees. This parameter is also "
                             "ignored when using a custom mesh subdivision scheme.");
          prm.declare_entry ("Phi periodic", "false",
                             Patterns::Bool (),
                             "Whether the shell should be periodic in the phi direction.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SphericalShell<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Spherical shell");
        {
          R0  = prm.get_double ("Inner radius");
          R1  = prm.get_double ("Outer radius");
          phi = prm.get_double ("Opening angle");
          n_cells_along_circumference = prm.get_integer ("Cells along circumference");
          initial_lateral_refinement = prm.get_integer ("Initial lateral refinement");
          R_values_list = Utilities::string_to_double(Utilities::split_string_list(prm.get("List of radial values")));
          n_slices = prm.get_integer ("Number of slices");

          if (prm.get ("Custom mesh subdivision") == "none")
            custom_mesh = none;
          else if (prm.get ("Custom mesh subdivision") == "list of radial values")
            custom_mesh = list;
          else if (prm.get ("Custom mesh subdivision") == "number of slices")
            custom_mesh = slices;
          else
            AssertThrow (false, ExcMessage ("Not a valid custom mesh subdivision scheme."));

          // Check that inner radius is less than outer radius
          AssertThrow (R0 < R1,
                       ExcMessage ("Inner radius must be less than outer radius."));

          // If we are using list of radial values for a custom mesh
          if (custom_mesh == list)
            {
              // Check that list is in ascending order
              for (unsigned int i = 1; i < R_values_list.size(); ++i)
                AssertThrow(R_values_list[i] > R_values_list[i-1],
                            ExcMessage("Radial values must be strictly ascending"));
              // Check that first value is not smaller than the inner radius
              AssertThrow(R_values_list[1] > R0,
                          ExcMessage("First value in List of radial values must be greater than inner radius"));
              // Check that last layer is not larger than the outer radius
              AssertThrow( *(R_values_list.end()-1) < R1,
                           ExcMessage("Last value in List of radial values must be less than outer radius"));
            }


          // If we are extruding the mesh according to a number of slices
          if (custom_mesh == slices)
            {
              AssertThrow (n_slices > 0, ExcMessage("You must set a positive number of slices for extrusion"));
            }

          periodic = prm.get_bool ("Phi periodic");
          if (periodic)
            {
              AssertThrow (dim == 2,  ExcMessage("Periodic boundaries in the spherical shell are only supported for 2d models."));
              AssertThrow (phi == 90, ExcMessage("Periodic boundaries in the spherical shell are only supported for an opening angle of 90 degrees."));
            }


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
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(SphericalShell,
                                   "spherical shell",
                                   "A geometry representing a spherical shell or a piece of it. "
                                   "Inner and outer radii are read from the parameter file "
                                   "in subsection 'Spherical shell'."
                                   "\n\n"
                                   "The spherical shell may be generated as per the original "
                                   "code (with respect to the inner and outer radius, and an "
                                   "initial number of cells along circumference) or following "
                                   "a custom mesh scheme: list of radial values or number of "
                                   "slices. A surface mesh is first generated and refined as "
                                   "desired, before it is extruded radially. A list of radial "
                                   "values subdivides the spherical shell at specified radii. "
                                   "The number of slices subdivides the spherical shell into N "
                                   "slices of equal thickness. The custom spherical shell only "
                                   "works with an opening angle of 360 degrees."
                                   "\n\n"
                                   "Despite the name, this geometry does not imply the use of "
                                   "a spherical coordinate system when used in 2d. Indeed, "
                                   "in 2d the geometry is simply an annulus in a Cartesian "
                                   "coordinate system and consequently would correspond to "
                                   "a cross section of the fluid filled space between two "
                                   "infinite cylinders where one has made the assumption that "
                                   "the velocity in direction of the cylinder axes is zero. "
                                   "This is consistent with the definition of what we consider "
                                   "the two-dimension case given in "
                                   "Section~\\ref{sec:methods:2d-models}."
                                   "\n\n"
                                   "The model assigns boundary indicators as follows: In 2d, "
                                   "inner and outer boundaries get boundary indicators zero "
                                   "and one, and if the opening angle set in the input file "
                                   "is less than 360, then left and right boundaries are "
                                   "assigned indicators two and three. These boundaries can "
                                   "also be referenced using the symbolic names `inner', `outer' "
                                   "and (if applicable) `left', `right'."
                                   "\n\n"
                                   "In 3d, inner and "
                                   "outer indicators are treated as in 2d. If the opening "
                                   "angle is chosen as 90 degrees, i.e., the domain is the "
                                   "intersection of a spherical shell and the first octant, "
                                   "then indicator 2 is at the face $x=0$, 3 at $y=0$, "
                                   "and 4 at $z=0$. These last three boundaries can then also "
                                   "be referred to as `east', `west' and `south' symbolically "
                                   "in input files.")
  }
}
