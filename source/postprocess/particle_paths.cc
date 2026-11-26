/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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

#include <aspect/postprocess/particle_paths.h>
#include <aspect/particle/manager.h>
#include <aspect/particle/teleporter/most_compressive_stress.h>
#include <aspect/utilities.h>

#include <deal.II/grid/reference_cell.h>

#include <fstream>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      void
      ParticlePathOutput<dim>::build_patches(const std::map<types::particle_index, std::vector<Point<dim>>> &paths)
      {
        dataset_names = {"particle_id", "segment"};
        patches.clear();

        for (const auto &path : paths)
          for (unsigned int s = 0; s + 1 < path.second.size(); ++s)
            {
              DataOutBase::Patch<1,dim> patch;
              patch.reference_cell = ReferenceCells::Line;
              patch.vertices[0] = path.second[s];
              patch.vertices[1] = path.second[s+1];
              patch.patch_index = patches.size();
              patch.data.reinit(dataset_names.size(), 2);
              for (unsigned int i = 0; i < 2; ++i)
                {
                  patch.data(0,i) = path.first;
                  patch.data(1,i) = s;
                }
              patches.push_back(std::move(patch));
            }
      }



      template <int dim>
      const std::vector<DataOutBase::Patch<1,dim>> &
      ParticlePathOutput<dim>::get_patches () const
      {
        return patches;
      }



      template <int dim>
      std::vector<std::string>
      ParticlePathOutput<dim>::get_dataset_names () const
      {
        return dataset_names;
      }
    }



    template <int dim>
    std::pair<std::string,std::string>
    ParticlePaths<dim>::execute (TableHandler &)
    {
      Particle::Teleporter::MostCompressiveStress<dim> *teleporter = nullptr;
      for (unsigned int i = 0; i < this->n_particle_managers() && teleporter == nullptr; ++i)
        teleporter = dynamic_cast<Particle::Teleporter::MostCompressiveStress<dim> *>
                     (this->get_particle_manager(i).get_teleporter());

      AssertThrow(teleporter != nullptr,
                  ExcMessage("The 'particle paths' postprocessor requires the 'most compressive stress' "
                             "particle teleporter."));

      const auto &paths = teleporter->get_particle_paths();

      unsigned int n_local_segments = 0;
      for (const auto &path : paths)
        n_local_segments += path.second.size() > 0 ? path.second.size() - 1 : 0;

      // Collective call, so every process takes the same branch below.
      if (Utilities::MPI::sum(n_local_segments, this->get_mpi_communicator()) == 0)
        {
          teleporter->clear_particle_paths();
          return {"", ""};
        }

      internal::ParticlePathOutput<dim> data_out;
      data_out.build_patches(paths);

      const std::string directory = "particle_paths/";
      const std::string prefix = "particle-paths-" + Utilities::int_to_string(output_file_number, 5);
      const double time = this->convert_output_to_years() ? this->get_time() / year_in_seconds
                          : this->get_time();

      const unsigned int my_id = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      DataOutBase::VtkFlags vtk_flags;
      vtk_flags.cycle = this->get_timestep_number();
      vtk_flags.time = time;
      data_out.set_flags(vtk_flags);

      {
        std::ofstream out(this->get_output_directory() + directory + prefix + "." +
                          Utilities::int_to_string(my_id, 4) + ".vtu");
        AssertThrow(out, ExcMessage("Unable to open particle path output file for writing."));
        data_out.write_vtu(out);
      }

      if (my_id == 0)
        {
          std::vector<std::string> filenames;
          for (unsigned int i = 0; i < n_processes; ++i)
            filenames.push_back(prefix + "." + Utilities::int_to_string(i, 4) + ".vtu");

          const std::string pvtu_name = directory + prefix + ".pvtu";
          std::ofstream pvtu_file(this->get_output_directory() + pvtu_name);
          data_out.write_pvtu_record(pvtu_file, filenames);

          times_and_pvtu_file_names.emplace_back(time, pvtu_name);
          std::ofstream pvd_file(this->get_output_directory() + "particle_paths.pvd");
          DataOutBase::write_pvd_record(pvd_file, times_and_pvtu_file_names);
        }

      ++output_file_number;
      teleporter->clear_particle_paths();

      return {"Writing particle paths:", this->get_output_directory() + directory + prefix};
    }



    template <int dim>
    void
    ParticlePaths<dim>::parse_parameters (ParameterHandler &)
    {
      Utilities::create_directory(this->get_output_directory() + "particle_paths/",
                                  this->get_mpi_communicator(),
                                  true);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
#define INSTANTIATE(dim) \
  template class ParticlePathOutput<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }

    ASPECT_REGISTER_POSTPROCESSOR(ParticlePaths,
                                  "particle paths",
                                  "A postprocessor that writes the paths along which the 'most compressive "
                                  "stress' particle teleporter moved its particles in each time step. Each "
                                  "path is written as a set of line segments in VTU format (with a .pvd "
                                  "time series file) for visualization in Paraview.")
  }
}
