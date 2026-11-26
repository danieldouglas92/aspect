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

#ifndef _aspect_postprocess_particle_paths_h
#define _aspect_postprocess_particle_paths_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * Writes particle paths as a list of line segments (one patch per segment)
       * so that they can be read by Paraview.
       */
      template <int dim>
      class ParticlePathOutput : public dealii::DataOutInterface<1,dim>
      {
        public:
          void build_patches(const std::map<types::particle_index, std::vector<Point<dim>>> &paths);

        private:
          const std::vector<DataOutBase::Patch<1,dim>> &
          get_patches () const override;

          std::vector<std::string>
          get_dataset_names () const override;

          std::vector<std::string> dataset_names;
          std::vector<DataOutBase::Patch<1,dim>> patches;
      };
    }

    /**
     * A postprocessor that writes the paths along which the 'most compressive stress'
     * particle teleporter moved its particles during each time step. The output
     * consists of one VTU file per process and time step, plus a .pvtu/.pvd file
     * set that can be loaded into Paraview.
     */
    template <int dim>
    class ParticlePaths : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        void parse_parameters (ParameterHandler &prm) override;

      private:
        std::vector<std::pair<double, std::string>> times_and_pvtu_file_names;
        unsigned int output_file_number = 0;
    };
  }
}

#endif
