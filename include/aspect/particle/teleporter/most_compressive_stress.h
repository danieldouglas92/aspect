/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_teleporter_most_compressive_stress_h
#define _aspect_particle_teleporter_most_compressive_stress_h
#include <aspect/particle/teleporter/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <aspect/global.h>
#include <aspect/parameters.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace Particle
  {
    namespace Teleporter
    {

      template <int dim>
      class MostCompressiveStress : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          void initialize () override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

          bool new_teleportation_step() override;

          void
          local_teleport_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                              const typename ParticleHandler<dim>::particle_iterator &end_particle,
                              const std::vector<Tensor<1,dim>> &vector_field) override;

          std::array<std::pair<double, Tensor<1,dim>>, dim>
          principal_stresses_and_directions(const MaterialModel::MaterialModelInputs<dim> &in,
                                            MaterialModel::MaterialModelOutputs<dim> &out,
                                            const SimulatorAccess<dim> &simulator_access,
                                            const unsigned int q);

          /**
           * Generate particles, and update materialmodel outputs to reflect the removal of melt.
           */
          void create_particles(const SimulatorAccess<dim> &simulator_access,
                                const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out);

          /**
           * Advect the particles along the most compressive stress direction.
           */
          void advect_particles(const SimulatorAccess<dim> &simulator_access);

        private:
          unsigned int allowed_number_of_jumps;
          double melt_volume_threshold;
          double dike_extraction_depth;

          ParticleHandler<dim> teleport_particle_handler;
          bool teleport_particle_handler_initialized = false;
          types::particle_index next_particle_id = 0;

      };
    }
  }
}

#endif
