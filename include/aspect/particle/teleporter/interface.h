/*
 Copyright (C) 2015 - 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_teleporter_interface_h
#define _aspect_particle_teleporter_interface_h

#include <aspect/particle/interface.h>
#include <aspect/global.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace Particle
  {
    namespace Teleporter
    {
      using namespace dealii::Particles;

      /**
       * An abstract class defining virtual methods for performing teleportation
       * of particle paths through the simulation velocity field.
       *
       * @ingroup ParticleTeleporter
       */
      template <int dim>
      class Interface : public ParticleInterfaceBase
      {
        public:
          /**
           * Perform a teleportation step of moving the particles of one cell
           * by the specified timestep dt. Implementations of this function
           * must update the particle location. Between calls to this function
           * the velocity at the updated particle positions is evaluated and
           * passed as input argument during the next call.
           *
           * @param [in] begin_particle An iterator to the first particle to be moved.
           * @param [in] end_particle An iterator to the last particle to be moved.
           * @param [in] old_velocities The velocities at t_n, i.e. before the
           * particle movement, for all particles between @p begin_particle
           * and @p end_particle at their current position.
           * @param [in] velocities The velocities at the particle positions
           * at t_{n+1}, i.e. after the particle movement. Note that this is
           * the velocity at the old positions, but at the new time. It is the
           * responsibility of this function to compute the new location of
           * the particles.
           * @param [in] dt The length of the teleportation timestep.
           */
          virtual
          void
          local_teleport_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                              const typename ParticleHandler<dim>::particle_iterator &end_particle,
                              const std::vector<Tensor<1,dim>> &vector_field) = 0;

          /**
           * This function is called at the end of every teleportation step.
           * In case of multi-step teleporters it signals the beginning of a
           * new teleportation step. The default implementation always returns
           * false, which is ok for single-step teleportation methods.
           *
           * @return This function returns true if the teleporter requires
           * another teleportation step. The particle manager will continue
           * to start new teleportation steps until this function returns false.
           */
          virtual bool new_teleportation_step();

          /**
           * Read teleportation related data for a particle specified by particle_id
           * from the data array. This function is called after transferring
           * a particle to the local domain during a teleportation step.
           *
           * @param [in] data A pointer into the data array. The pointer
           * marks the position where this function starts reading.
           * @param [in] particle An iterator pointing to the particle to read
           * the data for.
           * @return The updated position of the pointer into the data array.
           * The return value is @p data advanced by get_data_size() bytes.
           */
          virtual
          const void *
          read_data(const typename ParticleHandler<dim>::particle_iterator &particle,
                    const void *data);

          /**
           * Write teleportation related data to a vector for a particle
           * specified by particle_id. This function is called in cases where
           * particles leave the local domain during a teleportation step to
           * transfer this data to another process.
           *
           * @param [in] particle An iterator pointing to the particle to
           * write the data for.
           * @param [in] data A pointer into the array of teleportation data.
           * @return The updated position of the pointer into the data array.
           * The return value is @p data advanced by get_data_size() bytes.
           */
          virtual
          void *
          write_data(const typename ParticleHandler<dim>::particle_iterator &particle,
                     void *data) const;
      };


      /**
       * Return a list of names (separated by '|') of possible teleportation
       * classes for particles.
       */
      std::string
      teleporter_object_names ();


      /**
       * Register a particle teleporter so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle teleporter
       * @param description A text description of what this teleporter does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle teleporter wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle teleporter.
       *
       * @ingroup ParticleTeleporters
       */
      template <int dim>
      void
      register_particle_teleporter (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    std::unique_ptr<Interface<dim>> (*factory_function) ());

      /**
       * A function that given the name of a model returns a pointer to an
       * object that describes it. Ownership of the pointer is transferred to
       * the caller.

       * The model object returned is not yet initialized and has not
       * read its runtime parameters yet.

       * @ingroup ParticleTeleporter
       */
      template <int dim>
      std::unique_ptr<Interface<dim>>
      create_particle_teleporter (ParameterHandler &prm);


      /**
       * Declare the runtime parameters of the registered particle teleporters.
       *
       * @ingroup ParticleTeleporter
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);


      /**
       * For the current plugin subsystem, write a connection graph of all of the
       * plugins we know about, in the format that the
       * programs dot and neato understand. This allows for a visualization of
       * how all of the plugins that ASPECT knows about are interconnected, and
       * connect to other parts of the ASPECT code.
       *
       * @param output_stream The stream to write the output to.
       */
      template <int dim>
      void
      write_plugin_graph (std::ostream &output_stream);

      /**
       * Given a class name, a name, and a description for the parameter file
       * for a particle teleporter, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticleTeleporter
       */
#define ASPECT_REGISTER_PARTICLE_TELEPORTER(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_TELEPORTER_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Teleporter::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::Particle::Teleporter::register_particle_teleporter<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Teleporter::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Particle::Teleporter::register_particle_teleporter<3>, \
                                name, description); \
  }
    }
  }
}


#endif
