// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2019 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_depth_composition_volume_flux_statistics_h
#define _aspect_postprocess_depth_composition_volume_flux_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the depth composition volume flux
     * in and out of boundaries.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class DepthCompositionVolumeFluxStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the volumeflux of each composition at a list of specified depths.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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

      private:
        /**
         * The depths at which the compositional volume fluxes are calculated.
         */
        std::vector<double> flux_depths;

        /**
         * The number of depth zones for which the fluxes are evaluated.
         */
        unsigned int n_depth_zones;
    };
  }
}


#endif
