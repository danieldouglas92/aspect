# Copyright (C) 2011 - 2025 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file doc/COPYING.  If not see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(Aspect 3.1.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
  "Could not find a valid ASPECT build/installation directory. "
  "Please specify the directory where you are building ASPECT by passing\n"
  "   -D Aspect_DIR=<path to ASPECT>\n"
  "to cmake or by setting the environment variable ASPECT_DIR in your shell "
  "before calling cmake. See the section 'How to write a plugin' in the "
        "manual for more information.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

SET(TARGET "phase-transition-kinetics")
ADD_LIBRARY(${TARGET} SHARED phase-transition-kinetics.cc)
ASPECT_SETUP_PLUGIN(${TARGET})
