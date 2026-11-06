// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Opm::VeFluidState
 */
#ifndef OPM_VE_FLUID_STATE_HPP
#define OPM_VE_FLUID_STATE_HPP

#include <opm/material/fluidstates/FluidStateCompositionModules.hpp>
#include <opm/material/fluidstates/FluidStateDensityModules.hpp>
#include <opm/material/fluidstates/FluidStateEnthalpyModules.hpp>
#include <opm/material/fluidstates/FluidStateFugacityModules.hpp>
#include <opm/material/fluidstates/FluidStatePressureModules.hpp>
#include <opm/material/fluidstates/FluidStateSaturationModules.hpp>
#include <opm/material/fluidstates/FluidStateTemperatureModules.hpp>
#include <opm/material/fluidstates/FluidStateViscosityModules.hpp>
#include <opm/material/fluidstates/ModularFluidState.hpp>

namespace Opm
{

/*!
 * \brief Fluid state of the immiscible phases of the vertical-equilibrium model.
 *
 * The thermal variant keeps each phase at the temperature of its zone of the
 * layer and stores the enthalpies, the isothermal variant has one temperature.
 * Otherwise as ImmiscibleFluidState.
 *
 * \tparam ValueType The type of the values.
 * \tparam FluidSystem The fluid system.
 * \tparam thermal Whether each phase has its own temperature and an enthalpy.
 */
template <class ValueType, class FluidSystem, bool thermal = true>
class VeFluidState;

// specialization for the thermal case
template <class ValueType, class FluidSystem>
class VeFluidState<ValueType, FluidSystem, true>
    : public ModularFluidState<
          ValueType,
          FluidSystem::numPhases,
          FluidSystem::numComponents,
          FluidStateExplicitPressureModule<ValueType,
                                           FluidSystem::numPhases,
                                           VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitTemperatureModule<ValueType,
                                              FluidSystem::numPhases,
                                              VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateImmiscibleCompositionModule<ValueType,
                                                FluidSystem,
                                                VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitFugacityModule<ValueType,
                                           FluidSystem::numPhases,
                                           FluidSystem::numComponents,
                                           VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitSaturationModule<ValueType,
                                             FluidSystem::numPhases,
                                             VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitDensityModule<ValueType,
                                          FluidSystem::numPhases,
                                          VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitViscosityModule<ValueType,
                                            FluidSystem::numPhases,
                                            VeFluidState<ValueType, FluidSystem, true>>,
          FluidStateExplicitEnthalpyModule<ValueType,
                                           FluidSystem::numPhases,
                                           VeFluidState<ValueType, FluidSystem, true>>>
{
};

// specialization for the isothermal case
template <class ValueType, class FluidSystem>
class VeFluidState<ValueType, FluidSystem, false>
    : public ModularFluidState<
          ValueType,
          FluidSystem::numPhases,
          FluidSystem::numComponents,
          FluidStateExplicitPressureModule<ValueType,
                                           FluidSystem::numPhases,
                                           VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateEquilibriumTemperatureModule<ValueType,
                                                 FluidSystem::numPhases,
                                                 VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateImmiscibleCompositionModule<ValueType,
                                                FluidSystem,
                                                VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateExplicitFugacityModule<ValueType,
                                           FluidSystem::numPhases,
                                           FluidSystem::numComponents,
                                           VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateExplicitSaturationModule<ValueType,
                                             FluidSystem::numPhases,
                                             VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateExplicitDensityModule<ValueType,
                                          FluidSystem::numPhases,
                                          VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateExplicitViscosityModule<ValueType,
                                            FluidSystem::numPhases,
                                            VeFluidState<ValueType, FluidSystem, false>>,
          FluidStateNullEnthalpyModule<ValueType,
                                       FluidSystem::numPhases,
                                       VeFluidState<ValueType, FluidSystem, false>>>
{
};

} // namespace Opm

#endif // OPM_VE_FLUID_STATE_HPP
