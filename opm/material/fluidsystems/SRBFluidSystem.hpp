// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
* \file
*
* \copydoc Opm::SRBFluidSystem
*/
#ifndef OPM_SRB_SYSTEM_HPP
#define OPM_SRB_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/components/Brine.hpp>
#include <opm/material/components/SimpleHuDuanH2O.hpp>
#include <opm/material/components/H2.hpp>
#include <opm/material/components/Sulfate.hpp>
#include <opm/material/components/SRBacteria.hpp>

#include <opm/material/common/Valgrind.hpp>

#include <opm/material/binarycoefficients/Brine_H2.hpp>

namespace Opm {

/*!
* \brief A two-phase fluid system for sulfate-reducing bacteria.
*
* This fluid system contains hydrogen, brine, sulfate, and bacteria. 
*/
template <class Scalar>
class SRBFluidSystem : public BaseFluidSystem<Scalar, SRBFluidSystem<Scalar> >
{
    // Extrapolate tabulated properties
    static const bool extrapolate = true;

public:
    // Setup components
    using H2O = Opm::SimpleHuDuanH2O<Scalar>;
    using Brine = Opm::Brine<Scalar, H2O>;
    using H2 = Opm::H2<Scalar>;
    using SO4 = Opm::Sulfate<Scalar>;
    using Bact = Opm::SRBacteria<Scalar>;

    // The binary coefficients for brine and H2
    using BinaryCoeffBrineH2 = BinaryCoeff::Brine_H2<Scalar, H2O, H2>;

    // Phase information
    static constexpr int numPhases = 2;
    static constexpr int liquidPhaseIdx = 0;
    static constexpr int gasPhaseIdx = 1;

    // Component information
    static constexpr int numComponents = 4;
    static constexpr int BrineIdx = 0;
    static constexpr int H2Idx = 1;
    static constexpr int SO4Idx = 2;
    static constexpr int BactIdx = 3;

    // Do not use cache here
    template <class Evaluation>
    struct ParameterCache : public NullParameterCache<Evaluation>
    {};

    /* --------------------------------
    * Short, phase-related functions 
    * --------------------------------- */ 
    /*!
    * \copydoc BaseFluidSystem::phaseName
    */
    static const char* phaseName(unsigned phaseIdx)
    {
            static const char* name[] = {"l",  // liquid phase
                                        "g"};  // gas phase

            assert(0 <= phaseIdx && phaseIdx < 2);
            return name[phaseIdx];
    }

    /*!
    * \copydoc BaseFluidSystem::isCompressible
    */
    static bool isCompressible([[maybe_unused]] unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
    * \copydoc BaseFluidSystem::isIdealGas
    */
    static bool isIdealGas(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        
        return false;
    }

    /*!
    * \copydoc BaseFluidSystem::isIdealMixture
    */
    static bool isIdealMixture([[maybe_unused]] unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
    * \copydoc BaseFluidSystem::isLiquid
    */
    static bool isLiquid(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx == liquidPhaseIdx;
    }

    /* -----------------------------------
    * Short, component-related functions 
    * ------------------------------------ */
    /*!
    * \copydoc BaseFluidSystem::componentName
    */
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            Brine::name(),
            H2::name(),
            SO4::name(),
            Bact::name(),
        };

        assert(compIdx < numComponents);
        return name[compIdx];
    }

    /*!
    * \copydoc BaseFluidSystem::molarMass
    */
    static Scalar molarMass(unsigned compIdx)
    {
        assert(compIdx < numComponents);
        switch (compIdx) {
            case BrineIdx: return Brine::molarMass();
            case H2Idx: return H2::molarMass();
            case SO4Idx: return SO4::molarMass();
            case BactIdx: return Bact::molarMass();
        }
    }

    /* -------------------------
    * Thermodynamic relations 
    * -------------------------- */
    /*!
    * \copydoc BaseFluidSystem::init
    */
    static void init()
    {
        // Set brine salinity
        Brine::salinity = 0.1;
    }

    /*!
    * \copydoc BaseFluidSystem::density
    */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        // Check phase index
        assert(0 <= 0 && phaseIdx < numPhases);

        // Temperature and pressure
        const LhsEval& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));

        // Calculate density
        LhsEval dens;
        if (phaseIdx == liquidPhaseIdx) {
            // Density of brine with dissolved H2
            LhsEval xlH2 = min(1.0, max(0.0,  decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, H2Idx))));
            dens = liquidDensity_(temperature, pressure, xlH2);
        }
        else {
            // Check for gas phase index
            assert(phaseIdx == gasPhaseIdx);

            // Only H2 in gas phase
            dens = H2::gasDensity(temperature, pressure, extrapolate);
        }

        // Return
        Valgrind::CheckDefined(dens);
        return dens;
    }

    /*!
    * \copydoc BaseFluidSystem::viscosity
    */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        // Check phase index
        assert(0 <= 0 && phaseIdx < numPhases);

        // Temperature and pressure
        const LhsEval& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        
        // Calculate viscosity
        LhsEval visc;
        if (phaseIdx == liquidPhaseIdx) {
            // Assume only liquid viscosity = brine viscosity
            visc = Brine::liquidViscosity(temperature, pressure);
        }
        else {
            // Check for gas phase index
            assert(phaseIdx == gasPhaseIdx);

            // Assume gas viscosity = H2 viscosity
            visc = H2::gasViscosity(temperature, pressure);
        }

        // Return
        Valgrind::CheckDefined(visc);
        return visc;
    }
    
    /*!
    * \copydoc BaseFluidSystem::fugacityCoefficient
    */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        // Check indices
        assert(0 <= 0 && phaseIdx < numPhases);
        assert(0 <= 0 && compIdx < numComponents);

        // Temperature and pressure
        const LhsEval& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));
        
        // Calculate fugacity coefficient
        LhsEval phi;
        if (phaseIdx == gasPhaseIdx)
            // Assume ideal mixture for gas phase
            phi = 1.0;
        else {
            // Check phase index
            assert(phaseIdx == liquidPhaseIdx);

            // Only H2 in gas phase
            if (compIdx == H2Idx) {
                // Assuming only H2 in gas phase (yH2 = 1.0), phi = phi_g / xH2. Calculate xH2 first
                LhsEval xH2;
                Scalar salinity = Brine::salinity;
                BinaryCoeffBrineH2::calculateMoleFractions(temperature,
                                                   pressure,
                                                   salinity,
                                                   xH2,
                                                   extrapolate);

                // Use Helmholtz free energy EOS to calculate phi_g
                LhsEval phi_g = BinaryCoeffBrineH2::fugacityCoefficientH2(temperature, pressure, extrapolate);

                // Calculate phi
                phi = phi_g / xH2;
            }
            else {
                // Other components assumed zero
                phi = 1e-20;
            }
        }

        // Return
        Valgrind::CheckDefined(phi);
        return phi;
    }

    /*!
    * \copydoc BaseFluidSystem::diffusionCoefficient
    */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned compIdx)
    {
        // Check indices
        assert(0 <= 0 && phaseIdx < numPhases);
        assert(0 <= 0 && compIdx < numComponents);

        // Temperature and pressure
        const LhsEval& temperature = decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& pressure = decay<LhsEval>(fluidState.pressure(phaseIdx));

        // Calculate diffusion coefficient
        LhsEval dcoef;
        if (phaseIdx == liquidPhaseIdx) {
            // In liquid phase we only have brine-H2 diffusion at the moment
            if (compIdx == H2Idx) {
                // Diffusion coefficient of H2 in pure water according to Ferrell & Himmelbau, AIChE Journal, 13(4), 1967 (Eq.
                // 23)
                // Some intermediate calculations and definitions
                const Scalar vm = 28.45;  // molar volume at normal boiling point (20.271 K and 1 atm) in cm2/mol
                const Scalar sigma = 2.96 * 1e-8; // Lennard-Jones 6-12 potential in cm (1 Å = 1e-8 cm)
                const Scalar avogadro = 6.022e23; // Avogrado's number in mol^-1
                const Scalar alpha = sigma / pow((vm / avogadro), 1 / 3);  // Eq. (19)
                const Scalar lambda = 1.729; // quantum parameter [-]
                const LhsEval& mu_pure = H2O::liquidViscosity(temperature, pressure, extrapolate) * 1e3;  // water viscosity in cP

                // Diffusion coeff in pure water in cm2/s
                const LhsEval D_pure = ((4.8e-7 * temperature) / pow(mu_pure, alpha)) * pow((1 + pow(lambda, 2)) / vm, 0.6);

                // Diffusion coefficient in brine using Ratcliff and Holdcroft, J. G. Trans. Inst. Chem. Eng, 1963. OBS: Value
                // of n is noted as the recommended single value according to Akita, Ind. Eng. Chem. Fundam., 1981.
                const LhsEval& mu_brine = Brine::liquidViscosity(temperature, pressure) * 1e3; // Brine viscosity in cP
                const LhsEval log_D_brine = log10(D_pure) - 0.637 * log10(mu_brine / mu_pure);
                dcoef = pow(10.0, log_D_brine) * 1e-4;  // convert from cm2/s to m2/s
            }
            else
                // Zero diffusion for other components
                dcoef = 0.0;
        }
        else {
            // Check phase index
            assert(phaseIdx == gasPhaseIdx);

            // Only H2 in gas phase
            if (compIdx == H2Idx) {
                dcoef = BinaryCoeffBrineH2::gasDiffCoeff(temperature, pressure);
            }
            else
                // Zero diffusion for other components
                dcoef = 0.0;
        }

        // Return
        Valgrind::CheckDefined(dcoef);
        return dcoef;
    }

private:
    /*!
    * \brief Calculated the density of the aqueous solution where contributions of salinity and dissolved H2 is taken
    * into account.
    * 
    * \param T temperature [K]
    * \param pl liquid pressure [Pa]
    * \param xlH2 mole fraction H2 [-]
    */
    template <class LhsEval>
    static LhsEval liquidDensity_(const LhsEval& T,
                           const LhsEval& pl,
                           const LhsEval& xlH2)
    {
        // check input variables
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(xlH2);

        // check if pressure and temperature is valid
        if(!extrapolate && T < 273.15) {
            const std::string msg = 
                "Liquid density for Brine and H2 is only "
                "defined above 273.15K (is " + 
                std::to_string(getValue(T)) + "K)";
            throw NumericalProblem(msg);
        }
        if(!extrapolate && pl >= 2.5e8) {
            const std::string msg  =
                "Liquid density for Brine and H2 is only "
                "defined below 250MPa (is " +
                std::to_string(getValue(pl)) + "Pa)";
            throw NumericalProblem(msg);
        }

        // calculate individual contribution to density
        const LhsEval& rho_brine = Brine::liquidDensity(T, pl, extrapolate);
        const LhsEval& rho_pure = H2O::liquidDensity(T, pl, extrapolate);
        const LhsEval& rho_lH2 = liquidDensityWaterH2_(T, pl, xlH2);
        const LhsEval& contribH2 = rho_lH2 - rho_pure;

        return rho_brine + contribH2;
    }

    /*!
    * \brief Density of aqueous solution with dissolved H2. Formula from Li et al. (2018) and Garica, Lawrence Berkeley
    * National Laboratory, 2001.
    * 
    * \param temperature [K]
    * \param pl liquid pressure [Pa]
    * \param xlH2 mole fraction [-]
    */
    template <class LhsEval>
    static LhsEval liquidDensityWaterH2_(const LhsEval& temperature,
                                          const LhsEval& pl,
                                          const LhsEval& xlH2)
    {
        // molar masses
        Scalar M_H2 = H2::molarMass();
        Scalar M_H2O = H2O::molarMass();

        // density of pure water
        const LhsEval& rho_pure = H2O::liquidDensity(temperature, pl, extrapolate);

        // (apparent) molar volume of H2, Eq. (14) in Li et al. (2018)
        const LhsEval& A1 = 51.1904 - 0.208062*temperature + 3.4427e-4*(temperature*temperature);
        const LhsEval& A2 = -0.022;
        const LhsEval& V_phi = (A1 + A2 * (pl / 1e6)) / 1e6;  // pressure in [MPa] and Vphi in [m3/mol] (from [cm3/mol])

        // density of solution, Eq. (19) in Garcia (2001)
        const LhsEval xlH2O = 1.0 - xlH2;
        const LhsEval& M_T = M_H2O * xlH2O + M_H2 * xlH2;
        const LhsEval& rho_aq = 1 / (xlH2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));

        return rho_aq;

    }
};  // class SRBFluidSystem

}  // namespace Opm

#endif