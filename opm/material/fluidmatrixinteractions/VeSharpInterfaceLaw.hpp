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
 * \copydoc Opm::VeSharpInterfaceLaw
 */
#ifndef OPM_VE_SHARP_INTERFACE_LAW_HPP
#define OPM_VE_SHARP_INTERFACE_LAW_HPP

#include "VeSharpInterfaceLawParams.hpp"

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <type_traits>

namespace Opm
{

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Vertically upscaled saturation functions of the sharp-interface
 *        vertical-equilibrium (VE) model.
 *
 * A non-wetting phase plume of thickness \f$h\f$ lies below the top of a layer
 * of thickness \f$H\f$, with residual wetting phase \f$S_{wr}\f$ in the plume
 * and trapped non-wetting phase \f$S_{nr}\f$ down to the maximum plume
 * thickness \f$h_{max}\f$:
 * \f[
 * S_n H = h (1 - S_{wr}) + (h_{max} - h) S_{nr} \;.
 * \f]
 * The capillary pressure \f$(\rho_w - \rho_n) g h\f$ depends on the densities,
 * so this law returns the normalized capillary pressure \f$h/H\f$. With the
 * transmissivity fraction \f$F(z)\f$ above the relative depth \f$z\f$, the
 * relative permeabilities are
 * \f$k_{rw} = 1 - F(h_{max}/H) + k_{rw}^{end} (F(h_{max}/H) - F(h/H))\f$ and
 * \f$k_{rn}(F(h/H))\f$ of VeNonWettingRelpermModel.
 *
 * \tparam TraitsT The material traits.
 * \tparam ParamsT The parameter class.
 *
 * \see VeSharpInterfaceLawParams
 */
template <class TraitsT, class ParamsT = VeSharpInterfaceLawParams<TraitsT>>
class VeSharpInterfaceLaw : public TraitsT
{
public:
    using Traits = TraitsT; //!< The material traits
    using Params = ParamsT; //!< The parameter class
    using Scalar = typename Traits::Scalar; //!< The floating point type

    //! The number of fluid phases to which this material law applies.
    static constexpr int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The vertical-equilibrium sharp-interface law only applies "
                  "to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static constexpr bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static constexpr bool implementsTwoPhaseSatApi = true;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static constexpr bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static constexpr bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static constexpr bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static constexpr bool isCompositionDependent = false;

    /*!
     * \brief The normalized capillary pressures.
     *
     * \param values The capillary pressures of the phases.
     * \param params The material law parameters.
     * \param fs The fluid state.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        using Evaluation = std::remove_reference_t<decltype(values[0])>;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The saturations from the normalized pressure difference.
     *
     * \param values The saturations of the phases.
     * \param params The material law parameters.
     * \param fs The fluid state.
     */
    template <class Container, class FluidState>
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
        using Evaluation = std::remove_reference_t<decltype(values[0])>;

        values[Traits::wettingPhaseIdx] = Sw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = 1.0 - values[Traits::wettingPhaseIdx];
    }

    /*!
     * \brief The relative permeabilities.
     *
     * \param values The relative permeabilities of the phases.
     * \param params The material law parameters.
     * \param fs The fluid state.
     */
    template <class Container, class FluidState>
    static void
    relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        using Evaluation = std::remove_reference_t<decltype(values[0])>;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The normalized capillary pressure \f$h/H\f$.
     *
     * \param params The material law parameters.
     * \param fs The fluid state.
     * \return The normalized capillary pressure.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief The normalized capillary pressure \f$h/H\f$.
     *
     * \param params The material law parameters.
     * \param Sw The wetting phase saturation.
     * \return The normalized capillary pressure.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        return plumeFraction(params, 1.0 - Sw);
    }

    /*!
     * \brief The wetting phase saturation from the normalized capillary pressure.
     *
     * \param params The material law parameters.
     * \param pcnw The normalized capillary pressure.
     * \return The wetting phase saturation.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatPcnwInv(const Params& params, const Evaluation& pcnw)
    {
        return twoPhaseSatSw(params, pcnw);
    }

    /*!
     * \brief The wetting phase saturation from the normalized pressure difference.
     *
     * \param params The material law parameters.
     * \param fs The fluid state.
     * \return The wetting phase saturation.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation Sw(const Params& params, const FluidState& fs)
    {
        const Evaluation pc = decay<Evaluation>(fs.pressure(Traits::nonWettingPhaseIdx))
            - decay<Evaluation>(fs.pressure(Traits::wettingPhaseIdx));
        return twoPhaseSatSw(params, pc);
    }

    /*!
     * \brief The wetting phase saturation from the normalized capillary pressure.
     *
     * \param params The material law parameters.
     * \param pc The normalized capillary pressure.
     * \return The wetting phase saturation.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& params, const Evaluation& pc)
    {
        return 1.0 - nonWettingSaturation(params, max(min(pc, Scalar {1.0}), Scalar {0.0}));
    }

    /*!
     * \brief The non-wetting phase saturation from the normalized pressure difference.
     *
     * \param params The material law parameters.
     * \param fs The fluid state.
     * \return The non-wetting phase saturation.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    {
        return 1.0 - Sw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The non-wetting phase saturation from the normalized capillary pressure.
     *
     * \param params The material law parameters.
     * \param pc The normalized capillary pressure.
     * \return The non-wetting phase saturation.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& params, const Evaluation& pc)
    {
        return 1.0 - twoPhaseSatSw(params, pc);
    }

    /*!
     * \brief The relative permeability of the wetting phase.
     *
     * \param params The material law parameters.
     * \param fs The fluid state.
     * \return The relative permeability.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatKrw(params, Sw);
    }

    /*!
     * \brief The relative permeability of the wetting phase.
     *
     * \param params The material law parameters.
     * \param Sw The wetting phase saturation.
     * \return The relative permeability.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {
        const auto parts = twoPhaseSatKrwParts(params, Sw);
        return min(parts[0] + parts[1], Scalar {1.0});
    }

    /*!
     * \brief The parts of the wetting phase relative permeability in the
     *        trapped zone and in the brine zone below it.
     *
     * The sum is the relative permeability before it is limited to one.
     *
     * \param params The material law parameters.
     * \param Sw The wetting phase saturation.
     * \return The parts in the trapped zone and in the brine zone.
     */
    template <class Evaluation>
    static std::array<Evaluation, 2> twoPhaseSatKrwParts(const Params& params, const Evaluation& Sw)
    {
        const Evaluation h = plumeFraction(params, 1.0 - Sw);
        Evaluation hMax = maxPlumeFraction(params, 1.0 - Sw);
        if (scalarValue(h) > scalarValue(hMax)) {
            hMax = h;
        }

        const Evaluation plume = params.transmissivityFraction(h);
        const Evaluation trapped = params.transmissivityFraction(hMax);
        Evaluation brine = 1.0 - trapped;
        if (params.nonWettingRelpermModel() == VeNonWettingRelpermModel::Linear) {
            brine /= 1.0 - params.nonWettingRelpermThreshold();
        }
        return {params.wettingRelpermEndpoint() * (trapped - plume), brine};
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     *
     * \param params The material law parameters.
     * \param fs The fluid state.
     * \return The relative permeability.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = 1.0 - decay<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));
        return twoPhaseSatKrn(params, Sw);
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     *
     * \param params The material law parameters.
     * \param Sw The wetting phase saturation.
     * \return The relative permeability.
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    {
        const Evaluation g = params.transmissivityFraction(plumeFraction(params, 1.0 - Sw));
        const Scalar a = params.nonWettingRelpermThreshold();

        Evaluation kr = 0.0;
        switch (params.nonWettingRelpermModel()) {
        case VeNonWettingRelpermModel::Tabulated: {
            const auto& table = params.nonWettingRelpermTable();
            kr = table.eval(max(min(g, table.xMax()), table.xMin()));
            break;
        }
        case VeNonWettingRelpermModel::Sinusoidal:
            if (scalarValue(g) > a) {
                kr = sqrt(g * g - a * a);
            }
            break;
        case VeNonWettingRelpermModel::Linear:
            kr = (g - a) / (1.0 - a);
            break;
        case VeNonWettingRelpermModel::LinearThreshold:
            kr = g - a;
            break;
        }
        return min(params.nonWettingRelpermEndpoint() * max(kr, Scalar {0.0}), Scalar {1.0});
    }

    /*!
     * \brief The relative thickness \f$h/H\f$ of the mobile plume.
     *
     * \param params The material law parameters.
     * \param Sn The non-wetting phase saturation.
     * \return The relative plume thickness.
     */
    template <class Evaluation>
    static Evaluation plumeFraction(const Params& params, const Evaluation& Sn)
    {
        const Scalar swr = params.residualWettingSaturation();
        const Scalar snr = params.residualNonWettingSaturation();
        const Scalar hMax = params.maxPlumeFraction();

        if (scalarValue(Sn) >= hMax * (1.0 - swr)) {
            // drainage: no residual zone below the plume
            return min(Sn / (1.0 - swr), Scalar {1.0});
        }
        return max((Sn - hMax * snr) / (1.0 - swr - snr), Scalar {0.0});
    }

    /*!
     * \brief The relative maximum plume thickness \f$h_{max}/H\f$ of the history.
     *
     * The trapped zone holds at most the non-wetting phase of the cell, so
     * \f$h_{max} \le S_n H / S_{nr}\f$. This matters if the non-wetting
     * phase is removed from the trapped zone, e.g. by dissolution.
     *
     * \param params The material law parameters.
     * \param Sn The non-wetting phase saturation.
     * \return The relative maximum plume thickness.
     */
    template <class Evaluation>
    static Evaluation maxPlumeFraction(const Params& params, const Evaluation& Sn)
    {
        const Scalar snr = params.residualNonWettingSaturation();
        const Scalar hMax = params.maxPlumeFraction();
        if (snr > 0.0 && scalarValue(Sn) < hMax * snr) {
            return max(Sn / snr, Scalar {0.0});
        }
        return Evaluation {hMax};
    }

    /*!
     * \brief The relative thickness \f$(h_{max} - h)/H\f$ of the trapped zone.
     *
     * \param params The material law parameters.
     * \param Sn The non-wetting phase saturation.
     * \return The relative thickness of the trapped zone.
     */
    template <class Evaluation>
    static Evaluation residualFraction(const Params& params, const Evaluation& Sn)
    {
        return max(maxPlumeFraction(params, Sn) - plumeFraction(params, Sn), Scalar {0.0});
    }

    /*!
     * \brief The non-wetting phase saturation of a relative plume thickness.
     *
     * \param params The material law parameters.
     * \param plumeFraction The relative plume thickness \f$h/H\f$.
     * \return The non-wetting phase saturation.
     */
    template <class Evaluation>
    static Evaluation nonWettingSaturation(const Params& params, const Evaluation& plumeFraction)
    {
        const Scalar swr = params.residualWettingSaturation();
        const Scalar snr = params.residualNonWettingSaturation();
        const Scalar hMax = params.maxPlumeFraction();

        if (scalarValue(plumeFraction) >= hMax) {
            return plumeFraction * (1.0 - swr);
        }
        return plumeFraction * (1.0 - swr - snr) + hMax * snr;
    }

    /*!
     * \brief Reduce the maximum plume thickness after dissolution.
     *
     * The trapped non-wetting phase dissolves first, then the plume.
     *
     * \param params The material law parameters to update.
     * \param SnBefore The non-wetting phase saturation before dissolution.
     * \param SnAfter The non-wetting phase saturation after dissolution.
     */
    static void depleteResidual(Params& params, Scalar SnBefore, Scalar SnAfter)
    {
        if (SnAfter >= SnBefore) {
            return;
        }

        const Scalar swr = params.residualWettingSaturation();
        const Scalar snr = params.residualNonWettingSaturation();
        const Scalar h = plumeFraction(params, SnBefore);

        Scalar hMax = 0.0;
        if (snr > 0.0) {
            hMax = (SnAfter - h * (1.0 - swr - snr)) / snr;
        }
        if (hMax < h) {
            hMax = SnAfter / (1.0 - swr);
        }
        params.setMaxPlumeFraction(std::min(hMax, params.maxPlumeFraction()));
    }
};

} // namespace Opm

#endif // OPM_VE_SHARP_INTERFACE_LAW_HPP
