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
#include <type_traits>

namespace Opm
{

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Vertically upscaled capillary pressure and relative permeabilities
 *        of the sharp-interface vertical-equilibrium (VE) model.
 *
 * The non-wetting phase (CO2) forms a plume of thickness \f$h\f$ below the
 * top of a layer of thickness \f$H\f$. The wetting phase (brine) in the plume
 * is at its residual saturation \f$S_{wr}\f$. Between the plume and the
 * maximum plume thickness \f$h_{max}\f$ reached so far, the non-wetting phase
 * is trapped at its residual saturation \f$S_{nr}\f$. Below that, the layer
 * is filled with the wetting phase. The vertically averaged non-wetting
 * phase saturation thus is
 * \f[
 * S_n H = h (1 - S_{wr}) + (h_{max} - h) S_{nr} \;.
 * \f]
 *
 * The phases are in vertical equilibrium, so the capillary pressure at the
 * top of the layer is
 * \f[
 * p_{c} = (\rho_w - \rho_n) g H \hat{h} \;.
 * \f]
 * Only the relative plume thickness \f$\hat{h} = h/H\f$ is a saturation
 * function. Thus, the capillary pressure of this law is the normalized
 * capillary pressure \f$\hat{h}\f$, i.e., the capillary pressure divided by
 * \f$(\rho_w - \rho_n) g H\f$. The model scales it with the phase densities,
 * the gravity and the layer thickness. Likewise, the inverse functions expect
 * the normalized capillary pressure.
 *
 * The relative permeabilities are thickness averages weighted by the
 * permeability of the sub-layers (if any). With \f$F(z)\f$ the fraction of
 * the layer transmissivity above the relative depth \f$z\f$, \f$\hat{h} =
 * h/H\f$ and \f$\hat{h}_{max} = h_{max}/H\f$, the wetting phase relative
 * permeability is
 * \f[
 * k_{rw} = 1 - F(\hat{h}_{max}) + k_{rw}^{end} (F(\hat{h}_{max}) - F(\hat{h}))
 * \f]
 * and the non-wetting phase relative permeability is a function of the
 * plume fraction \f$g = F(\hat{h})\f$ given by VeNonWettingRelpermModel.
 * Both are limited to one.
 *
 * Directional relative permeabilities are obtained by using one parameter
 * object per direction.
 *
 *\see VeSharpInterfaceLawParams
 */
template <class TraitsT, class ParamsT = VeSharpInterfaceLawParams<TraitsT>>
class VeSharpInterfaceLaw : public TraitsT
{
public:
    using Traits = TraitsT;
    using Params = ParamsT;
    using Scalar = typename Traits::Scalar;

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
     * \brief The normalized capillary pressure-saturation curves.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        using Evaluation = std::remove_reference_t<decltype(values[0])>;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their normalized pressure difference.
     */
    template <class Container, class FluidState>
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
        using Evaluation = std::remove_reference_t<decltype(values[0])>;

        values[Traits::wettingPhaseIdx] = Sw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = 1.0 - values[Traits::wettingPhaseIdx];
    }

    /*!
     * \brief The relative permeability-saturation curves.
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
     * \brief The normalized capillary pressure at the top of the layer, i.e.,
     *        the relative plume thickness \f$h/H\f$.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatPcnw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        return plumeFraction(params, 1.0 - Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnwInv(const Params& params, const Evaluation& pcnw)
    {
        return twoPhaseSatSw(params, pcnw);
    }

    /*!
     * \brief The wetting phase saturation from the normalized pressure
     *        difference of the fluid state.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation Sw(const Params& params, const FluidState& fs)
    {
        const Evaluation pc = decay<Evaluation>(fs.pressure(Traits::nonWettingPhaseIdx))
            - decay<Evaluation>(fs.pressure(Traits::wettingPhaseIdx));
        return twoPhaseSatSw(params, pc);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& params, const Evaluation& pc)
    {
        return 1.0 - nonWettingSaturation(params, max(min(pc, Scalar {1.0}), Scalar {0.0}));
    }

    /*!
     * \brief The non-wetting phase saturation from the normalized pressure
     *        difference of the fluid state.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    {
        return 1.0 - Sw<FluidState, Evaluation>(params, fs);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& params, const Evaluation& pc)
    {
        return 1.0 - twoPhaseSatSw(params, pc);
    }

    /*!
     * \brief The relative permeability of the wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {
        const Evaluation h = plumeFraction(params, 1.0 - Sw);
        Evaluation hMax = h;
        if (scalarValue(h) < params.maxPlumeFraction()) {
            hMax = params.maxPlumeFraction();
        }

        const Evaluation plume = params.transmissivityFraction(h);
        const Evaluation trapped = params.transmissivityFraction(hMax);
        Evaluation mobile = 1.0 - trapped;
        if (params.nonWettingRelpermModel() == VeNonWettingRelpermModel::Linear) {
            mobile /= 1.0 - params.nonWettingRelpermThreshold();
        }
        return min(mobile + params.wettingRelpermEndpoint() * (trapped - plume), Scalar {1.0});
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation Sw = 1.0 - decay<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));
        return twoPhaseSatKrn(params, Sw);
    }

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
     * \brief The relative thickness \f$(h_{max} - h)/H\f$ of the zone with
     *        trapped non-wetting phase below the plume.
     */
    template <class Evaluation>
    static Evaluation residualFraction(const Params& params, const Evaluation& Sn)
    {
        return max(params.maxPlumeFraction() - plumeFraction(params, Sn), Scalar {0.0});
    }

    /*!
     * \brief The non-wetting phase saturation for a given relative plume
     *        thickness \f$h/H\f$ and the history of the parameter object.
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
     * \brief Reduce the maximum plume thickness after non-wetting phase has
     *        been removed by dissolution.
     *
     * The trapped non-wetting phase is removed first, i.e., the plume
     * thickness is kept unchanged while the residual zone shrinks. Once the
     * residual zone is used up, the plume itself shrinks.
     *
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
