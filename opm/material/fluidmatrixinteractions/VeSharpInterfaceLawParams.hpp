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
 * \copydoc Opm::VeSharpInterfaceLawParams
 */
#ifndef OPM_VE_SHARP_INTERFACE_LAW_PARAMS_HPP
#define OPM_VE_SHARP_INTERFACE_LAW_PARAMS_HPP

#include <opm/material/common/EnsureFinalized.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace Opm {

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Model for the vertically upscaled relative permeability of the
 *        non-wetting phase in the plume.
 *
 * The argument of all models is the plume fraction \f$g\f$, i.e., the
 * permeability-weighted fraction of the layer thickness occupied by the mobile
 * plume, and \f$a\f$ is the threshold set by
 * VeSharpInterfaceLawParams::setNonWettingRelpermThreshold().
 */
enum class VeNonWettingRelpermModel {
    //! \f$k_{rn} = k_{rn}^{end} f(g)\f$ with a tabulated function \f$f\f$
    Tabulated,
    //! \f$k_{rn} = k_{rn}^{end} \sqrt{g^2 - a^2}\f$ for \f$g > a\f$
    Sinusoidal,
    //! \f$k_{rn} = k_{rn}^{end} (g - a)/(1 - a)\f$ for \f$g > a\f$, mobile
    //! wetting phase scaled by \f$1/(1 - a)\f$
    Linear,
    //! \f$k_{rn} = k_{rn}^{end} (g - a)\f$ for \f$g > a\f$
    LinearThreshold,
};

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the
 *        vertical-equilibrium sharp-interface law.
 *
 * Besides the static parameters this class holds the history of the
 * cell, i.e., the maximum relative plume thickness \f$\hat{h}_{max} =
 * h_{max}/H\f$ which has been reached so far. Below the mobile plume and
 * down to \f$h_{max}\f$ the non-wetting phase is trapped at its residual
 * saturation.
 *
 *\see VeSharpInterfaceLaw
 */
template <class TraitsT>
class VeSharpInterfaceLawParams : public EnsureFinalized
{
    using Scalar = typename TraitsT::Scalar;

public:
    using EnsureFinalized::finalize;

    using Traits = TraitsT;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        if (residualWettingSaturation_ < 0.0 || residualNonWettingSaturation_ < 0.0
            || residualWettingSaturation_ + residualNonWettingSaturation_ >= 1.0)
        {
            throw std::invalid_argument("The residual saturations must be non-negative "
                                        "and sum up to less than one");
        }
        if (nonWettingRelpermThreshold_ < 0.0 || nonWettingRelpermThreshold_ >= 1.0) {
            throw std::invalid_argument("The relative permeability threshold must be in [0, 1)");
        }
        if (nonWettingRelpermModel_ == VeNonWettingRelpermModel::Tabulated
            && nonWettingRelpermTable_.numSamples() < 2)
        {
            throw std::invalid_argument("The tabulated relative permeability model "
                                        "requires a table with at least two entries");
        }

        EnsureFinalized::finalize();
    }

    /*!
     * \brief Returns the residual saturation of the wetting phase in the plume.
     */
    Scalar residualWettingSaturation() const
    { EnsureFinalized::check(); return residualWettingSaturation_; }

    /*!
     * \brief Set the residual saturation of the wetting phase in the plume.
     */
    void setResidualWettingSaturation(Scalar value)
    { residualWettingSaturation_ = value; }

    /*!
     * \brief Returns the residual saturation of the non-wetting phase below the plume.
     */
    Scalar residualNonWettingSaturation() const
    { EnsureFinalized::check(); return residualNonWettingSaturation_; }

    /*!
     * \brief Set the residual saturation of the non-wetting phase below the plume.
     */
    void setResidualNonWettingSaturation(Scalar value)
    { residualNonWettingSaturation_ = value; }

    /*!
     * \brief Returns the wetting phase relative permeability in the residual zone.
     */
    Scalar wettingRelpermEndpoint() const
    { EnsureFinalized::check(); return wettingRelpermEndpoint_; }

    /*!
     * \brief Set the wetting phase relative permeability in the residual zone.
     */
    void setWettingRelpermEndpoint(Scalar value)
    { wettingRelpermEndpoint_ = value; }

    /*!
     * \brief Returns the non-wetting phase relative permeability in the plume.
     */
    Scalar nonWettingRelpermEndpoint() const
    { EnsureFinalized::check(); return nonWettingRelpermEndpoint_; }

    /*!
     * \brief Set the non-wetting phase relative permeability in the plume.
     */
    void setNonWettingRelpermEndpoint(Scalar value)
    { nonWettingRelpermEndpoint_ = value; }

    /*!
     * \brief Returns the model for the non-wetting phase relative permeability.
     */
    VeNonWettingRelpermModel nonWettingRelpermModel() const
    { EnsureFinalized::check(); return nonWettingRelpermModel_; }

    /*!
     * \brief Set the model for the non-wetting phase relative permeability.
     */
    void setNonWettingRelpermModel(VeNonWettingRelpermModel value)
    { nonWettingRelpermModel_ = value; }

    /*!
     * \brief Returns the plume fraction below which the non-wetting phase is immobile.
     */
    Scalar nonWettingRelpermThreshold() const
    { EnsureFinalized::check(); return nonWettingRelpermThreshold_; }

    /*!
     * \brief Set the plume fraction below which the non-wetting phase is immobile.
     */
    void setNonWettingRelpermThreshold(Scalar value)
    { nonWettingRelpermThreshold_ = value; }

    /*!
     * \brief Returns the tabulated non-wetting phase relative permeability
     *        as a function of the plume fraction.
     */
    const TabulatedFunction& nonWettingRelpermTable() const
    { EnsureFinalized::check(); return nonWettingRelpermTable_; }

    /*!
     * \brief Set the tabulated non-wetting phase relative permeability
     *        as a function of the plume fraction.
     *
     * The values are multiplied with the end point relative permeability.
     * Plume fractions outside of the table range are clamped to the table.
     */
    template <class Container>
    void setNonWettingRelpermTable(const Container& plumeFraction,
                                   const Container& relperm)
    { nonWettingRelpermTable_.setXYContainers(plumeFraction, relperm); }

    /*!
     * \brief Returns true if the layer consists of sub-layers with different
     *        permeabilities.
     */
    bool hasSubLayers() const
    { EnsureFinalized::check(); return subLayerTransmissivity_.numSamples() > 0; }

    /*!
     * \brief Returns the fraction of the layer transmissivity located above a
     *        given relative depth below the top of the layer.
     *
     * The relative depth is the depth divided by the layer thickness.
     */
    template <class Evaluation>
    Evaluation transmissivityFraction(const Evaluation& relativeDepth) const
    {
        EnsureFinalized::check();
        if (!hasSubLayers()) {
            return relativeDepth;
        }
        return subLayerTransmissivity_.eval(relativeDepth, /*extrapolate=*/true);
    }

    /*!
     * \brief Set the sub-layers of the layer.
     *
     * \param thickness The thicknesses of the sub-layers ordered from top to
     *                  bottom. They are normalized by their sum.
     * \param permeability The permeabilities of the sub-layers.
     */
    template <class Container>
    void setSubLayers(const Container& thickness, const Container& permeability)
    {
        if (thickness.size() != permeability.size() || thickness.empty()) {
            throw std::invalid_argument("Each sub-layer requires a thickness and a permeability");
        }

        std::vector<Scalar> depth(thickness.size() + 1, 0.0);
        std::vector<Scalar> transmissivity(thickness.size() + 1, 0.0);
        for (std::size_t i = 0; i < thickness.size(); ++i) {
            if (thickness[i] <= 0.0 || permeability[i] <= 0.0) {
                throw std::invalid_argument("The sub-layer thicknesses and "
                                            "permeabilities must be positive");
            }
            depth[i + 1] = depth[i] + thickness[i];
            transmissivity[i + 1] = transmissivity[i] + thickness[i]*permeability[i];
        }
        const Scalar totalDepth = depth.back();
        const Scalar totalTransmissivity = transmissivity.back();
        for (std::size_t i = 0; i < depth.size(); ++i) {
            depth[i] /= totalDepth;
            transmissivity[i] /= totalTransmissivity;
        }
        subLayerTransmissivity_.setXYContainers(depth, transmissivity, /*sortInputs=*/false);
    }

    /*!
     * \brief Returns the maximum relative plume thickness reached so far.
     */
    Scalar maxPlumeFraction() const
    { return maxPlumeFraction_; }

    /*!
     * \brief Set the maximum relative plume thickness reached so far.
     */
    void setMaxPlumeFraction(Scalar value)
    { maxPlumeFraction_ = std::clamp(value, Scalar{0.0}, Scalar{1.0}); }

    /*!
     * \brief Returns the maximum non-wetting phase saturation reached so far.
     */
    Scalar maxNonWettingSaturation() const
    { return maxPlumeFraction_*(1.0 - residualWettingSaturation_); }

    /*!
     * \brief Update the history with the non-wetting phase saturation of a
     *        converged time step.
     *
     * \return True if the maximum plume thickness has changed.
     */
    bool update(Scalar Sn)
    {
        const Scalar plumeFraction = Sn/(1.0 - residualWettingSaturation_);
        if (plumeFraction <= maxPlumeFraction_) {
            return false;
        }
        setMaxPlumeFraction(plumeFraction);
        return true;
    }

private:
    Scalar residualWettingSaturation_{0.0};
    Scalar residualNonWettingSaturation_{0.0};
    Scalar wettingRelpermEndpoint_{1.0};
    Scalar nonWettingRelpermEndpoint_{1.0};
    VeNonWettingRelpermModel nonWettingRelpermModel_{VeNonWettingRelpermModel::LinearThreshold};
    Scalar nonWettingRelpermThreshold_{0.0};
    TabulatedFunction nonWettingRelpermTable_{};
    TabulatedFunction subLayerTransmissivity_{};
    Scalar maxPlumeFraction_{0.0};
};

} // namespace Opm

#endif // OPM_VE_SHARP_INTERFACE_LAW_PARAMS_HPP
