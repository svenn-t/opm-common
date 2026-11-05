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

namespace Opm
{

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Models for the vertically upscaled non-wetting phase relative permeability.
 *
 * \f$h\f$ is the relative plume thickness, \f$g\f$ the permeability-weighted
 * plume fraction and \f$a\f$ the threshold of
 * VeSharpInterfaceLawParams::setNonWettingRelpermThreshold().
 */
enum class VeNonWettingRelpermModel {
    //! \f$k_{rn} = k_{rn}^{end} f(h)\f$ with a tabulated upscaled function \f$f\f$
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
 * \brief Material parameters of the vertical-equilibrium sharp-interface law.
 *
 * Also holds the history of the cell, i.e., the maximum relative plume
 * thickness \f$h_{max}/H\f$ reached so far.
 *
 * \tparam TraitsT The material traits.
 *
 * \see VeSharpInterfaceLaw
 */
template <class TraitsT>
class VeSharpInterfaceLawParams : public EnsureFinalized
{
    using Scalar = typename TraitsT::Scalar;

public:
    using EnsureFinalized::finalize;

    using Traits = TraitsT; //!< The material traits
    using TabulatedFunction = Tabulated1DFunction<Scalar>; //!< Type of the tables

    /*!
     * \brief Check the parameters and mark the object as finalized.
     *
     * \throw std::invalid_argument If a parameter is out of range.
     */
    void finalize()
    {
        if (residualWettingSaturation_ < 0.0 || residualNonWettingSaturation_ < 0.0
            || residualWettingSaturation_ + residualNonWettingSaturation_ >= 1.0) {
            throw std::invalid_argument("The residual saturations must be non-negative "
                                        "and sum up to less than one");
        }
        if (nonWettingRelpermThreshold_ < 0.0 || nonWettingRelpermThreshold_ >= 1.0) {
            throw std::invalid_argument("The relative permeability threshold must be in [0, 1)");
        }
        if (nonWettingRelpermModel_ == VeNonWettingRelpermModel::Tabulated
            && nonWettingRelpermTable_.numSamples() < 2) {
            throw std::invalid_argument("The tabulated relative permeability model "
                                        "requires a table with at least two entries");
        }

        EnsureFinalized::finalize();
    }

    /*!
     * \brief Returns the residual saturation of the wetting phase in the plume.
     *
     * \return The residual saturation.
     */
    Scalar residualWettingSaturation() const
    {
        EnsureFinalized::check();
        return residualWettingSaturation_;
    }

    /*!
     * \brief Set the residual saturation of the wetting phase in the plume.
     *
     * \param value The residual saturation.
     */
    void setResidualWettingSaturation(Scalar value)
    {
        residualWettingSaturation_ = value;
    }

    /*!
     * \brief Returns the residual saturation of the non-wetting phase below the plume.
     *
     * \return The residual saturation.
     */
    Scalar residualNonWettingSaturation() const
    {
        EnsureFinalized::check();
        return residualNonWettingSaturation_;
    }

    /*!
     * \brief Set the residual saturation of the non-wetting phase below the plume.
     *
     * \param value The residual saturation.
     */
    void setResidualNonWettingSaturation(Scalar value)
    {
        residualNonWettingSaturation_ = value;
    }

    /*!
     * \brief Returns the wetting phase relative permeability in the residual zone.
     *
     * \return The relative permeability.
     */
    Scalar wettingRelpermEndpoint() const
    {
        EnsureFinalized::check();
        return wettingRelpermEndpoint_;
    }

    /*!
     * \brief Set the wetting phase relative permeability in the residual zone.
     *
     * \param value The relative permeability.
     */
    void setWettingRelpermEndpoint(Scalar value)
    {
        wettingRelpermEndpoint_ = value;
    }

    /*!
     * \brief Returns the non-wetting phase relative permeability in the plume.
     *
     * \return The relative permeability.
     */
    Scalar nonWettingRelpermEndpoint() const
    {
        EnsureFinalized::check();
        return nonWettingRelpermEndpoint_;
    }

    /*!
     * \brief Set the non-wetting phase relative permeability in the plume.
     *
     * \param value The relative permeability.
     */
    void setNonWettingRelpermEndpoint(Scalar value)
    {
        nonWettingRelpermEndpoint_ = value;
    }

    /*!
     * \brief Returns the model for the non-wetting phase relative permeability.
     *
     * \return The model.
     */
    VeNonWettingRelpermModel nonWettingRelpermModel() const
    {
        EnsureFinalized::check();
        return nonWettingRelpermModel_;
    }

    /*!
     * \brief Set the model for the non-wetting phase relative permeability.
     *
     * \param value The model.
     */
    void setNonWettingRelpermModel(VeNonWettingRelpermModel value)
    {
        nonWettingRelpermModel_ = value;
    }

    /*!
     * \brief Returns the plume fraction below which the non-wetting phase is immobile.
     *
     * \return The threshold plume fraction.
     */
    Scalar nonWettingRelpermThreshold() const
    {
        EnsureFinalized::check();
        return nonWettingRelpermThreshold_;
    }

    /*!
     * \brief Set the plume fraction below which the non-wetting phase is immobile.
     *
     * \param value The threshold plume fraction in [0, 1).
     */
    void setNonWettingRelpermThreshold(Scalar value)
    {
        nonWettingRelpermThreshold_ = value;
    }

    /*!
     * \brief Returns the table of the non-wetting phase relative permeability.
     *
     * \return The relative permeability as a function of the plume fraction.
     */
    const TabulatedFunction& nonWettingRelpermTable() const
    {
        EnsureFinalized::check();
        return nonWettingRelpermTable_;
    }

    /*!
     * \brief Set the table of the non-wetting phase relative permeability.
     *
     * \param plumeFraction The plume fractions, clamped to the table range.
     * \param relperm The relative permeabilities, scaled by the end point.
     */
    template <class Container>
    void setNonWettingRelpermTable(const Container& plumeFraction, const Container& relperm)
    {
        nonWettingRelpermTable_.setXYContainers(plumeFraction, relperm);
    }

    /*!
     * \brief Returns whether the layer consists of sub-layers.
     *
     * \return True if sub-layers are set.
     */
    bool hasSubLayers() const
    {
        EnsureFinalized::check();
        return subLayerTransmissivity_.numSamples() > 0;
    }

    /*!
     * \brief Returns the fraction of the layer transmissivity above a depth.
     *
     * \param relativeDepth The depth below the top divided by the layer thickness.
     * \return The transmissivity fraction.
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
     * \param thickness The thicknesses from top to bottom, normalized by their sum.
     * \param permeability The permeabilities of the sub-layers.
     * \throw std::invalid_argument If the sizes differ or a value is not positive.
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
            transmissivity[i + 1] = transmissivity[i] + thickness[i] * permeability[i];
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
     *
     * \return The maximum plume thickness divided by the layer thickness.
     */
    Scalar maxPlumeFraction() const
    {
        return maxPlumeFraction_;
    }

    /*!
     * \brief Set the maximum relative plume thickness reached so far.
     *
     * \param value The maximum plume thickness divided by the layer thickness.
     */
    void setMaxPlumeFraction(Scalar value)
    {
        maxPlumeFraction_ = std::clamp(value, Scalar {0.0}, Scalar {1.0});
    }

    /*!
     * \brief Returns the maximum non-wetting phase saturation reached so far.
     *
     * \return The maximum saturation.
     */
    Scalar maxNonWettingSaturation() const
    {
        return maxPlumeFraction_ * (1.0 - residualWettingSaturation_);
    }

    /*!
     * \brief Update the history with the saturation of a converged time step.
     *
     * The maximum plume thickness grows with the plume and shrinks if the
     * trapped zone holds more than the non-wetting phase of the cell.
     *
     * \param Sn The non-wetting phase saturation.
     * \return True if the maximum plume thickness has changed.
     */
    bool update(Scalar Sn)
    {
        // the trapped zone holds at most the non-wetting phase of the cell
        if (Sn < maxPlumeFraction_ * residualNonWettingSaturation_) {
            setMaxPlumeFraction(Sn / residualNonWettingSaturation_);
            return true;
        }

        const Scalar plumeFraction = Sn / (1.0 - residualWettingSaturation_);
        if (plumeFraction <= maxPlumeFraction_) {
            return false;
        }
        setMaxPlumeFraction(plumeFraction);
        return true;
    }

private:
    Scalar residualWettingSaturation_ {0.0}; //!< Residual wetting saturation
    Scalar residualNonWettingSaturation_ {0.0}; //!< Residual non-wetting saturation
    Scalar wettingRelpermEndpoint_ {1.0}; //!< Wetting relperm in the residual zone
    Scalar nonWettingRelpermEndpoint_ {1.0}; //!< Non-wetting relperm in the plume
    //! Non-wetting relperm model
    VeNonWettingRelpermModel nonWettingRelpermModel_ {VeNonWettingRelpermModel::LinearThreshold};
    Scalar nonWettingRelpermThreshold_ {0.0}; //!< Plume fraction threshold
    TabulatedFunction nonWettingRelpermTable_ {}; //!< Tabulated non-wetting relperm
    TabulatedFunction subLayerTransmissivity_ {}; //!< Transmissivity fraction over depth
    Scalar maxPlumeFraction_ {0.0}; //!< History: maximum relative plume thickness
};

} // namespace Opm

#endif // OPM_VE_SHARP_INTERFACE_LAW_PARAMS_HPP
