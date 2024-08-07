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
 * \copydoc Opm::PiecewiseLinearTwoPhaseMaterialParams
 */
#ifndef OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <algorithm>
#include <cassert>
#include <vector>
#include <stdexcept>

#include <opm/common/ErrorMacros.hpp>
#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm
{
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for a two-phase material law which
 *        uses a table and piecewise constant interpolation.
 */
template <class TraitsT, class VectorT = std::vector<typename TraitsT::Scalar>>
class PiecewiseLinearTwoPhaseMaterialParams : public EnsureFinalized
{
    using Scalar = typename TraitsT::Scalar;

public:
    using ValueVector = VectorT;

    using Traits = TraitsT;

    PiecewiseLinearTwoPhaseMaterialParams()
    {
    }

    PiecewiseLinearTwoPhaseMaterialParams(ValueVector SwPcwnSamples,
                                          ValueVector pcwnSamples,
                                          ValueVector SwKrwSamples,
                                          ValueVector krwSamples,
                                          ValueVector SwKrnSamples,
                                          ValueVector krnSamples)
        : SwPcwnSamples_(SwPcwnSamples)
        , SwKrwSamples_(SwKrwSamples)
        , SwKrnSamples_(SwKrnSamples)
        , pcwnSamples_(pcwnSamples)
        , krwSamples_(krwSamples)
        , krnSamples_(krnSamples)
    {
        finalize();
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        EnsureFinalized ::finalize();

        // revert the order of the sampling points if they were given
        // in reverse direction.
        if (SwPcwnSamples_.front() > SwPcwnSamples_.back())
            swapOrderIfPossibleThrowOtherwise_(SwPcwnSamples_, pcwnSamples_);

        if (SwKrwSamples_.front() > SwKrwSamples_.back())
            swapOrderIfPossibleThrowOtherwise_(SwKrwSamples_, krwSamples_);

        if (SwKrnSamples_.front() > SwKrnSamples_.back())
            swapOrderIfPossibleThrowOtherwise_(SwKrnSamples_, krnSamples_);
    }

    /*!
     * \brief Return the wetting-phase saturation values of all sampling points.
     */
    const ValueVector& SwKrwSamples() const
    {
        EnsureFinalized::check();
        return SwKrwSamples_;
    }

    /*!
     * \brief Return the wetting-phase saturation values of all sampling points.
     */
    const ValueVector& SwKrnSamples() const
    {
        EnsureFinalized::check();
        return SwKrnSamples_;
    }

    /*!
     * \brief Return the wetting-phase saturation values of all sampling points.
     */
    const ValueVector& SwPcwnSamples() const
    {
        EnsureFinalized::check();
        return SwPcwnSamples_;
    }

    /*!
     * \brief Return the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& pcwnSamples() const
    {
        EnsureFinalized::check();
        return pcwnSamples_;
    }

    /*!
     * \brief Set the sampling points for the capillary pressure curve.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    template <class Container>
    void setPcnwSamples(const Container& SwValues, const Container& values)
    {
        assert(SwValues.size() == values.size());

        size_t n = SwValues.size();
        SwPcwnSamples_.resize(n);
        pcwnSamples_.resize(n);

        std::copy(SwValues.begin(), SwValues.end(), SwPcwnSamples_.begin());
        std::copy(values.begin(), values.end(), pcwnSamples_.begin());
    }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& krwSamples() const
    {
        EnsureFinalized::check();
        return krwSamples_;
    }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    template <class Container>
    void setKrwSamples(const Container& SwValues, const Container& values)
    {
        assert(SwValues.size() == values.size());

        size_t n = SwValues.size();
        SwKrwSamples_.resize(n);
        krwSamples_.resize(n);

        std::copy(SwValues.begin(), SwValues.end(), SwKrwSamples_.begin());
        std::copy(values.begin(), values.end(), krwSamples_.begin());
    }

    /*!
     * \brief Return the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    const ValueVector& krnSamples() const
    {
        EnsureFinalized::check();
        return krnSamples_;
    }

    /*!
     * \brief Set the sampling points for the relative permeability
     *        curve of the non-wetting phase.
     *
     * This curve is assumed to depend on the wetting phase saturation
     */
    template <class Container>
    void setKrnSamples(const Container& SwValues, const Container& values)
    {
        assert(SwValues.size() == values.size());

        size_t n = SwValues.size();
        SwKrnSamples_.resize(n);
        krnSamples_.resize(n);

        std::copy(SwValues.begin(), SwValues.end(), SwKrnSamples_.begin());
        std::copy(values.begin(), values.end(), krnSamples_.begin());
    }

private:
    void swapOrderIfPossibleThrowOtherwise_(ValueVector& swValues, ValueVector& values) const
    {
        // TODO: comparing saturation values to the actual values we sample from looks strange
        // TODO: yet changing to swValues.back() breaks tests
        if (swValues.front() > values.back()) {
            // Reverting the order involves swapping which only works for non-consts.
            // The const expr ensures we can create constant parameter views.
            if constexpr (!std::is_const_v<typename ValueVector::value_type> && !std::is_const_v<ValueVector>) {
                for (unsigned origSampleIdx = 0; origSampleIdx < swValues.size() / 2; ++origSampleIdx) {
                    size_t newSampleIdx = swValues.size() - origSampleIdx - 1;

                    std::swap(swValues[origSampleIdx], swValues[newSampleIdx]);
                    std::swap(values[origSampleIdx], values[newSampleIdx]);
                }
            }
            else{
                OPM_THROW(std::logic_error, "Saturation values in interpolation table provided in wrong order, but table is immutable");
            }
        }
    }

    ValueVector SwPcwnSamples_;
    ValueVector SwKrwSamples_;
    ValueVector SwKrnSamples_;
    ValueVector pcwnSamples_;
    ValueVector krwSamples_;
    ValueVector krnSamples_;
};
} // namespace Opm

#endif
