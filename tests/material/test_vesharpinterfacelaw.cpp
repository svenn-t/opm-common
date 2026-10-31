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
 * \brief Unit test for the vertical-equilibrium sharp-interface material law.
 *
 * The reference values are computed by hand from the formulas in the
 * documentation of VeSharpInterfaceLaw.
 */
#include "config.h"

#define BOOST_TEST_MODULE VeSharpInterfaceLaw
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/VeSharpInterfaceLaw.hpp>

#include <array>
#include <stdexcept>
#include <vector>

namespace
{

using Scalar = double;
using Traits = Opm::TwoPhaseMaterialTraits<Scalar, /*wettingPhaseIdx=*/0, /*nonWettingPhaseIdx=*/1>;
using Law = Opm::VeSharpInterfaceLaw<Traits>;
using Params = Law::Params;
using Evaluation = Opm::DenseAd::Evaluation<Scalar, 1>;

constexpr Scalar tolerance = 1e-12;

// Minimal two-phase fluid state providing what the law needs
template <class ValueT>
struct FluidState
{
    using ValueType = ValueT;

    const ValueType& saturation(int phaseIdx) const
    {
        return saturation_[phaseIdx];
    }
    const ValueType& pressure(int phaseIdx) const
    {
        return pressure_[phaseIdx];
    }

    std::array<ValueType, 2> saturation_ {};
    std::array<ValueType, 2> pressure_ {};
};

Params
makeParams(Opm::VeNonWettingRelpermModel model = Opm::VeNonWettingRelpermModel::LinearThreshold,
           Scalar threshold = 0.0)
{
    Params params;
    params.setResidualWettingSaturation(0.2);
    params.setResidualNonWettingSaturation(0.25);
    params.setWettingRelpermEndpoint(0.3);
    params.setNonWettingRelpermEndpoint(0.8);
    params.setNonWettingRelpermModel(model);
    params.setNonWettingRelpermThreshold(threshold);
    if (model == Opm::VeNonWettingRelpermModel::Tabulated) {
        params.setNonWettingRelpermTable(std::vector<Scalar> {0.0, 0.5, 1.0},
                                         std::vector<Scalar> {0.0, 0.25, 1.0});
    }
    params.finalize();
    return params;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(Drainage)
{
    const Params params = makeParams();
    const Scalar Sn = 0.4;

    // h/H = Sn/(1 - Swr) = 0.5
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, Sn), 0.5, tolerance);
    BOOST_CHECK_SMALL(Law::residualFraction(params, Sn), tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 1.0 - Sn), 0.8 * 0.5, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, 1.0 - Sn), 0.5, tolerance);

    // the normalized capillary pressure is h/H
    const Scalar pc = 0.5;
    BOOST_CHECK_CLOSE(Law::twoPhaseSatPcnw(params, 1.0 - Sn), pc, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatSw(params, pc), 1.0 - Sn, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatPcnwInv(params, pc), 1.0 - Sn, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatSn(params, pc), Sn, tolerance);

    // a completely filled layer
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, 0.8), 1.0, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 0.2), 0.8, tolerance);
    BOOST_CHECK_SMALL(Law::twoPhaseSatKrw(params, 0.2), tolerance);
}

BOOST_AUTO_TEST_CASE(Imbibition)
{
    Params params = makeParams();
    BOOST_CHECK(params.update(0.4));
    BOOST_CHECK(!params.update(0.3));
    BOOST_CHECK_CLOSE(params.maxPlumeFraction(), 0.5, tolerance);
    BOOST_CHECK_CLOSE(params.maxNonWettingSaturation(), 0.4, tolerance);

    // h/H = (Sn - hmax/H*Snr)/(1 - Swr - Snr) = 0.125/0.55
    const Scalar Sn = 0.25;
    const Scalar h = 0.125 / 0.55;
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, Sn), h, tolerance);
    BOOST_CHECK_CLOSE(Law::residualFraction(params, Sn), 0.5 - h, tolerance);
    BOOST_CHECK_CLOSE(Law::nonWettingSaturation(params, h), Sn, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 1.0 - Sn), 0.8 * h, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, 1.0 - Sn), 0.5 + 0.3 * (0.5 - h), tolerance);

    BOOST_CHECK_CLOSE(Law::twoPhaseSatPcnw(params, 1.0 - Sn), h, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatSw(params, h), 1.0 - Sn, tolerance);

    // all non-wetting phase is trapped
    BOOST_CHECK_SMALL(Law::plumeFraction(params, 0.125), tolerance);
    BOOST_CHECK_SMALL(Law::plumeFraction(params, 0.1), tolerance);
    BOOST_CHECK_SMALL(Law::twoPhaseSatKrn(params, 1.0 - 0.125), tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, 1.0 - 0.125), 0.5 + 0.3 * 0.5, tolerance);

    // the mapping is continuous at the maximum saturation
    const Scalar SnMax = params.maxNonWettingSaturation();
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, SnMax), 0.5, tolerance);
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, SnMax - 1e-10), 0.5, 1e-6);

    // the derivative has a kink at the maximum saturation
    const Evaluation SnBelow = Evaluation::createVariable(SnMax - 1e-6, 0);
    const Evaluation SnAbove = Evaluation::createVariable(SnMax + 1e-6, 0);
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, SnBelow).derivative(0), 1.0 / 0.55, tolerance);
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, SnAbove).derivative(0), 1.0 / 0.8, tolerance);
}

BOOST_AUTO_TEST_CASE(Monotonicity)
{
    for (const Scalar SnMax : {0.0, 0.4, 0.8}) {
        Params params = makeParams();
        params.update(SnMax);

        Scalar krnOld = Law::twoPhaseSatKrn(params, 1.0);
        Scalar krwOld = Law::twoPhaseSatKrw(params, 1.0);
        Scalar pcOld = Law::twoPhaseSatPcnw(params, 1.0);
        for (int i = 1; i <= 80; ++i) {
            const Scalar Sw = 1.0 - i * 0.01;
            const Scalar krn = Law::twoPhaseSatKrn(params, Sw);
            const Scalar krw = Law::twoPhaseSatKrw(params, Sw);
            const Scalar pc = Law::twoPhaseSatPcnw(params, Sw);
            BOOST_CHECK_GE(krn, krnOld);
            BOOST_CHECK_LE(krw, krwOld);
            BOOST_CHECK_GE(pc, pcOld);
            BOOST_CHECK_GE(krw, 0.0);
            BOOST_CHECK_LE(krw, 1.0);
            BOOST_CHECK_GE(krn, 0.0);
            BOOST_CHECK_LE(krn, 1.0);
            krnOld = krn;
            krwOld = krw;
            pcOld = pc;
        }
    }
}

BOOST_AUTO_TEST_CASE(NonWettingRelpermModels)
{
    using Model = Opm::VeNonWettingRelpermModel;

    // plume fraction 0.5 and 0.05
    const Scalar Sw = 0.6;
    const Scalar SwSmall = 0.96;

    {
        const Params params = makeParams(Model::Sinusoidal, 0.1);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, Sw), 0.8 * std::sqrt(0.25 - 0.01), tolerance);
        BOOST_CHECK_SMALL(Law::twoPhaseSatKrn(params, SwSmall), tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, Sw), 0.5, tolerance);
    }
    {
        const Params params = makeParams(Model::Linear, 0.1);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, Sw), 0.8 * 0.4 / 0.9, tolerance);
        BOOST_CHECK_SMALL(Law::twoPhaseSatKrn(params, SwSmall), tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, Sw), 0.5 / 0.9, tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, 1.0), 1.0, tolerance);
    }
    {
        const Params params = makeParams(Model::LinearThreshold, 0.1);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, Sw), 0.8 * 0.4, tolerance);
        BOOST_CHECK_SMALL(Law::twoPhaseSatKrn(params, SwSmall), tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, Sw), 0.5, tolerance);
    }
    {
        const Params params = makeParams(Model::Tabulated);
        // plume fractions 0.75 and 0.25
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 0.4), 0.8 * 0.625, tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 0.8), 0.8 * 0.125, tolerance);
        BOOST_CHECK_SMALL(Law::twoPhaseSatKrn(params, 1.0), tolerance);
        BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 0.2), 0.8, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(SubLayers)
{
    Params params = makeParams();
    // thickness fractions 0.2 and 0.8, transmissivities 5 and 4
    params.setSubLayers(std::vector<Scalar> {10.0, 40.0}, std::vector<Scalar> {5e-13, 1e-13});
    params.finalize();
    BOOST_CHECK(params.hasSubLayers());
    BOOST_CHECK_CLOSE(params.transmissivityFraction(0.2), 5.0 / 9.0, tolerance);
    BOOST_CHECK_CLOSE(params.transmissivityFraction(1.0), 1.0, tolerance);

    // h/H = 0.5 gives F = 5/9 + 0.3/0.8*4/9 = 6.5/9
    const Scalar Sw = 0.6;
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, Sw), 0.8 * 6.5 / 9.0, tolerance);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrw(params, Sw), 2.5 / 9.0, tolerance);

    // the residual zone between h/H = 0.1 and 0.5 is weighted with the
    // transmissivity as well
    params.update(0.4);
    const Scalar Sn = Law::nonWettingSaturation(params, 0.1);
    BOOST_CHECK_CLOSE(Law::twoPhaseSatKrn(params, 1.0 - Sn), 0.8 * 2.5 / 9.0, tolerance);
    BOOST_CHECK_CLOSE(
        Law::twoPhaseSatKrw(params, 1.0 - Sn), 2.5 / 9.0 + 0.3 * 4.0 / 9.0, tolerance);
}

BOOST_AUTO_TEST_CASE(DepleteResidual)
{
    Params params = makeParams();
    params.update(0.4);

    // the residual zone shrinks while the plume stays unchanged
    const Scalar h = 0.125 / 0.55;
    Law::depleteResidual(params, 0.25, 0.2);
    BOOST_CHECK_CLOSE(params.maxPlumeFraction(), 0.3, tolerance);
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, 0.2), h, tolerance);

    // the residual zone is used up and the plume shrinks
    Law::depleteResidual(params, 0.2, 0.1);
    BOOST_CHECK_CLOSE(params.maxPlumeFraction(), 0.125, tolerance);
    BOOST_CHECK_CLOSE(Law::plumeFraction(params, 0.1), 0.125, tolerance);
    BOOST_CHECK_SMALL(Law::residualFraction(params, 0.1), tolerance);

    // no dissolution keeps the history
    Law::depleteResidual(params, 0.1, 0.1);
    BOOST_CHECK_CLOSE(params.maxPlumeFraction(), 0.125, tolerance);
}

BOOST_AUTO_TEST_CASE(FluidStateApi)
{
    Params params = makeParams();
    params.update(0.4);

    FluidState<Evaluation> fs;
    const Scalar Sn = 0.25;
    fs.saturation_[1] = Evaluation::createVariable(Sn, 0);
    fs.saturation_[0] = 1.0 - fs.saturation_[1];

    // the normalized capillary pressure is h/H
    std::array<Evaluation, 2> pc;
    Law::capillaryPressures(pc, params, fs);
    const Scalar h = 0.125 / 0.55;
    BOOST_CHECK_SMALL(pc[0].value(), tolerance);
    BOOST_CHECK_CLOSE(pc[1].value(), h, tolerance);
    BOOST_CHECK_CLOSE(pc[1].derivative(0), 1.0 / 0.55, tolerance);

    std::array<Evaluation, 2> kr;
    Law::relativePermeabilities(kr, params, fs);
    BOOST_CHECK_CLOSE(kr[0].value(), Law::twoPhaseSatKrw(params, 1.0 - Sn), tolerance);
    BOOST_CHECK_CLOSE(kr[1].value(), Law::twoPhaseSatKrn(params, 1.0 - Sn), tolerance);
    BOOST_CHECK_CLOSE(kr[1].derivative(0), 0.8 / 0.55, tolerance);
    BOOST_CHECK_CLOSE(kr[0].derivative(0), -0.3 / 0.55, tolerance);

    fs.pressure_[0] = 0.0;
    fs.pressure_[1] = pc[1].value();
    std::array<Evaluation, 2> S;
    Law::saturations(S, params, fs);
    BOOST_CHECK_CLOSE(S[1].value(), Sn, tolerance);
}

BOOST_AUTO_TEST_CASE(InvalidParameters)
{
    Params params = makeParams();
    params.setResidualNonWettingSaturation(0.8);
    BOOST_CHECK_THROW(params.finalize(), std::invalid_argument);

    params = makeParams();
    params.setNonWettingRelpermModel(Opm::VeNonWettingRelpermModel::Tabulated);
    BOOST_CHECK_THROW(params.finalize(), std::invalid_argument);

    params = makeParams();
    BOOST_CHECK_THROW(params.setSubLayers(std::vector<Scalar> {1.0}, std::vector<Scalar> {}),
                      std::invalid_argument);
}
