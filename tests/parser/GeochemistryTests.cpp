// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 Equinor ASA.

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
#define BOOST_TEST_MODULE GeochemistryTests

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/SpeciesConfig.hpp>

using namespace Opm;

static Deck createGeochemDeck()
{
    return Parser{}.parseString(R"(
        RUNSPEC
        GEOCHEM
        test.json /
        )");
}

static Deck createSpeciesDeck()
{
    return Parser{}.parseString(R"(
        DIMENS
        3 3 3/
        TABDIMS
        /
        EQLDIMS
        /

        GRID
        DX
        27*1.0 /
        DY
        27*1.0 /
        DZ
        27*1.0 /
        TOPS
        9*10.0 /

        PROPS
        SPECIES
        CA NA /

        SOLUTION
        SBLKCA
        27*1e-3 /
        SVDPNA
        5.0  1e-4
        10.0 1e-4 /
        )");
}

BOOST_AUTO_TEST_CASE(GeochemDeck) {
    const auto deck = createGeochemDeck();
    Runspec runspec(deck);

    const auto& geochem = runspec.geochem();
    std::string file_name = "test.json";
    BOOST_CHECK(geochem.enabled());
    BOOST_CHECK_EQUAL(geochem.geochem_file_name(), file_name);
}

BOOST_AUTO_TEST_CASE(SpeciesConfigDeck) {
    const auto deck = createSpeciesDeck();
    EclipseState state(deck);

    const SpeciesConfig& species = state.species();
    BOOST_CHECK_EQUAL(species.size(), 2U);

    const auto& ca_species = species["CA"];
    BOOST_CHECK_EQUAL(ca_species.name, "CA");
    BOOST_CHECK(ca_species.concentration.has_value());
    BOOST_CHECK(!ca_species.svdp.has_value());
    for (const auto& elem : ca_species.concentration.value()) {
        BOOST_CHECK_CLOSE(elem, 1e-3, 1e-5);
    }

    const auto& na_species = species[1];
    BOOST_CHECK_EQUAL(na_species.name, "NA");
    BOOST_CHECK(!na_species.concentration.has_value());
    BOOST_CHECK(na_species.svdp.has_value());
    BOOST_CHECK_EQUAL(na_species.svdp.value().numColumns(), 2U);
}

