/*
Copyright 2023 NORCE.

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
#include <opm/input/eclipse/EclipseState/BactPara.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/B.hpp>

Opm::BactPara::BactPara() :
    m_max_growth_rate( ParserKeywords::BACTPARA::MAX_GROWTH_RATE::defaultValue ),
    m_half_velocity_gas( ParserKeywords::BACTPARA::HALF_VELOCITY_GAS::defaultValue ),
    m_yield_coefficient( ParserKeywords::BACTPARA::YIELD_COEFFICIENT::defaultValue ),
    m_stoichiometric_coefficient( ParserKeywords::BACTPARA::STOICHIOMETRIC_COEFFICIENT::defaultValue ),
    m_decay_coefficient( ParserKeywords::BACTPARA::DECAY_COEFFICIENT::defaultValue )
{}

Opm::BactPara::BactPara( const Opm::Deck& deck ) 
    : BactPara()
{
    using namespace Opm::ParserKeywords;

    if (!deck.hasKeyword<BACTPARA>())
        return;
    
    const auto& keyword = deck.get<BACTPARA>().back();
    const auto& record = keyword.getRecord(0);
    this->m_max_growth_rate = record.getItem<BACTPARA::MAX_GROWTH_RATE>().get< double >(0);
    this->m_half_velocity_gas = record.getItem<BACTPARA::HALF_VELOCITY_GAS>().get< double >(0);
    this->m_yield_coefficient = record.getItem<BACTPARA::YIELD_COEFFICIENT>().get< double >(0);
    this->m_stoichiometric_coefficient = record.getItem<BACTPARA::STOICHIOMETRIC_COEFFICIENT>().get< int >(0);
    this->m_decay_coefficient = record.getItem<BACTPARA::DECAY_COEFFICIENT>().get< double >(0);
}