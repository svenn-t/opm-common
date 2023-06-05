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
#ifndef OPM_PARSER_BACTPARA_HPP
#define	OPM_PARSER_BACTPARA_HPP

namespace Opm {

class Deck;

class BactPara {
public:

    BactPara();
    explicit BactPara(const Deck& deck);

    static BactPara serializationTestObject();

    double getMaxGrowthRate() const {
        return m_max_growth_rate;
    }

    double getHalfVelocityGas() const {
        return m_half_velocity_gas;
    }

    double getYieldCoefficient() const {
        return m_yield_coefficient;
    }

    double getStoichiometricCoefficient() const {
        return m_stoichiometric_coefficient;
    }

    double getDecayCoefficient() const {
        return m_decay_coefficient;
    }

    bool operator==(const BactPara& data) const
    {
        return this->getMaxGrowthRate() == data.getMaxGrowthRate() &&
               this->getHalfVelocityGas() == data.getHalfVelocityGas() &&
               this->getYieldCoefficient() == data.getYieldCoefficient() &&
               this->getStoichiometricCoefficient() == data.getStoichiometricCoefficient() &&
               this->getDecayCoefficient() == data.getDecayCoefficient();
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(m_max_growth_rate);
        serializer(m_half_velocity_gas);
        serializer(m_yield_coefficient);
        serializer(m_stoichiometric_coefficient);
        serializer(m_decay_coefficient);
    }

private:
    double m_max_growth_rate;
    double m_half_velocity_gas;
    double m_yield_coefficient;
    double m_stoichiometric_coefficient;
    double m_decay_coefficient;

};

}  // namespace Opm
#endif