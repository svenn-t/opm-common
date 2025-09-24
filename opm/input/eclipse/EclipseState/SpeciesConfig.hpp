/*
  Copyright (C) 2020 Equinor

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

#ifndef SPECIES_CONFIG_HPP
#define SPECIES_CONFIG_HPP

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SpeciesVdTable.hpp>

#include <optional>
#include <string>

namespace Opm {

class SpeciesConfig {
public:
    struct SpeciesEntry {
        std::string name;
        std::optional<std::vector<double>> concentration;
        std::optional<SpeciesVdTable> svdp;

        SpeciesEntry() = default;

        SpeciesEntry(const std::string& name_, std::vector<double> concentration_)
            : name(name_)
            , concentration(concentration_)
        {}

        SpeciesEntry(const std::string& name_, SpeciesVdTable svdp_)
            : name(name_)
            , svdp(svdp_)
        {}

        SpeciesEntry(const std::string& name_)
            : name(name_)
        {}

        bool operator==(const SpeciesEntry& data) const {
            return this->name == data.name &&
                   this->concentration == data.concentration &&
                   this->svdp == data.svdp;
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(name);
            serializer(concentration);
            serializer(svdp);
        }
    }; // struct SpeciesEntry

    SpeciesConfig() = default;
    SpeciesConfig(const Deck& deck);

    static SpeciesConfig serializationTestObject();

    size_t size() const;
    bool empty() const;
    const std::vector<SpeciesEntry>::const_iterator begin() const;
    const std::vector<SpeciesEntry>::const_iterator end() const;
    const SpeciesEntry& operator[](const std::string& name) const;
    const SpeciesEntry& operator[](std::size_t index) const;
    bool operator==(const SpeciesConfig& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(species);
    }

private:
    std::vector<SpeciesEntry> species;
}; // class SpeciesConfig

} // namespec Opm
#endif