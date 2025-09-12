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

#include <config.h>

#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/InfoLogger.hpp>

#include <opm/input/eclipse/EclipseState/SpeciesConfig.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>

#include <fmt/format.h>

#include <algorithm>

namespace Opm {

SpeciesConfig::SpeciesConfig(const Deck& deck)
{
    using SPECIES = ParserKeywords::SPECIES;
    if (deck.hasKeyword<SPECIES>()) {
        const auto& keyword = deck.get<SPECIES>().back();
        const auto& item = keyword.getRecord(0).getItem<SPECIES::data>();
        double inv_volume = 1.0;

        OpmLog::info( keyword.location().format("\nInitializing species from {keyword} in {file} line {line}") );
        InfoLogger logger("Species tables", 3);
        for (std::size_t i = 0; i < item.getData<std::string>().size(); ++i) {
            const auto& species_name = item.getTrimmedString(i);
            std::string sblk_name = "SBLK" + species_name;
            std::string svdp_name = "SVDP" + species_name;

            if (deck.hasKeyword(sblk_name)) {
                const auto& sblk_keyword = deck[sblk_name].back();
                auto sblk_conc = sblk_keyword.getRecord(0).getItem(0).getData<double>();
                logger(sblk_keyword.location().format("Loading species concentration from {keyword} in {file} line {line}"));

                std::transform(sblk_conc.begin(), sblk_conc.end(), sblk_conc.begin(),
                               [inv_volume](const auto& c) { return c * inv_volume; });
                
                this->species.emplace_back(species_name, std::move(sblk_conc));
            }
            else if (deck.hasKeyword(svdp_name)) {
                const auto& svdp_keyword = deck[svdp_name].back();
                auto svdp_table = svdp_keyword.getRecord(0).getItem(0);
                logger(svdp_keyword.location().format("Loading species concentration from {keyword} in {file} line {line}"));

                this->species.emplace_back(species_name, SpeciesVdTable(svdp_table, inv_volume, species.size()));
            }
            else {
                this->species.emplace_back(species_name);
            }
        }
    }
}

SpeciesConfig SpeciesConfig::serializationTestObject()
{
    SpeciesConfig result;
    result.species = {{"test", {1.0}}};

    return result;
}

const SpeciesConfig::SpeciesEntry& SpeciesConfig::operator[](std::size_t index) const {
    return this->species.at(index);
}

const SpeciesConfig::SpeciesEntry& SpeciesConfig::operator[](const std::string& name) const {
    auto iter = std::find_if(this->species.begin(), this->species.end(), 
                                [&name](const SpeciesEntry& single_species) 
                                { return single_species.name == name;}
                            );

    if (iter == this->species.end())
        throw std::logic_error(fmt::format("No such species: {}", name));

    return *iter;
}

std::size_t SpeciesConfig::size() const {
    return this->species.size();
}

bool SpeciesConfig::empty() const {
    return this->species.empty();
}

const std::vector<SpeciesConfig::SpeciesEntry>::const_iterator SpeciesConfig::begin() const {
    return this->species.begin();
}

const std::vector<SpeciesConfig::SpeciesEntry>::const_iterator SpeciesConfig::end() const {
    return this->species.end();
}

bool SpeciesConfig::operator==(const SpeciesConfig& other) const {
    return this->species == other.species;
}

} // namespace Opm