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

#include <opm/input/eclipse/EclipseState/Geochemistry/GenericSpeciesConfig.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>

#include <fmt/format.h>

#include <algorithm>

namespace Opm {

void GenericSpeciesConfig::initFromXBLK(const DeckKeyword& sblk_keyword,
                            const std::string& species_name,
                            InfoLogger& logger)
{
    double inv_volume = 1.0;
    auto sblk_conc = sblk_keyword.getRecord(0).getItem(0).getData<double>();
    logger(sblk_keyword.location().format("Loading species concentration from {keyword} in {file} line {line}"));

    std::transform(sblk_conc.begin(), sblk_conc.end(), sblk_conc.begin(),
                    [inv_volume](const auto& c) { return c * inv_volume; });

    this->species.emplace_back(species_name, std::move(sblk_conc));
}

void GenericSpeciesConfig::initFromXVDP(const DeckKeyword& svdp_keyword,
                      const std::string& species_name,
                      InfoLogger& logger)
{
    double inv_volume = 1.0;
    auto svdp_table = svdp_keyword.getRecord(0).getItem(0);
    logger(svdp_keyword.location().format("Loading species concentration from {keyword} in {file} line {line}"));

    this->species.emplace_back(species_name, SpeciesVdTable(svdp_table, inv_volume, species.size()));
}

void GenericSpeciesConfig::initEmpty(const std::string& species_name)
{
    this->species.emplace_back(species_name);
}

GenericSpeciesConfig GenericSpeciesConfig::serializationTestObject()
{
    GenericSpeciesConfig result;
    result.species = {{"test", {1.0}}};

    return result;
}

const GenericSpeciesConfig::SpeciesEntry& GenericSpeciesConfig::operator[](std::size_t index) const {
    return this->species.at(index);
}

const GenericSpeciesConfig::SpeciesEntry& GenericSpeciesConfig::operator[](const std::string& name) const {
    auto iter = std::find_if(this->species.begin(), this->species.end(),
                                [&name](const SpeciesEntry& single_species)
                                { return single_species.name == name;}
                            );

    if (iter == this->species.end())
        throw std::logic_error(fmt::format("No such species: {}", name));

    return *iter;
}

std::size_t GenericSpeciesConfig::size() const {
    return this->species.size();
}

bool GenericSpeciesConfig::empty() const {
    return this->species.empty();
}

const std::vector<GenericSpeciesConfig::SpeciesEntry>::const_iterator GenericSpeciesConfig::begin() const {
    return this->species.begin();
}

const std::vector<GenericSpeciesConfig::SpeciesEntry>::const_iterator GenericSpeciesConfig::end() const {
    return this->species.end();
}

bool GenericSpeciesConfig::operator==(const GenericSpeciesConfig& other) const {
    return this->species == other.species;
}

} // namespace Opm