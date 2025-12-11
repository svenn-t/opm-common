/*
  Copyright (C) 2025 Equinor

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
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/InfoLogger.hpp>

#include <opm/input/eclipse/EclipseState/Geochemistry/SpeciesConfig.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>

namespace Opm {

SpeciesConfig::SpeciesConfig(const Deck& deck) {
    using SPECIES = ParserKeywords::SPECIES;
    if (deck.hasKeyword<SPECIES>()) {
        const auto& keyword = deck.get<SPECIES>().back();
        const auto& item = keyword.getRecord(0).getItem<SPECIES::data>();

        OpmLog::info( keyword.location().format("\nInitializing species from {keyword} in {file} line {line}") );
        InfoLogger logger("Species tables", 3);
        for (std::size_t i = 0; i < item.getData<std::string>().size(); ++i) {
            const auto& species_name = item.getTrimmedString(i);
            std::string sblk_name = "SBLK" + species_name;
            std::string svdp_name = "SVDP" + species_name;

            if (deck.hasKeyword(sblk_name)) {
                const auto& sblk_keyword = deck[sblk_name].back();
                this->initFromXBLK(sblk_keyword, species_name, logger);
            }
            else if (deck.hasKeyword(svdp_name)) {
                const auto& svdp_keyword = deck[svdp_name].back();
                this->initFromXVDP(svdp_keyword, species_name, logger);
            }
            else {
                this->initEmpty(species_name);
            }
        }
    }
}

}  // namespace Opm