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

#include <opm/input/eclipse/EclipseState/Geochemistry/MineralConfig.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>

namespace Opm {

MineralConfig::MineralConfig(const Deck& deck) {
    using MINERAL = ParserKeywords::MINERAL;
    if (deck.hasKeyword<MINERAL>()) {
        const auto& keyword = deck.get<MINERAL>().back();
        const auto& item = keyword.getRecord(0).getItem<MINERAL::data>();

        OpmLog::info( keyword.location().format("\nInitializing species from {keyword} in {file} line {line}") );
        InfoLogger logger("Species tables", 3);
        for (std::size_t i = 0; i < item.getData<std::string>().size(); ++i) {
            const auto& species_name = item.getTrimmedString(i);
            std::string mblk_name = "MBLK" + species_name;
            std::string mvdp_name = "MVDP" + species_name;

            if (deck.hasKeyword(mblk_name)) {
                const auto& mblk_keyword = deck[mblk_name].back();
                this->initFromXBLK(mblk_keyword, species_name, logger);
            }
            else if (deck.hasKeyword(mvdp_name)) {
                const auto& mvdp_keyword = deck[mvdp_name].back();
                this->initFromXVDP(mvdp_keyword, species_name, logger);
            }
            else {
                this->initEmpty(species_name);
            }
        }
    }
}

}  // namespace Opm