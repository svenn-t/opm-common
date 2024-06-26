/*
  Copyright 2016 Statoil ASA.

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

#ifndef OPM_REGION_CACHE_HPP
#define OPM_REGION_CACHE_HPP

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Opm {
    class Schedule;
    class EclipseGrid;
    class FieldPropsManager;
} // namespace Opm

namespace Opm { namespace out {
    class RegionCache {
    public:
        RegionCache() = default;
        RegionCache(const std::set<std::string>& fip_regions,
                    const FieldPropsManager&     fp,
                    const EclipseGrid&           grid,
                    const Schedule&              schedule);

        void buildCache(const std::set<std::string>& fip_regions,
                        const FieldPropsManager&     fp,
                        const EclipseGrid&           grid,
                        const Schedule&              schedule);

        const std::vector<std::pair<std::string, std::size_t>>&
        connections(const std::string& region_name, int region_id) const;

        // A well is assigned to the region_id of its first connection.
        std::vector<std::string> wells(const std::string& region_name, int region_id) const;

    private:
        using RegID = std::pair<std::string, int>;            // { Region set, region ID }
        using WellConn = std::pair<std::string, std::size_t>; // { Well name, cell ID }

        std::vector<WellConn> connections_empty{};
        std::map<RegID, std::vector<WellConn>> connection_map{};
        std::map<RegID, std::vector<std::string>> well_map{};
    };
}} // namespace Opm::out

#endif // OPM_REGION_CACHE_HPP
