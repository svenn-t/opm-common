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

#ifndef OPM_OUTPUT_CELLS_HPP
#define OPM_OUTPUT_CELLS_HPP

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <utility>
#include <variant>
#include <vector>

namespace Opm { namespace data {

    // The 3D data which are saved to file are assembled in one large
    // container.  In the container the data is tagged with an element from
    // the TargetType enum which specifies the vector's intended output
    // destination.
    //
    //RESTART_SOLUTION: Cell-based quantities that are output to the
    //  SOLUTION section of the restart file.  ECLIPSE-compatible names.
    //
    //RESTART_AUXILIARY: Fields with extra information, not required
    //  for restart.  Examples of this include fluid in place values or
    //  evaluations of relative permeability. Will end up in the
    //  restart file.  Deprecated and will be removed.
    //
    //SUMMARY: Fields which are added only to serve as input data for
    //  calculations of summary results. The Summary implementation can
    //  use data with any tag value, but if it is tagged as SUMMARY it
    //  will not be output anywhere else.
    //
    //INIT: Fields which should go to the INIT file.
    //
    //RESTART_OPM_EXTENDED: Cell-based quantities that are specific to
    //  OPM-Flow.  Output only to extended OPM restart files.  Specifically
    //  not output to ECLIPSE-compatible restart files.

    enum class TargetType
    {
        RESTART_SOLUTION,
        RESTART_AUXILIARY,
        RESTART_TRACER_SOLUTION,
        SUMMARY,
        INIT,
        RESTART_OPM_EXTENDED,
    };

    /// Small struct that keeps track of data for output to restart/summary
    /// files.
    struct CellData
    {
        /// Dimension of the data to write
        UnitSystem::measure dim{UnitSystem::measure::identity};

        /// File output destination
        TargetType target{TargetType::RESTART_SOLUTION};

        CellData() = default;
        explicit CellData(UnitSystem::measure m,
                          std::vector<double> x,
                          TargetType          dest)
            : dim    { m }
            , target { dest }
            , data_  { std::move(x) }
        {}

        explicit CellData(std::vector<int> x,
                          TargetType          dest)
            : dim    { UnitSystem::measure::identity }
            , target { dest }
            , data_  { std::move(x) }
        {}

        bool operator==(const CellData& cell2) const
        {
            return (dim    == cell2.dim)
                && (target == cell2.target)
                && (data_   == cell2.data_);
        }

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->dim);
            serializer(this->data_);
            serializer(this->target);
        }

        static CellData serializationTestObject()
        {
            return CellData {
                UnitSystem::measure::runtime,
                std::vector<double>{1.0, 2.0, 3.0},
                TargetType::RESTART_OPM_EXTENDED
            };
        }

        template<class T>
        std::vector<T>& data()
        {
            return std::get<std::vector<T>>(data_);
        }

        template<class T>
        const std::vector<T>& data() const
        {
            return std::get<std::vector<T>>(data_);
        }

        template<class Visitor>
        void visit(Visitor&& visit) const
        {
            std::visit(std::forward<Visitor>(visit), data_);
        }

    private:
        /// Per-cell solution values
        using DataVector = std::variant<std::monostate,
                                        std::vector<double>,
                                        std::vector<int>>;

        DataVector data_{};

    };

}} // namespace Opm::data

#endif //OPM_OUTPUT_CELLS_HPP
