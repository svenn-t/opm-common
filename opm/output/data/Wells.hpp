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

#ifndef OPM_OUTPUT_WELLS_HPP
#define OPM_OUTPUT_WELLS_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <opm/json/JsonObject.hpp>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace Opm { namespace data {

    class Rates {
        /* Methods are defined inline for performance, as the actual *work* done
         * is trivial, but somewhat frequent (typically once per time step per
         * completion per well).
         *
         * To add a new rate type, add an entry in the enum with the correct
         * shift, and if needed, increase the size type. Add a member variable
         * and a new case in get_ref.
         */

        public:
            Rates() = default;
            enum class opt : uint32_t {
                wat               = (1 << 0),
                oil               = (1 << 1),
                gas               = (1 << 2),
                polymer           = (1 << 3),
                solvent           = (1 << 4),
                energy            = (1 << 5),
                dissolved_gas     = (1 << 6),
                vaporized_oil     = (1 << 7),
                reservoir_water   = (1 << 8),
                reservoir_oil     = (1 << 9),
                reservoir_gas     = (1 << 10),
                productivity_index_water = (1 << 11),
                productivity_index_oil   = (1 << 12),
                productivity_index_gas   = (1 << 13),
                well_potential_water   = (1 << 14),
                well_potential_oil     = (1 << 15),
                well_potential_gas     = (1 << 16),
                brine            = (1 << 17),
                alq              = (1 << 18),
                tracer           = (1 << 19),
                microbial        = (1 << 20),
                oxygen           = (1 << 21),
                urea             = (1 << 22),
                vaporized_water  = (1 << 23),
                mass_gas         = (1 << 24),
                mass_wat         = (1 << 25)
            };

            using enum_size = std::underlying_type< opt >::type;

            /// Query if a value is set.
            inline bool has( opt ) const;

            /// Read the value indicated by m. Throws an exception if
            /// if the requested value is unset.
            inline double get( opt m ) const;
            /// Read the value indicated by m. Returns a default value if
            /// the requested value is unset.
            inline double get( opt m, double default_value ) const;
            inline double get( opt m, double default_value , const std::string& tracer_name ) const;
            /// Set the value specified by m. Throws an exception if multiple
            /// values are requested. Returns a self-reference to support
            /// chaining.
            inline Rates& set( opt m, double value );
            inline Rates& set( opt m, double value , const std::string& tracer_name );

            /// Returns true if any of the rates oil, gas, water is nonzero
            inline bool flowing() const;

            template <class MessageBufferType>
            void write(MessageBufferType& buffer) const;
            template <class MessageBufferType>
            void read(MessageBufferType& buffer);

            bool operator==(const Rates& rat2) const;

            inline void init_json(Json::JsonObject& json_data) const;

            template<class Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(mask);
                serializer(wat);
                serializer(oil);
                serializer(gas);
                serializer(polymer);
                serializer(solvent);
                serializer(energy);
                serializer(dissolved_gas);
                serializer(vaporized_oil);
                serializer(reservoir_water);
                serializer(reservoir_oil);
                serializer(reservoir_gas);
                serializer(productivity_index_water);
                serializer(productivity_index_oil);
                serializer(productivity_index_gas);
                serializer(well_potential_water);
                serializer(well_potential_oil);
                serializer(well_potential_gas);
                serializer(brine);
                serializer(alq);
                serializer(tracer);
                serializer(microbial);
                serializer(oxygen);
                serializer(urea);
                serializer(vaporized_water);
                serializer(mass_gas);
                serializer(mass_wat);
            }

            static Rates serializationTestObject()
            {
                Rates rat1;
                rat1.set(opt::wat, 1.0);
                rat1.set(opt::oil, 2.0);
                rat1.set(opt::gas, 3.0);
                rat1.set(opt::polymer, 4.0);
                rat1.set(opt::solvent, 5.0);
                rat1.set(opt::energy, 6.0);
                rat1.set(opt::dissolved_gas, 7.0);
                rat1.set(opt::vaporized_oil, 8.0);
                rat1.set(opt::reservoir_water, 9.0);
                rat1.set(opt::reservoir_oil, 10.0);
                rat1.set(opt::reservoir_gas, 11.0);
                rat1.set(opt::productivity_index_water, 12.0);
                rat1.set(opt::productivity_index_oil, 13.0);
                rat1.set(opt::productivity_index_gas, 14.0);
                rat1.set(opt::well_potential_water, 15.0);
                rat1.set(opt::well_potential_oil, 16.0);
                rat1.set(opt::well_potential_gas, 17.0);
                rat1.set(opt::brine, 18.0);
                rat1.set(opt::alq, 19.0);
                rat1.set(opt::microbial, 20.0);
                rat1.set(opt::oxygen, 21.0);
                rat1.set(opt::urea, 22.0);
                rat1.set(opt::vaporized_water, 23.0);
                rat1.set(opt::mass_gas, 24.0);
                rat1.set(opt::mass_wat, 25.0);
                rat1.tracer.insert({"test_tracer", 1.0});

                return rat1;
            }

        private:
            double& get_ref( opt );
            double& get_ref( opt, const std::string& tracer_name );
            const double& get_ref( opt ) const;
            const double& get_ref( opt, const std::string& tracer_name ) const;

            opt mask = static_cast< opt >( 0 );

            double wat = 0.0;
            double oil = 0.0;
            double gas = 0.0;
            double polymer = 0.0;
            double solvent = 0.0;
            double energy = 0.0;
            double dissolved_gas = 0.0;
            double vaporized_oil = 0.0;
            double reservoir_water = 0.0;
            double reservoir_oil = 0.0;
            double reservoir_gas = 0.0;
            double productivity_index_water = 0.0;
            double productivity_index_oil = 0.0;
            double productivity_index_gas = 0.0;
            double well_potential_water = 0.0;
            double well_potential_oil = 0.0;
            double well_potential_gas = 0.0;
            double brine = 0.0;
            double alq = 0.0;
            std::map<std::string, double> tracer{};
            double microbial = 0.0;
            double oxygen = 0.0;
            double urea = 0.0;
            double vaporized_water = 0.0;
            double mass_gas = 0.0;
            double mass_wat = 0.0;
    };

    struct ConnectionFiltrate
    {
        double rate;
        double total;
        double skin_factor;
        double thickness;
        double perm;
        double poro;
        double radius;
        double area_of_flow;

        template<class Serializer>
        void serializeOp(Serializer& serializer) {
            serializer(rate);
            serializer(total);
            serializer(skin_factor);
            serializer(thickness);
            serializer(perm);
            serializer(poro);
            serializer(radius);
            serializer(area_of_flow);
        }

        bool operator==(const ConnectionFiltrate& filtrate) const
        {
            return this->rate == filtrate.rate &&
                   this->total == filtrate.total &&
                   this->skin_factor == filtrate.skin_factor &&
                   this->thickness == filtrate.thickness &&
                   this->perm == filtrate.perm &&
                   this->poro == filtrate.poro &&
                   this->radius == filtrate.radius &&
                   this->area_of_flow == filtrate.area_of_flow;
        }

        static ConnectionFiltrate serializationTestObject()
        {
            return {0.8, 100., -1., 2., 1.e-9,
                    0.3, 0.05, 0.8};
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);
    };

    /// Connection Level Fracturing Statistics
    struct ConnectionFracturing
    {
        /// Statistics collection for a single quantity
        struct Statistics
        {
            /// Arithmetic average.
            double avg{};

            /// Maximum value.
            double max{};

            /// Minimum value.
            double min{};

            /// Unbiased sample standard deviation.
            ///
            /// Usable only if sample size is at least two.
            double stdev{};

            /// Create a serialization test object.
            static Statistics serializationTestObject()
            {
                return {
                    12.34, 56.78, 9.10, 11.12
                };
            }

            /// Convert between byte array and object representation.
            ///
            /// \tparam Serializer Byte array conversion protocol.
            ///
            /// \param[in,out] serializer Byte array conversion object.
            template <class Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(this->avg);
                serializer(this->max);
                serializer(this->min);
                serializer(this->stdev);
            }

            /// Equality predicate.
            ///
            /// \param[in] that Object against which \code *this \endcode
            /// will be tested for equality.
            ///
            /// \return Whether or not \code *this \endcode is the same as
            /// \p that.
            bool operator==(const Statistics& that) const
            {
                return (this->avg == that.avg)
                    && (this->max == that.max)
                    && (this->min == that.min)
                    && (this->stdev == that.stdev)
                    ;
            }

            /// MPI communication protocol--serialisation operation
            template <class MessageBufferType>
            void write(MessageBufferType& buffer) const
            {
                buffer.write(this->avg);
                buffer.write(this->max);
                buffer.write(this->min);
                buffer.write(this->stdev);
            }

            /// MPI communication protocol--deserialisation operation
            template <class MessageBufferType>
            void read(MessageBufferType& buffer)
            {
                buffer.read(this->avg);
                buffer.read(this->max);
                buffer.read(this->min);
                buffer.read(this->stdev);
            }
        };

        /// Sample size.
        ///
        /// Expected to be the same for each quantiy.
        std::size_t numCells{};

        /// Statistical measures for connection's fracture pressures.
        Statistics press{};

        /// Statistical measures for connection's fracture fracture flow rate.
        Statistics rate{};

        /// Statistical measures for connection's fracture fracture width.
        Statistics width{};

        /// Create a serialisation test object.
        static ConnectionFracturing serializationTestObject()
        {
            auto fract = ConnectionFracturing{};

            fract.numCells = 123;
            fract.press = Statistics::serializationTestObject();
            fract.rate = Statistics::serializationTestObject();
            fract.width = Statistics::serializationTestObject();

            return fract;
        }

        /// Convert between byte array and object representation.
        ///
        /// \tparam Serializer Byte array conversion protocol.
        ///
        /// \param[in,out] serializer Byte array conversion object.
        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->numCells);
            serializer(this->press);
            serializer(this->rate);
            serializer(this->width);
        }

        /// Equality predicate.
        ///
        /// \param[in] that Object against which \code *this \endcode will
        /// be tested for equality.
        ///
        /// \return Whether or not \code *this \endcode is the same as \p
        /// that.
        bool operator==(const ConnectionFracturing& that) const
        {
            return (this->numCells == that.numCells)
                && (this->press == that.press)
                && (this->rate == that.rate)
                && (this->width == that.width)
                ;
        }

        /// MPI communication protocol--serialisation operation
        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const
        {
            buffer.write(this->numCells);
            buffer.write(this->press);
            buffer.write(this->rate);
            buffer.write(this->width);
        }

        /// MPI communication protocol--deserialisation operation
        template <class MessageBufferType>
        void read(MessageBufferType& buffer)
        {
            buffer.read(this->numCells);
            buffer.read(this->press);
            buffer.read(this->rate);
            buffer.read(this->width);
        }
    };

    struct Connection
    {
        using global_index = std::size_t;
        static const constexpr int restart_size = 6;

        global_index index{};
        Rates rates{};
        double pressure{};
        double reservoir_rate{};
        double cell_pressure{};
        double cell_saturation_water{};
        double cell_saturation_gas{};
        double effective_Kh{};
        double trans_factor{};
        double d_factor{};
        double compact_mult{1.0}; // Rock compaction transmissibility multiplier (ROCKTAB)

        ConnectionFiltrate filtrate{};

        /// Connection level fracturing statistics.
        ConnectionFracturing fract{};

        bool operator==(const Connection& conn2) const
        {
            return (index == conn2.index)
                && (rates == conn2.rates)
                && (pressure == conn2.pressure)
                && (reservoir_rate == conn2.reservoir_rate)
                && (cell_pressure == conn2.cell_pressure)
                && (cell_saturation_water == conn2.cell_saturation_water)
                && (cell_saturation_gas == conn2.cell_saturation_gas)
                && (effective_Kh == conn2.effective_Kh)
                && (trans_factor == conn2.trans_factor)
                && (d_factor == conn2.d_factor)
                && (compact_mult == conn2.compact_mult)
                && (filtrate == conn2.filtrate)
                && (this->fract == conn2.fract)
                ;
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;
        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        inline void init_json(Json::JsonObject& json_data) const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(index);
            serializer(rates);
            serializer(pressure);
            serializer(reservoir_rate);
            serializer(cell_pressure);
            serializer(cell_saturation_water);
            serializer(cell_saturation_gas);
            serializer(effective_Kh);
            serializer(trans_factor);
            serializer(d_factor);
            serializer(compact_mult);
            serializer(filtrate);
            serializer(this->fract);
        }

        static Connection serializationTestObject()
        {
            return Connection {
                1, Rates::serializationTestObject(),
                2.0, 3.0, 4.0, 5.0,
                6.0, 7.0, 8.0, 9.0, 0.987,
                ConnectionFiltrate::serializationTestObject(),
                ConnectionFracturing::serializationTestObject()
            };
        }
    };

    class SegmentPressures
    {
    public:
        enum class Value : std::size_t {
            Pressure, PDrop, PDropHydrostatic, PDropAccel, PDropFriction,
        };

        double& operator[](const Value i)
        {
            return this->values_[this->index(i)];
        }

        double operator[](const Value i) const
        {
            return this->values_[this->index(i)];
        }

        bool operator==(const SegmentPressures& segpres2) const
        {
            return this->values_ == segpres2.values_;
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const
        {
            for (const auto& value : this->values_) {
                buffer.write(value);
            }
        }

        template <class MessageBufferType>
        void read(MessageBufferType& buffer)
        {
            for (auto& value : this->values_) {
                buffer.read(value);
            }
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(values_);
        }

        static SegmentPressures serializationTestObject()
        {
            SegmentPressures spres;
            spres[Value::Pressure] = 1.0;
            spres[Value::PDrop] = 2.0;
            spres[Value::PDropHydrostatic] = 3.0;
            spres[Value::PDropAccel] = 4.0;
            spres[Value::PDropFriction] = 5.0;

            return spres;
        }

    private:
        constexpr static std::size_t numvals = 5;

        std::array<double, numvals> values_ = {0};

        std::size_t index(const Value ix) const
        {
            return static_cast<std::size_t>(ix);
        }
    };

    template <typename Items>
    class QuantityCollection
    {
    public:
        using Item = typename Items::Item;

        void clear()
        {
            this->has_ = static_cast<unsigned char>(0);
            this->value_.fill(0.0);
        }

        constexpr bool has(const Item p) const
        {
            const auto i = this->index(p);

            return (i < Size) && this->hasItem(i);
        }

        bool operator==(const QuantityCollection& that) const
        {
            return (this->has_   == that.has_)
                && (this->value_ == that.value_);
        }

        double get(const Item p) const
        {
            if (! this->has(p)) {
                throw std::invalid_argument {
                    "Request for Unset Item Value for " + Items::itemName(p)
                };
            }

            return this->value_[ this->index(p) ];
        }

        QuantityCollection& set(const Item p, const double value)
        {
            const auto i = this->index(p);

            if (i >= Size) {
                throw std::invalid_argument {
                    "Cannot Assign Item Value for Unsupported Item '"
                    + Items::itemName(p) + '\''
                };
            }

            this->has_ |= 1 << i;
            this->value_[i] = value;

            return *this;
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const
        {
            buffer.write(this->has_);

            for (const auto& x : this->value_) {
                buffer.write(x);
            }
        }

        template <class MessageBufferType>
        void read(MessageBufferType& buffer)
        {
            this->clear();
            buffer.read(this->has_);

            for (auto& x : this->value_) {
                buffer.read(x);
            }
        }

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->has_);
            serializer(this->value_);
        }

        static QuantityCollection serializationTestObject()
        {
            auto quant = QuantityCollection{};

            for (const auto& [item, value] : Items::serializationTestItems()) {
                quant.set(item, value);
            }

            return quant;
        }

    private:
        enum { Size = static_cast<std::size_t>(Item::NumItems) };

        static_assert(Size <= static_cast<std::size_t>(CHAR_BIT),
                      "Number of items must not exceed CHAR_BIT");

        /// Whether or not item has a defined value.  We use the bottom
        /// 'Size' bits.
        unsigned char has_{};

        /// Numerical value of each item.
        std::array<double, Size> value_{};

        constexpr std::size_t index(const Item p) const noexcept
        {
            return static_cast<std::size_t>(p);
        }

        bool hasItem(const std::size_t i) const
        {
            return (this->has_ & (1 << i)) != 0;
        }
    };

    struct PhaseItems
    {
        enum class Item {
            Oil, Gas, Water,

            // -- Must be last enumerator --
            NumItems,
        };

        static std::string itemName(const Item p)
        {
            switch (p) {
            case Item::Oil:   return "Oil";
            case Item::Gas:   return "Gas";
            case Item::Water: return "Water";

            case Item::NumItems:
                return "Out of bounds (NumItems)";
            }

            return "Unknown (" + std::to_string(static_cast<int>(p)) + ')';
        }

        static auto serializationTestItems()
        {
            return std::vector {
                std::pair { Item::Oil  , 1.0 },
                std::pair { Item::Gas  , 7.0 },
                std::pair { Item::Water, 2.9 },
            };
        }
    };

    struct DensityItems
    {
        enum class Item {
            Oil, Gas, Water, Mixture, MixtureWithExponents,

            // -- Must be last enumerator --
            NumItems,
        };

        static std::string itemName(const Item p)
        {
            switch (p) {
            case Item::Oil:                  return "Oil";
            case Item::Gas:                  return "Gas";
            case Item::Water:                return "Water";
            case Item::Mixture:              return "Mixture";
            case Item::MixtureWithExponents: return "MixtureWithExponents";

            case Item::NumItems:
                return "Out of bounds (NumItems)";
            }

            return "Unknown (" + std::to_string(static_cast<int>(p)) + ')';
        }

        static auto serializationTestItems()
        {
            return std::vector {
                std::pair { Item::Oil                 , 876.54 },
                std::pair { Item::Gas                 , 321.09 },
                std::pair { Item::Water               , 987.65 },
                std::pair { Item::Mixture             , 975.31 },
                std::pair { Item::MixtureWithExponents, 765.43 },
            };
        }
    };

    using SegmentPhaseQuantity = QuantityCollection<PhaseItems>;
    using SegmentPhaseDensity = QuantityCollection<DensityItems>;

    struct Segment
    {
        Rates rates{};
        SegmentPressures pressures{};
        SegmentPhaseQuantity velocity{};
        SegmentPhaseQuantity holdup{};
        SegmentPhaseQuantity viscosity{};
        SegmentPhaseDensity density{};
        std::size_t segNumber{};

        bool operator==(const Segment& seg2) const
        {
            return (rates == seg2.rates)
                && (pressures == seg2.pressures)
                && (velocity == seg2.velocity)
                && (holdup == seg2.holdup)
                && (viscosity == seg2.viscosity)
                && (density == seg2.density)
                && (segNumber == seg2.segNumber);
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->rates);
            serializer(this->pressures);
            serializer(this->velocity);
            serializer(this->holdup);
            serializer(this->viscosity);
            serializer(this->density);
            serializer(this->segNumber);
        }

        static Segment serializationTestObject()
        {
            return {
                Rates::serializationTestObject(),
                SegmentPressures::serializationTestObject(),
                SegmentPhaseQuantity::serializationTestObject(), // velocity
                SegmentPhaseQuantity::serializationTestObject(), // holdup
                SegmentPhaseQuantity::serializationTestObject(), // viscosity
                SegmentPhaseDensity::serializationTestObject(),  // density
                10
            };
        }
    };

    struct CurrentControl
    {
        bool isProducer{true};

        ::Opm::WellProducerCMode prod {
            ::Opm::WellProducerCMode::CMODE_UNDEFINED
        };

        ::Opm::WellInjectorCMode inj {
            ::Opm::WellInjectorCMode::CMODE_UNDEFINED
        };

        bool operator==(const CurrentControl& rhs) const
        {
            return (this->isProducer == rhs.isProducer)
                && ((this->isProducer && (this->prod == rhs.prod)) ||
                    (!this->isProducer && (this->inj == rhs.inj)));
        }

        void init_json(Json::JsonObject& json_data) const
        {
            if (this->inj == ::Opm::WellInjectorCMode::CMODE_UNDEFINED)
                json_data.add_item("inj", "CMODE_UNDEFINED");
            else
                json_data.add_item("inj", ::Opm::WellInjectorCMode2String(this->inj));

            if (this->prod == ::Opm::WellProducerCMode::CMODE_UNDEFINED)
                json_data.add_item("prod", "CMODE_UNDEFINED");
            else
                json_data.add_item("prod", ::Opm::WellProducerCMode2String(this->prod));
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(isProducer);
            serializer(prod);
            serializer(inj);
        }

        static CurrentControl serializationTestObject()
        {
          return CurrentControl{false,
                                ::Opm::WellProducerCMode::BHP,
                                ::Opm::WellInjectorCMode::GRUP
                 };
        }
    };

    class WellBlockAvgPress
    {
    public:
        enum class Quantity { WBP, WBP4, WBP5, WBP9 };

        double& operator[](const Quantity q)
        {
            return this->wbp_[static_cast<std::size_t>(q)];
        }

        double operator[](const Quantity q) const
        {
            return this->wbp_[static_cast<std::size_t>(q)];
        }

        bool operator==(const WellBlockAvgPress& that) const
        {
            return this->wbp_ == that.wbp_;
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->wbp_);
        }

        static WellBlockAvgPress serializationTestObject()
        {
            auto wbp = WellBlockAvgPress{};

            wbp[Quantity::WBP]  = 17.29;
            wbp[Quantity::WBP4] =  2.718;
            wbp[Quantity::WBP5] =  3.1415;
            wbp[Quantity::WBP9] =  1.618;

            return wbp;
        }

    private:
        static constexpr auto NumQuantities =
            static_cast<std::size_t>(Quantity::WBP9) + 1;

        std::array<double, NumQuantities> wbp_{};
    };

    struct WellFiltrate
    {
        double rate{0.};
        double total{0.};
        double concentration{0.};

        template<class Serializer>
        void serializeOp(Serializer& serializer) {
            serializer(rate);
            serializer(total);
            serializer(concentration);
        }

        bool operator==(const WellFiltrate& filtrate) const {
           return this->rate == filtrate.rate
              && this->total == filtrate.total
              && this->concentration == filtrate.concentration;
        }

        static WellFiltrate serializationTestObject() {
            WellFiltrate res;
            res.rate = 1.;
            res.total = 10.;
            res.concentration = 0.;
            return res;
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);
    };

    struct WellControlLimitItems
    {
        enum class Item {
            Bhp, OilRate, WaterRate, GasRate, ResVRate, LiquidRate,

            // -- Must be last enumerator --
            NumItems,
        };

        static std::string itemName(const Item p)
        {
            switch (p) {
            case Item::Bhp:        return "Bhp";
            case Item::OilRate:    return "OilRate";
            case Item::WaterRate:  return "WaterRate";
            case Item::GasRate:    return "GasRate";
            case Item::ResVRate:   return "ResVRate";
            case Item::LiquidRate: return "LiquidRate";

            case Item::NumItems:
                return "Out of bounds (NumItems)";
            }

            return "Unknown (" + std::to_string(static_cast<int>(p)) + ')';
        }

        static auto serializationTestItems()
        {
            return std::vector {
                std::pair { Item::Bhp       , 321.09 },
                std::pair { Item::OilRate   , 987.65 },
                std::pair { Item::WaterRate , 975.31 },
                std::pair { Item::GasRate   , 765.43 },
                std::pair { Item::ResVRate  , 876.54 },
                std::pair { Item::LiquidRate,  54.32 },
            };
        }
    };

    using WellControlLimits = QuantityCollection<WellControlLimitItems>;

    struct Well
    {
        Rates rates{};

        double bhp{0.0};
        double thp{0.0};
        double temperature{0.0};
        int control{0};
        double efficiency_scaling_factor{1.0};

        WellFiltrate filtrate;

        ::Opm::WellStatus dynamicStatus { Opm::WellStatus::OPEN };

        std::vector<Connection> connections{};
        std::unordered_map<std::size_t, Segment> segments{};
        CurrentControl current_control{};
        GuideRateValue guide_rates{};
        WellControlLimits limits{};

        inline bool flowing() const noexcept;

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        inline void init_json(Json::JsonObject& json_data) const;

        const Connection*
        find_connection(const Connection::global_index connection_grid_index) const
        {
            auto connection = std::find_if(this->connections.begin(),
                                           this->connections.end(),
                                           [connection_grid_index](const Connection& c)
                                           { return c.index == connection_grid_index; });

            if (connection == this->connections.end()) {
                return nullptr;
            }

            return &*connection;
        }

        Connection*
        find_connection(const Connection::global_index connection_grid_index)
        {
            auto connection = std::find_if(this->connections.begin(),
                                           this->connections.end(),
                                           [connection_grid_index](const Connection& c)
                                           { return c.index == connection_grid_index; });

            if (connection == this->connections.end()) {
                return nullptr;
            }

            return &*connection;
        }

        bool operator==(const Well& well2) const
        {
            return (this->rates == well2.rates)
                && (this->bhp == well2.bhp)
                && (this->thp == well2.thp)
                && (this->temperature == well2.temperature)
                && (this->filtrate == well2.filtrate)
                && (this->control == well2.control)
                && (this->dynamicStatus == well2.dynamicStatus)
                && (this->connections == well2.connections)
                && (this->segments == well2.segments)
                && (this->current_control == well2.current_control)
                && (this->guide_rates == well2.guide_rates)
                && (this->limits == well2.limits)
                ;
        }

        bool operator!=(const Well& well2) const
        {
            return !(*this == well2);
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(rates);
            serializer(bhp);
            serializer(thp);
            serializer(temperature);
            serializer(control);
            serializer(efficiency_scaling_factor);
            serializer(filtrate);
            serializer(dynamicStatus);
            serializer(connections);
            serializer(segments);
            serializer(current_control);
            serializer(guide_rates);
            serializer(limits);
        }

        static Well serializationTestObject()
        {
            return Well {
                Rates::serializationTestObject(),
                1.0,
                2.0,
                3.0,
                4,
                5.0,
                WellFiltrate::serializationTestObject(),
                ::Opm::WellStatus::SHUT,
                {Connection::serializationTestObject()},
                {{0, Segment::serializationTestObject()}},
                CurrentControl::serializationTestObject(),
                GuideRateValue::serializationTestObject(),
                WellControlLimits::serializationTestObject()
            };
        }
    };

    class Wells: public std::map<std::string , Well> {
    public:

        double get(const std::string& well_name , Rates::opt m) const {
            const auto& well = this->find( well_name );
            if( well == this->end() ) return 0.0;

            return well->second.rates.get( m, 0.0 );
        }

        double get(const std::string& well_name , Rates::opt m, const std::string& tracer_name) const {
            const auto& well = this->find( well_name );
            if( well == this->end() ) return 0.0;

            return well->second.rates.get( m, 0.0, tracer_name);
        }

        double get(const std::string& well_name , Connection::global_index connection_grid_index, Rates::opt m) const {
            const auto& witr = this->find( well_name );
            if( witr == this->end() ) return 0.0;

            const auto& well = witr->second;
            const auto& connection = std::find_if( well.connections.begin() ,
                                                   well.connections.end() ,
                                                   [=]( const Connection& c ) {
                                                        return c.index == connection_grid_index; });

            if( connection == well.connections.end() )
                return 0.0;

            return connection->rates.get( m, 0.0 );
        }

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const {
            unsigned int size = this->size();
            buffer.write(size);
            for (const auto& witr : *this) {
                const std::string& name = witr.first;
                buffer.write(name);
                const Well& well = witr.second;
                well.write(buffer);
            }
        }

        template <class MessageBufferType>
        void read(MessageBufferType& buffer) {
            unsigned int size;
            buffer.read(size);
            for (size_t i = 0; i < size; ++i) {
                std::string name;
                buffer.read(name);
                Well well;
                well.read(buffer);
                auto result = this->emplace(name, well);
                // In case there was already an entry for the well we want to insert, then result.second == false.
                // Then we check if this entry is the same as the one we want to insert.
                if (!result.second && result.first->second != well) {
                    OPM_THROW(std::runtime_error, "Received different output data for well " + name + " from more than one process, the output of this simulation will be wrong!");
                } else if (!result.second) {
                    OpmLog::warning("Received consistently duplicated output data for well " + name + " from more than one process - this might be problematic!");
                }
            }
        }

        void init_json(Json::JsonObject& json_data) const {
            for (const auto& [wname, well] : *this) {
                auto json_well = json_data.add_object(wname);
                well.init_json(json_well);
            }
        }


        Json::JsonObject json() const {
            Json::JsonObject json_data;
            this->init_json(json_data);
            return json_data;
        }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(static_cast<std::map<std::string,Well>&>(*this));
        }

        static Wells serializationTestObject()
        {
            Wells w;
            w.insert({"test_well", Well::serializationTestObject()});

            return w;
        }
    };

    struct WellBlockAveragePressures
    {
        std::unordered_map<std::string, WellBlockAvgPress> values{};

        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const;

        template <class MessageBufferType>
        void read(MessageBufferType& buffer);

        bool operator==(const WellBlockAveragePressures& that) const
        {
            return this->values == that.values;
        }

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->values);
        }

        static WellBlockAveragePressures serializationTestObject()
        {
            return {
                { { "I-45", WellBlockAvgPress::serializationTestObject() } },
            };
        }
    };

    /* IMPLEMENTATIONS */

    inline bool Rates::has( opt m ) const {
        const auto mand = static_cast< enum_size >( this->mask )
                        & static_cast< enum_size >( m );

        return static_cast< opt >( mand ) == m;
    }

    inline double Rates::get( opt m ) const {
        if( !this->has( m ) )
            throw std::invalid_argument( "Uninitialized value." );

        return this->get_ref( m );
    }

    inline double Rates::get( opt m, double default_value ) const {
        if( !this->has( m ) ) return default_value;

        return this->get_ref( m );
    }

    inline double Rates::get( opt m, double default_value, const std::string& tracer_name) const {
        if( !this->has( m ) ) return default_value;

        if( m == opt::tracer && this->tracer.find(tracer_name) == this->tracer.end()) return default_value;

        return this->get_ref( m, tracer_name);
    }

    inline Rates& Rates::set( opt m, double value ) {
        this->get_ref( m ) = value;

        /* mask |= m */
        this->mask = static_cast< opt >(
                        static_cast< enum_size >( this->mask ) |
                        static_cast< enum_size >( m )
                    );

        return *this;
    }

    inline Rates& Rates::set( opt m, double value , const std::string& tracer_name ) {
        this->get_ref( m , tracer_name) = value;

        /* mask |= m */
        this->mask = static_cast< opt >(
                        static_cast< enum_size >( this->mask ) |
                        static_cast< enum_size >( m )
                    );

        return *this;
    }

    inline bool Rates::operator==(const Rates& rate) const
    {
      return mask == rate.mask &&
             wat == rate.wat &&
             oil == rate.oil &&
             gas == rate.gas &&
             polymer == rate.polymer &&
             solvent == rate.solvent &&
             energy == rate.energy &&
             dissolved_gas == rate.dissolved_gas &&
             vaporized_oil == rate.vaporized_oil &&
             reservoir_water == rate.reservoir_water &&
             reservoir_oil == rate.reservoir_oil &&
             reservoir_gas == rate.reservoir_gas &&
             productivity_index_water == rate.productivity_index_water &&
             productivity_index_gas == rate.productivity_index_gas &&
             productivity_index_oil == rate.productivity_index_oil &&
             well_potential_water == rate.well_potential_water &&
             well_potential_oil == rate.well_potential_oil &&
             well_potential_gas == rate.well_potential_gas &&
             brine == rate.brine &&
             alq == rate.alq &&
             tracer == rate.tracer &&
             microbial == rate.microbial &&
             oxygen == rate.oxygen &&
             urea == rate.urea &&
             vaporized_water == rate.vaporized_water &&
             mass_gas == rate.mass_gas &&
             mass_wat == rate.mass_wat;
    }


    /*
     * To avoid error-prone and repetitve work when extending rates with new
     * values, the get+set methods use this helper get_ref to determine what
     * member to manipulate. To add a new option, just add another case
     * corresponding to the enum entry in Rates to this function.
     *
     * This is an implementation detail and understanding this has no
     * significant impact on correct use of the class.
     */
    inline const double& Rates::get_ref( opt m ) const {
        switch( m ) {
            case opt::wat: return this->wat;
            case opt::oil: return this->oil;
            case opt::gas: return this->gas;
            case opt::polymer: return this->polymer;
            case opt::solvent: return this->solvent;
            case opt::energy: return this->energy;
            case opt::dissolved_gas: return this->dissolved_gas;
            case opt::vaporized_oil: return this->vaporized_oil;
            case opt::reservoir_water: return this->reservoir_water;
            case opt::reservoir_oil: return this->reservoir_oil;
            case opt::reservoir_gas: return this->reservoir_gas;
            case opt::productivity_index_water: return this->productivity_index_water;
            case opt::productivity_index_oil: return this->productivity_index_oil;
            case opt::productivity_index_gas: return this->productivity_index_gas;
            case opt::well_potential_water: return this->well_potential_water;
            case opt::well_potential_oil: return this->well_potential_oil;
            case opt::well_potential_gas: return this->well_potential_gas;
            case opt::brine: return this->brine;
            case opt::alq: return this->alq;
            case opt::tracer: /* Should _not_ be called with tracer argument */
                break;
            case opt::microbial: return this->microbial;
            case opt::oxygen: return this->oxygen;
            case opt::urea: return this->urea;
            case opt::vaporized_water: return this->vaporized_water;
            case opt::mass_gas: return this->mass_gas;
            case opt::mass_wat: return this->mass_wat;
        }

        throw std::invalid_argument(
                "Unknown value type '"
                + std::to_string( static_cast< enum_size >( m ) )
                + "'" );

    }

    inline const double& Rates::get_ref( opt m, const std::string& tracer_name ) const {
        if (m != opt::tracer)
            throw std::logic_error("Logic error - should be called with tracer argument");

        return this->tracer.at(tracer_name);
    }

    inline double& Rates::get_ref( opt m ) {
        return const_cast< double& >(
                static_cast< const Rates* >( this )->get_ref( m )
                );
    }

    inline double& Rates::get_ref( opt m, const std::string& tracer_name ) {
        if (m == opt::tracer) this->tracer.emplace(tracer_name, 0.0);
        return this->tracer.at(tracer_name);
    }

    void Rates::init_json(Json::JsonObject& json_data) const {

        if (this->has(opt::wat))
            json_data.add_item("wat", this->get(opt::wat));

        if (this->has(opt::oil))
            json_data.add_item("oil", this->get(opt::oil));

        if (this->has(opt::gas))
            json_data.add_item("gas", this->get(opt::gas));

    }

    bool inline Rates::flowing() const {
        return ((this->wat != 0) ||
                (this->oil != 0) ||
                (this->gas != 0));
    }

    inline bool Well::flowing() const noexcept {
        return this->rates.flowing();
    }

    template <class MessageBufferType>
    void Rates::write(MessageBufferType& buffer) const {
            buffer.write(this->mask);
            buffer.write(this->wat);
            buffer.write(this->oil);
            buffer.write(this->gas);
            buffer.write(this->polymer);
            buffer.write(this->solvent);
            buffer.write(this->energy);
            buffer.write(this->dissolved_gas);
            buffer.write(this->vaporized_oil);
            buffer.write(this->reservoir_water);
            buffer.write(this->reservoir_oil);
            buffer.write(this->reservoir_gas);
            buffer.write(this->productivity_index_water);
            buffer.write(this->productivity_index_oil);
            buffer.write(this->productivity_index_gas);
            buffer.write(this->well_potential_water);
            buffer.write(this->well_potential_oil);
            buffer.write(this->well_potential_gas);
            buffer.write(this->brine);
            buffer.write(this->alq);

            //tracer:
            unsigned int size = this->tracer.size();
            buffer.write(size);
            for (const auto& [name, rate] : this->tracer) {
                buffer.write(name);
                buffer.write(rate);
            }

            buffer.write(this->microbial);
            buffer.write(this->oxygen);
            buffer.write(this->urea);
            buffer.write(this->vaporized_water);
            buffer.write(this->mass_gas);
            buffer.write(this->mass_wat);
    }

    template <class MessageBufferType>
    void ConnectionFiltrate::write(MessageBufferType& buffer) const {
        buffer.write(this->rate);
        buffer.write(this->total);
        buffer.write(this->skin_factor);
        buffer.write(this->thickness);
        buffer.write(this->perm);
        buffer.write(this->poro);
        buffer.write(this->radius);
        buffer.write(this->area_of_flow);
    }

    template <class MessageBufferType>
    void Connection::write(MessageBufferType& buffer) const {
            buffer.write(this->index);
            this->rates.write(buffer);
            buffer.write(this->pressure);
            buffer.write(this->reservoir_rate);
            buffer.write(this->cell_pressure);
            buffer.write(this->cell_saturation_water);
            buffer.write(this->cell_saturation_gas);
            buffer.write(this->effective_Kh);
            buffer.write(this->trans_factor);
            buffer.write(this->d_factor);
            buffer.write(this->compact_mult);
            this->filtrate.write(buffer);
            this->fract.write(buffer);
    }

    void Connection::init_json(Json::JsonObject& json_data) const {
        auto json_rates = json_data.add_object("rates");
        this->rates.init_json(json_rates);

        json_data.add_item("global_index", static_cast<int>(this->index));
        json_data.add_item("pressure", this->pressure);
        json_data.add_item("reservoir_rate", this->reservoir_rate);
        json_data.add_item("cell_pressure", this->cell_pressure);
        json_data.add_item("swat", this->cell_saturation_water);
        json_data.add_item("sgas", this->cell_saturation_gas);
        json_data.add_item("Kh", this->effective_Kh);
        json_data.add_item("trans_factor", this->trans_factor);
        json_data.add_item("d_factor", this->d_factor);
        json_data.add_item("compact_mult", this->compact_mult);
    }

    template <class MessageBufferType>
    void Segment::write(MessageBufferType& buffer) const
    {
        buffer.write(this->segNumber);
        this->rates.write(buffer);
        this->pressures.write(buffer);
        this->velocity.write(buffer);
        this->holdup.write(buffer);
        this->viscosity.write(buffer);
        this->density.write(buffer);
    }

    template <class MessageBufferType>
    void CurrentControl::write(MessageBufferType& buffer) const
    {
        buffer.write(this->isProducer);
        if (this->isProducer) {
            buffer.write(this->prod);
        }
        else {
            buffer.write(this->inj);
        }
    }

    template <class MessageBufferType>
    void WellBlockAvgPress::write(MessageBufferType& buffer) const
    {
        for (const auto& quantity : this->wbp_) {
            buffer.write(quantity);
        }
    }

    template <class MessageBufferType>
    void WellFiltrate::write(MessageBufferType& buffer) const
    {
        buffer.write(this->rate);
        buffer.write(this->total);
        buffer.write(this->concentration);
    }

    template <class MessageBufferType>
    void Well::write(MessageBufferType& buffer) const
    {
        this->rates.write(buffer);

        buffer.write(this->bhp);
        buffer.write(this->thp);
        buffer.write(this->temperature);
        buffer.write(this->control);
        buffer.write(this->efficiency_scaling_factor);

        this->filtrate.write(buffer);

        {
            const auto status = ::Opm::WellStatus2String(this->dynamicStatus);
            buffer.write(status);
        }

        {
            const unsigned int size = this->connections.size();
            buffer.write(size);

            for (const Connection& comp : this->connections) {
                comp.write(buffer);
            }
        }

        {
            const auto nSeg =
                static_cast<unsigned int>(this->segments.size());
            buffer.write(nSeg);

            for (const auto& seg : this->segments) {
                seg.second.write(buffer);
            }
        }

        this->current_control.write(buffer);
        this->guide_rates.write(buffer);
        this->limits.write(buffer);
    }

    template <class MessageBufferType>
    void WellBlockAveragePressures::write(MessageBufferType& buffer) const
    {
        buffer.write(this->values.size());

        for (const auto& [well, value] : this->values) {
            buffer.write(well);
            value.write(buffer);
        }
    }

    template <class MessageBufferType>
    void Rates::read(MessageBufferType& buffer) {
            buffer.read(this->mask);
            buffer.read(this->wat);
            buffer.read(this->oil);
            buffer.read(this->gas);
            buffer.read(this->polymer);
            buffer.read(this->solvent);
            buffer.read(this->energy);
            buffer.read(this->dissolved_gas);
            buffer.read(this->vaporized_oil);
            buffer.read(this->reservoir_water);
            buffer.read(this->reservoir_oil);
            buffer.read(this->reservoir_gas);
            buffer.read(this->productivity_index_water);
            buffer.read(this->productivity_index_oil);
            buffer.read(this->productivity_index_gas);
            buffer.read(this->well_potential_water);
            buffer.read(this->well_potential_oil);
            buffer.read(this->well_potential_gas);
            buffer.read(this->brine);
            buffer.read(this->alq);

            //tracer:
            unsigned int size;
            buffer.read(size);
            for (size_t i = 0; i < size; ++i) {
                std::string tracer_name;
                buffer.read(tracer_name);
                double tracer_rate;
                buffer.read(tracer_rate);
                this->tracer.emplace(tracer_name, tracer_rate);
            }

            buffer.read(this->microbial);
            buffer.read(this->oxygen);
            buffer.read(this->urea);
            buffer.read(this->vaporized_water);
            buffer.read(this->mass_gas);
            buffer.read(this->mass_wat);
    }

    template <class MessageBufferType>
    void ConnectionFiltrate::read(MessageBufferType& buffer) {
        buffer.read(this->rate);
        buffer.read(this->total);
        buffer.read(this->skin_factor);
        buffer.read(this->thickness);
        buffer.read(this->perm);
        buffer.read(this->poro);
        buffer.read(this->radius);
        buffer.read(this->area_of_flow);
    }

   template <class MessageBufferType>
   void Connection::read(MessageBufferType& buffer) {
            buffer.read(this->index);
            this->rates.read(buffer);
            buffer.read(this->pressure);
            buffer.read(this->reservoir_rate);
            buffer.read(this->cell_pressure);
            buffer.read(this->cell_saturation_water);
            buffer.read(this->cell_saturation_gas);
            buffer.read(this->effective_Kh);
            buffer.read(this->trans_factor);
            buffer.read(this->d_factor);
            buffer.read(this->compact_mult);
            this->filtrate.read(buffer);
            this->fract.read(buffer);
   }

    template <class MessageBufferType>
    void Segment::read(MessageBufferType& buffer)
    {
        buffer.read(this->segNumber);
        this->rates.read(buffer);
        this->pressures.read(buffer);
        this->velocity.read(buffer);
        this->holdup.read(buffer);
        this->viscosity.read(buffer);
        this->density.read(buffer);
    }

    template <class MessageBufferType>
    void CurrentControl::read(MessageBufferType& buffer)
    {
        buffer.read(this->isProducer);
        if (this->isProducer) {
            buffer.read(this->prod);
        }
        else {
            buffer.read(this->inj);
        }
    }

    template <class MessageBufferType>
    void WellBlockAvgPress::read(MessageBufferType& buffer)
    {
        for (auto& quantity : this->wbp_) {
            buffer.read(quantity);
        }
    }

    template <class MessageBufferType>
    void WellFiltrate::read(MessageBufferType& buffer)
    {
        buffer.read(this->rate);
        buffer.read(this->total);
        buffer.read(this->concentration);
    }

    template <class MessageBufferType>
    void Well::read(MessageBufferType& buffer)
    {
        this->rates.read(buffer);

        buffer.read(this->bhp);
        buffer.read(this->thp);
        buffer.read(this->temperature);
        buffer.read(this->control);
        buffer.read(this->efficiency_scaling_factor);

        this->filtrate.read(buffer);

        {
            auto status = std::string{};
            buffer.read(status);
            this->dynamicStatus = ::Opm::WellStatusFromString(status);
        }

        // Connection information
        {
            unsigned int size = 0;
            buffer.read(size);

            this->connections.resize(size);
            for (auto& connection : this->connections) {
                connection.read(buffer);
            }
        }

        // Segment information (if applicable)
        const auto nSeg = [&buffer]() -> unsigned int
        {
            auto n = 0u;
            buffer.read(n);

            return n;
        }();

        for (auto segID = 0*nSeg; segID < nSeg; ++segID) {
            auto seg = Segment{};
            seg.read(buffer);

            const auto segNumber = seg.segNumber;
            this->segments.emplace(segNumber, std::move(seg));
        }

        this->current_control.read(buffer);
        this->guide_rates.read(buffer);
        this->limits.read(buffer);
    }

    template <class MessageBufferType>
    void WellBlockAveragePressures::read(MessageBufferType& buffer)
    {
        const auto numWells = [&buffer, this]()
        {
            auto size = 0*this->values.size();
            buffer.read(size);

            return size;
        }();

        auto wellName = std::string{};
        for (auto well = 0*numWells; well < numWells; ++well) {
            buffer.read(wellName);

            this->values[wellName].read(buffer);
        }
    }

    void Well::init_json(Json::JsonObject& json_data) const {
        auto json_connections = json_data.add_array("connections");
        for (const auto& conn : this->connections) {
            auto json_conn = json_connections.add_object();
            conn.init_json(json_conn);
        }
        auto json_rates = json_data.add_object("rates");
        this->rates.init_json(json_rates);

        json_data.add_item("bhp", this->bhp);
        json_data.add_item("thp", this->thp);
        json_data.add_item("temperature", this->temperature);
        json_data.add_item("status", ::Opm::WellStatus2String(this->dynamicStatus));

        auto json_control = json_data.add_object("control");
        this->current_control.init_json(json_control);

        auto json_guiderate = json_data.add_object("guiderate");
        this->guide_rates.init_json(json_guiderate);
    }

}} // Opm::data

#endif //OPM_OUTPUT_WELLS_HPP
