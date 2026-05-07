/*
    Copyright 2026 NORCE.

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
#ifndef OPM_SALTARRAY_HPP
#define OPM_SALTARRAY_HPP

#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <opm/material/components/CaIon.hpp>
#include <opm/material/components/ClIon.hpp>
#include <opm/material/components/KIon.hpp>
#include <opm/material/components/MgIon.hpp>
#include <opm/material/components/NaIon.hpp>
#include <opm/material/components/SO4Ion.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <numeric>

namespace Opm
{
enum class SaltIndex : std::size_t
{
    NA, K, CA, MG, CL, SO4,
    NCOMP // DO NOT CHANGE! Must be the last enum to set size of SaltArray!
};

enum class IonType
{
    Cation, Anion
};

inline SaltIndex
saltIndexfromString(const std::string& s)
{
    if (s == "NA") {
        return SaltIndex::NA;
    }
    if (s == "K") {
        return SaltIndex::K;
    }
    if (s == "CA") {
        return SaltIndex::CA;
    }
    if (s == "MG") {
        return SaltIndex::MG;
    }
    if (s == "CL") {
        return SaltIndex::CL;
    }
    if (s == "SO4") {
        return SaltIndex::SO4;
    }
    throw std::invalid_argument("SaltIndex not valid!");
}

inline IonType
categorizeIon(const SaltIndex s)
{
    switch (s) {
    case SaltIndex::CA:
    case SaltIndex::NA:
    case SaltIndex::MG:
    case SaltIndex::K:
        return IonType::Cation;
    case SaltIndex::CL:
    case SaltIndex::SO4:
        return IonType::Anion;
    default:
        throw std::runtime_error("Unknown SaltIndex");
    }
}

inline short
ionStrength(const SaltIndex s)
{
    switch (s) {
    case SaltIndex::CA:
        return 2;
    case SaltIndex::NA:
        return 1;
    case SaltIndex::MG:
        return 2;
    case SaltIndex::K:
        return 1;
    case SaltIndex::CL:
        return -1;
    case SaltIndex::SO4:
        return -2;
    default:
        throw std::runtime_error("Unknown SaltIndex");
    }
}

template <class Scalar>
Scalar
saltMolarMass(const SaltIndex ind)
{
    switch (ind) {
    case SaltIndex::NA:
        return NaIon<Scalar>::molarMass();
    case SaltIndex::K:
        return KIon<Scalar>::molarMass();
    case SaltIndex::CA:
        return CaIon<Scalar>::molarMass();
    case SaltIndex::CL:
        return ClIon<Scalar>::molarMass();
    case SaltIndex::MG:
        return MgIon<Scalar>::molarMass();
    case SaltIndex::SO4:
        return SO4Ion<Scalar>::molarMass();
    default:
        throw std::invalid_argument("SaltIndex not valid!");
    }
}

template <class T>
class SaltArray
{
public:
    static constexpr std::size_t ncomp = static_cast<std::size_t>(SaltIndex::NCOMP);
    using SaltIndexPair = std::pair<std::vector<SaltIndex>, std::vector<SaltIndex> >;

    SaltArray() = default;

    void assign(const DeckRecord& record)
    {
        assert(record.size() == saltData_.size());
        T sum = 1.0;
        for (const auto& item : record) {
            auto index = saltIndexfromString(item.name());
            auto mMxMolal = saltMolarMass<T>(index) * static_cast<T>(item.get<double>(0));
            (*this)[index] = mMxMolal;
            sum += mMxMolal;
        }
        for (auto& elem : saltData_) {
            elem /= sum;
        }
    }

    const T& operator[](const SaltIndex ind) const
    {
        return saltData_[static_cast<std::size_t>(ind)];
    }

    T& operator[](const SaltIndex ind)
    {
        return saltData_[static_cast<std::size_t>(ind)];
    }

    bool operator==(const SaltArray& other) const
    {
        return saltData_ == other.saltData_;
    }

    SaltArray& operator=(const SaltArray& other) = default;

    template <class U>
    SaltArray& operator=(const SaltArray<U>& other)
    {
        std::transform(other.begin(),
                       other.end(),
                       saltData_.begin(),
                       [](const U& val) {
                           return static_cast<T>(val);
                       });
        return *this;
    }

    auto begin()
    {
        return saltData_.begin();
    }

    auto end()
    {
        return saltData_.end();
    }

    [[nodiscard]] auto begin() const
    {
        return saltData_.begin();
    }

    [[nodiscard]] auto end() const
    {
        return saltData_.end();
    }

    [[nodiscard]] constexpr std::size_t size() const
    {
        return saltData_.size();
    }

    T sum() const
    {
        return std::accumulate(begin(), end(), T{});
    }

    void clear() noexcept
    {
        saltData_.fill(T{});
    }

    [[nodiscard]] bool any_nonzero() const noexcept
    {
        return std::any_of(begin(),
                           end(),
                           [](const T& val) {
                               return val != T{};
                           });
    }

    [[nodiscard]] SaltIndexPair cation_anion_pair() const
    {
        std::vector<std::pair<short, SaltIndex> > cations;
        std::vector<std::pair<short, SaltIndex> > anions;
        cations.reserve(size());
        anions.reserve(size());
        for (std::size_t i = 0; i < size(); ++i) {
            auto ind = static_cast<SaltIndex>(i);
            auto ionStr = ionStrength(ind);
            if ((*this)[ind] != T{}) {
                if (categorizeIon(ind) == IonType::Cation) {
                    cations.emplace_back(ionStr, ind);

                } else {
                    anions.emplace_back(ionStr, ind);
                }
            }
        }

        // Sort by ion strength: descending for cations, ascending for anions (most to less
        // negative)
        std::ranges::sort(cations, std::ranges::greater{}, &std::pair<short, SaltIndex>::first);
        std::ranges::sort(anions, {}, &std::pair<short, SaltIndex>::first);

        // Return only cations and anions SaltIndex
        std::vector<SaltIndex> cationIndex;
        std::vector<SaltIndex> anionIndex;
        cationIndex.reserve(cations.size());
        anionIndex.reserve(anions.size());
        std::ranges::transform(cations,
                               std::back_inserter(cationIndex),
                               &std::pair<short, SaltIndex>::second);
        std::ranges::transform(anions,
                               std::back_inserter(anionIndex),
                               &std::pair<short, SaltIndex>::second);

        return {cationIndex, anionIndex};
    }

    SaltArray<T> to_molality() const
    {
        SaltArray<T> molalityArray;
        for (std::size_t i = 0; i < saltData_.size(); ++i) {
            auto sIdx = static_cast<SaltIndex>(i);
            molalityArray[sIdx] = (*this)[sIdx] / saltMolarMass<T>(sIdx);
        }

        T s = 1.0 - sum();
        for (auto& elem : molalityArray) {
            elem /= s;
        }
        return molalityArray;
    }

private:
    std::array<T, ncomp> saltData_{};
};

}
#endif //OPM_SALTARRAY_HPP
