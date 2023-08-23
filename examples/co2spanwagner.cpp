// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
This file is part of the Open Porous Media project (OPM).

OPM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

OPM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OPM.  If not, see <http://www.gnu.org/licenses/>.

Consult the COPYING file in the top-level source directory of this
module for the precise wording of the license and the list of
copyright holders.
*/
/*!
* \file
*
* \brief Simple CLI program to calculate properties in CO2SpanWagner class
*
*/
#include "config.h"
#include <opm/material/components/CO2SpanWagner.hpp>

#include <iostream>

int main(int argc, char **argv)
{
    // Display help
    bool help = false;
    for (int i = 1; i < argc; ++i) {
        std::string tmp = argv[i];
        help = help || (tmp  == "--h") || (tmp  == "--help");
    }
    if (argc < 4 || help) {
        std::cout << "usage: co2spanwagner prop p t" << std::endl;
        std::cout << std::endl;
        std::cout << "arguments : " << std::endl;
        std::cout << "prop\t\tProperty = {density [kg/m3], enthalpy [J/kg], viscosity [Pa*s]}" << std::endl;
        std::cout << "p   \t\tPressure [Pa]" << std::endl;
        std::cout << "t   \t\tTemperature [K]" << std::endl;
        return EXIT_FAILURE;
    }

    // Parse input arguments
    std::string prop = argv[1];
    double p = atof(argv[2]);
    double T = atof(argv[3]);

    // Span-Wagner class
    using CO2 = Opm::CO2SpanWagner<double>;

    // Calculate property
    double value;
    if (prop == "density") {
        value = CO2::gasDensity(T, p);
    }
    else if (prop == "enthalpy") {
        value = CO2::gasEnthalpy(T, p);
    }
    else if (prop == "viscosity") {
        value = CO2::gasViscosity(T, p);
    }
    else {
        throw std::runtime_error("prop {" + prop + "} not recognized!");
    }
    // Print result
    std::cout << value << std::endl;

    return 0;

}