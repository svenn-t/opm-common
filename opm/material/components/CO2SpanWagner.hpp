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
* \ingroup Components
*
* \copydoc Opm::CO2SpanWagner
*
*/
#ifndef OPM_CO2_SPAN_WAGNER_HPP
#define OPM_CO2_SPAN_WAGNER_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/densead/Math.hpp>

namespace Opm {

template <class Scalar>
class CO2SpanWagner : public Component<Scalar, CO2SpanWagner<Scalar> >
{
    using IdealGas = Opm::IdealGas<Scalar>;
public:
    /*!
    * \brief A human readable name for the \f$H_2\f$.
    */
    static std::string name()
    { return "CO2"; }

    /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$.
    */
    static constexpr Scalar molarMass()
    { return 44.009e-3; }

    /*!
    * \brief Returns the critical temperature \f$\mathrm{[K]}\f$.
    */
    static Scalar criticalTemperature()
    { return 304.128; }

    /*!
    * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$.
    */
    static Scalar criticalPressure()
    { return 7.377e6; }

    /*!
    * \brief Returns the critical density \f$\mathrm{[kg/m^3]}\f$.
    */
    static Scalar criticalDensity()
    { return 467.6; }

    /*!
    * \brief Returns the temperature \f$\mathrm{[K]}\f$ at triple point.
    */
    static Scalar tripleTemperature()
    { return 216.592; }

    /*!
    * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at triple point.
    */
    static Scalar triplePressure()
    { return 0.51795e6; }

    /*!
     * \brief Critical volume in \f$\mathrm{[m3/kmol]}\f$.
     */
    static Scalar criticalVolume()
    { return 9.412e-5; }

    /*!
    * \brief The vapor pressure in [N/m^2] of pure CO2
    *        at a given temperature.
    *
    * See:
    *
    * R. Span and W. Wagner: A New Equation of State for Carbon
    * Dioxide Covering the Fluid Region from the Triple‐Point
    * Temperature to 1100 K at Pressures up to 800 MPa. Journal of
    * Physical and Chemical Reference Data, 25 (6), pp. 1509-1596,
    * 1996
    */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& T)
    {
        static constexpr Scalar a[4] =
            { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
        static constexpr Scalar t[4] =
            { 1.0, 1.5, 2.0, 4.0 };

        // this is on page 1524 of the reference
        Evaluation exponent = 0;
        Evaluation Tred = T/criticalTemperature();
        for (int i = 0; i < 4; ++i)
            exponent += a[i]*pow(1 - Tred, t[i]);
        exponent *= 1.0/Tred;

        return exp(exponent)*criticalPressure();
    }

    /*!
    * \brief Returns true iff the gas phase is assumed to be compressible
    */
    static bool gasIsCompressible()
    { return true; }

    /*!
    * \brief Returns true iff the gas phase is assumed to be ideal
    */
    static bool gasIsIdeal()
    { return false; }

    /*!
    * \brief Calculate reduced density (rho/rho_crit) from pressure and temperature
    * 
    * \param pg gas phase pressure [Pa]
    * \param temperature temperature [K]
    */
    template <class Evaluation> 
    static Evaluation reducedDensity(const Evaluation& temperature, const Evaluation& pg)
    {
        // Generate a bracket for bisection
        auto [rho_red_min, rho_red_max] = reducedDensityBracket_(temperature, pg);

        // Bisection loop
        Evaluation rho_red;
        Evaluation fmid;
        Evaluation fmin = rootFindingObj_(rho_red_min, temperature, pg);
        for (int iteration = 1; iteration < 100; ++iteration) {
            // New midpoint and its obj. value
            rho_red = (rho_red_min + rho_red_max) / 2;
            fmid = rootFindingObj_(rho_red, temperature, pg);

            // Check if midpoint fulfills f=0 or interval is sufficiently small
            if (Opm::abs(fmid) < 1e-6 || Opm::abs((rho_red_max - rho_red_min) / 2) < 1e-6) {
                return rho_red;
            }

            // Else we repeat with midpoint being either xmin or xmax (depending on the signs)
            else if ((Opm::getValue(fmid) > 0.0 && Opm::getValue(fmin) < 0.0) ||
                (Opm::getValue(fmid) < 0.0 && Opm::getValue(fmin) > 0.0)) {
                // fmid has same sign as fmax so we set xmid as the new xmax
                rho_red_max = rho_red;
            }
            else {
                // fmid has same sign as fmin so we set xmid as the new xmin
                rho_red_min = rho_red;
                fmin = fmid;
            }
        }
        throw std::runtime_error("No reduced density could be found for current pressure using bisection!");
    }

    /*!
    * \brief The density of CO2 at a given pressure and temperature [kg/m^3].
    */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature,
                                 const Evaluation& pressure)
    {
        return reducedDensity(temperature, pressure) * criticalDensity();
    }

    /*!
    * \brief Specific internal energy of CO2 [J/kg].
    */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        // Intermediate calculations
        Evaluation rho_red = reducedDensity(temperature, pressure);
        Evaluation T_red = criticalTemperature() / temperature;
        Scalar R = IdealGas::R;

        Evaluation dphi0_dTred = derivIdealHelmholtzWrtRecipRedTemp(T_red);
        Evaluation dphir_dTred = derivResHelmholtzWrtRecipRedTemp(T_red, rho_red);

        // Table 3
        return  R * criticalTemperature() * (dphi0_dTred + dphir_dTred) / molarMass();
    }

    /*!
    * \brief Specific enthalpy of gaseous CO2 [J/kg].
    */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    {
        // Intermediate calculations
        Evaluation rho_red = reducedDensity(temperature, pressure);
        Evaluation T_red = criticalTemperature() / temperature;
        Scalar R = IdealGas::R;

        Evaluation dphi0_dTred = derivIdealHelmholtzWrtRecipRedTemp(T_red);
        Evaluation dphir_dTred = derivResHelmholtzWrtRecipRedTemp(T_red, rho_red);
        Evaluation dphir_dRhoRed = derivResHelmholtzWrtRedRho(T_red, rho_red);

        // Table 3
        return R * temperature * (1 + T_red * (dphi0_dTred + dphir_dTred) + rho_red * dphir_dRhoRed) / molarMass();
    }

    /*!
    * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of CO2. This is equivalent to the
    * partial derivative of the specific enthalpy wrt temperature.
    * 
    * \param temperature temperature of component in \f$\mathrm{[K]}\f$
    * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    *
    */
    template <class Evaluation>
    static const Evaluation gasHeatCapacity(Evaluation temperature,
                                            Evaluation pressure)
    {
        // Reduced variables
        Evaluation rho_red = reducedDensity(temperature, pressure);
        Evaluation T_red = criticalTemperature() / temperature;
        
        // Intermediate calculations
        Evaluation dphir_dRhoRed = derivResHelmholtzWrtRedRho(T_red, rho_red);
        Evaluation d2phir_dTred_dRhoRed = secDerivResHelmholtzWrtRecipRedTempAndRedRho(T_red, rho_red);
        Evaluation d2phir_d2RhoRed = secDerivResHelmholtzWrtRedRho(T_red, rho_red);
        Evaluation d2phi0_d2Tred = secDerivIdealHelmholtzWrtRecipRedTemp(T_red);
        Evaluation d2phir_d2Tred = secDerivResHelmholtzWrtRecipRedTemp(T_red, rho_red);
        Scalar R = IdealGas::R;
        
        Evaluation numerator = 1 + rho_red * dphir_dRhoRed - rho_red * T_red * d2phir_dTred_dRhoRed;
        Evaluation denominator = 1 + 2 * rho_red * dphir_dRhoRed + rho_red * rho_red * d2phir_d2RhoRed;

        // Table 3
        return R * (-T_red * T_red * (d2phi0_d2Tred + d2phir_d2Tred) + (numerator * numerator / denominator)) / molarMass();  // [J/(kg*K)]
    }

    /*!
    * \brief Specific isochoric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of CO2.
    * 
    * \param temperature temperature of component in \f$\mathrm{[K]}\f$
    * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    *
    */
    template <class Evaluation>
    static const Evaluation gasIsochoricHeatCapacity(Evaluation temperature,
                                                     Evaluation pressure)
    {
        // Reduced variables
        Evaluation rho_red = reducedDensity(temperature, pressure);
        Evaluation T_red = criticalTemperature() / temperature;

        // Intermediate calculations
        Evaluation d2phi0_d2Tred = secDerivIdealHelmholtzWrtRecipRedTemp(T_red);
        Evaluation d2phir_d2Tred = secDerivResHelmholtzWrtRecipRedTemp(T_red, rho_red);

        // Table 3
        Scalar R = IdealGas::R;
     
        return -R * T_red * T_red * (d2phi0_d2Tred + d2phir_d2Tred) / molarMass();  // [J/(mol*K)]
    }

    /*!
    * \brief The dynamic viscosity [Pa s] of CO2.
    *
    * Equations given in: - Vesovic et al., 1990
    *                        - Fenhour etl al., 1998
    */
    template <class Evaluation>
    static Evaluation gasViscosity(Evaluation temperature,
                                   const Evaluation& pressure)
    {
        constexpr Scalar a0 = 0.235156;
        constexpr Scalar a1 = -0.491266;
        constexpr Scalar a2 = 5.211155e-2;
        constexpr Scalar a3 = 5.347906e-2;
        constexpr Scalar a4 = -1.537102e-2;

        constexpr Scalar d11 = 0.4071119e-2;
        constexpr Scalar d21 = 0.7198037e-4;
        constexpr Scalar d64 = 0.2411697e-16;
        constexpr Scalar d81 = 0.2971072e-22;
        constexpr Scalar d82 = -0.1627888e-22;

        constexpr Scalar ESP = 251.196;

        if(temperature < 275.) // regularization
            temperature = 275.0;
        Evaluation TStar = temperature/ESP;

        // mu0: viscosity in zero-density limit
        const Evaluation logTStar = log(TStar);
        Evaluation SigmaStar = exp(a0 + logTStar*(a1 + logTStar*(a2 + logTStar*(a3 + logTStar*a4))));

        Evaluation mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        const Evaluation rho = gasDensity(temperature, pressure); // CO2 mass density [kg/m^3]

        // dmu : excess viscosity at elevated density
        Evaluation dmu =
            d11*rho
            + d21*rho*rho
            + d64*pow(rho, 6.0)/(TStar*TStar*TStar)
            + d81*pow(rho, 8.0)
            + d82*pow(rho, 8.0)/TStar;

        return (mu0 + dmu)/1.0e6; // conversion to [Pa s]
    }

private:
    // 
    // Parameter values need in the ideal-gas contribution to the reduced Helmholtz free energy (Table 27)
    // 
    static constexpr Scalar a0_[8] = {8.37304456, -3.70454304, 2.5, 1.99427042, 0.62105248, 0.41195293, 1.04028922, 
        0.08327678};
    static constexpr Scalar theta0_[5] = {3.15163, 6.11190, 6.77708, 11.32384, 27.08792};

    // 
    // Parameter values needed in the residual contribution to the reduced Helmholtz free energy (Table 31)
    // 
    // Main parameters
    static constexpr Scalar n_[42] = {0.38856823203161, 2.9385475942740, -5.5867188534934, -0.76753199592477, 
        0.31729005580416, 0.54803315897767, 0.12279411220335, 2.1658961543220, 1.5841735109724, -0.23132705405503, 
        0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678, 
        -0.12175860225246, -0.37055685270086, -0.016775879700426, -0.11960736637987, -0.045619362508778, 
        0.035612789270346, -0.0074427727132052, -0.0017395704902432, -0.021810121289527, 0.024332166559236, 
        -0.037440133423463, 0.14338715756878, -0.13491969083286, -0.023151225053480, 0.012363125492901,
        0.0021058321972940, -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320, 
        26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179, -0.66642276540751, 0.72608632349897, 
        0.055068668612842}; 
    static constexpr Scalar t_[39] = {0., 0.75, 1., 2., 0.75, 2., 0.75, 1.5, 1.5, 2.5, 0., 1.5, 2., 0., 1., 2., 3., 6., 
        3., 6., 8., 6., 0., 7., 12., 16., 22., 24., 16., 24., 8., 2., 28., 14., 1., 0., 1., 3., 3.};
    static constexpr Scalar d_[39] = {1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 
        5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3};
    static constexpr Scalar c_[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6};
    
    // Gaussian term
    static constexpr Scalar alpha_[5] = {25., 25., 25., 15., 20.};
    static constexpr Scalar beta_[5] = {325., 300., 300., 275., 275.};
    static constexpr Scalar epsilon_[5] = {1., 1., 1., 1., 1.};
    static constexpr Scalar gamma_[5] = {1.16, 1.19, 1.19, 1.25, 1.22};

    // Modified polynomial term
    static constexpr Scalar A_[3] = {0.7, 0.7, 0.7};
    static constexpr Scalar B_[3] = {0.3, 0.3, 1.0};
    static constexpr Scalar C_[3] = {10., 10., 12.5};
    static constexpr Scalar D_[3] = {275., 275., 275.};
    static constexpr Scalar a_[3] = {3.5, 3.5, 3.0};
    static constexpr Scalar b_[3] = {0.875, 0.925, 0.875};
    static constexpr Scalar betap_[3] = {0.3, 0.3, 0.3};  // beta in table

    // /////////////////////////
    // IDEAL-GAS PART HELMHOLTZ
    // /////////////////////////
    /*!
    * \brief The ideal-gas part of Helmholtz energy. Eq. (6.3) in Span & Wagner (1996).
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation idealGasHelmholtz_(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // Terms not in sum
        Evaluation s1 = log(rho_red) + a0_[0] + a0_[1] * T_red + a0_[2] * log(T_red);

        // Terms in sum
        Evaluation s2 = 0.0;
        for (int i = 3; i < 8; ++i) {
            s2 += a0_[i] * log(1 - exp(-theta0_[i-3] * T_red));
        }

        // Return total
        Evaluation s = s1 + s2;
        return s;
    }

    // ///////////////////////////////////////////
    // DERIVATIVES OF IDEAL HELMOLTZ (TABLE 28)
    // //////////////////////////////////////////
    /*!
    * \brief Derivative of the ideal-gas part of Helmholtz energy wrt to reciprocal reduced temperature.
    * 
    * \param T_red reciprocal reduced temperature [-]
    */
    template <class Evaluation> 
    static Evaluation derivIdealHelmholtzWrtRecipRedTemp(const Evaluation& T_red)
    {
        // Terms not in sum
        Evaluation s1 = a0_[1] + (a0_[2] / T_red);

        // Terms in sum
        Evaluation s2 = 0.0;
        for (int i = 3; i < 8; ++i) {
            s2 += a0_[i] * theta0_[i-3] * (1 / (1 - exp(-theta0_[i-3] * T_red)) - 1);
        }

        // Return total
        Evaluation s = s1 + s2;
        return s;
    }

    /*!
    * \brief Second derivative of the ideal-gas part of Helmholtz energy wrt to reciprocal reduced temperature.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation secDerivIdealHelmholtzWrtRecipRedTemp(const Evaluation& T_red)
    {
        // Terms not in sum
        Evaluation s1 = -a0_[2] / (T_red * T_red);
        
        // Terms in sum
        Evaluation s2 = 0.0;
        Evaluation interm;
        for (int i = 3; i < 8; ++i) {
            // Intermediate calculation
            interm = 1 - exp(-theta0_[i-3] * T_red);

            s2 += a0_[i] * (theta0_[i-3] * theta0_[i-3]) * exp(-theta0_[i-3] * T_red) / (interm * interm);
        }

        // Return total
        Evaluation s = s1 - s2;
        return s;
    }

    // ////////////////////////////////////
    // RESIDUAL PART HELMHOLTZ (TABLE 32)
    // ////////////////////////////////////
    /*!
    * \brief The residual part of Helmholtz energy. Eq. (6.5) in Span & Wagner (1996).
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation residualHelmholtz_(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 += n_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i]);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * pow(T_red, t_[i]) * pow(rho_red, d_[i]) * exp(-pow(rho_red, c_[i-7]));
        }

        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            s3 += n_[i] * pow(T_red, t_[i]) * pow(rho_red, d_[i]) * exp_term;
        }

        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation Psi;
        for (int i = 39; i < 42; ++i) {
            // Delta term
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);

            // Sum
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            s4 += n_[i] * rho_red * pow(Delta, b_[i-39]) * Psi;
        }
        
        // Return total sum
        return s1 + s2 + s3 + s4;
    }

    // /////////////////////////////////////////////
    // DERIVATIVES OF RESIDUAL HELMOLTZ (TABLE 32)
    // /////////////////////////////////////////////
    /*!
    * \brief Derivative of the residual part of Helmholtz energy wrt. reduced density.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation derivResHelmholtzWrtRedRho(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term 
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 += d_[i] * n_[i] * pow(rho_red, d_[i] - 1) * pow(T_red, t_[i]);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * pow(T_red, t_[i]) * pow(rho_red, d_[i] - 1) * exp(-pow(rho_red, c_[i-7])) *
                (d_[i] - c_[i-7] * pow(rho_red, c_[i-7]));
        }
        
        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            s3 += n_[i] * pow(T_red, t_[i]) * pow(rho_red, d_[i]) * exp_term *
                ((d_[i] / rho_red) - 2 * alpha_[i-34] * (rho_red - epsilon_[i-34]));
        }

        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation dDeltaPowB_ddelta;
        Evaluation Psi;
        Evaluation dPsi_ddelta;
        for (int i = 39; i < 42; ++i) {
            // Delta function terms
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);
            dDeltaPowB_ddelta = dDistPowB_dRhoRed(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);

            // Psi function terms
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_ddelta = dPsi_dRhoRed(T_red, rho_red, C_[i-39], D_[i-39]);

            // Sum
            s4 += n_[i] * (pow(Delta, b_[i-39]) * (Psi + rho_red * dPsi_ddelta) + dDeltaPowB_ddelta * rho_red * Psi);
        }

        // Return total sum
        return s1 + s2 + s3 + s4;
    }

    /*!
    * \brief Second derivative of the residual part of Helmholtz energy wrt. reduced density.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation secDerivResHelmholtzWrtRedRho(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term 
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 += d_[i] * (d_[i] - 1) * n_[i] * pow(rho_red, d_[i] - 2) * pow(T_red, t_[i]);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * exp(-pow(rho_red, c_[i-7])) * pow(rho_red, d_[i] - 2) * pow(T_red, t_[i]) * 
                ((d_[i] - c_[i-7] * pow(rho_red, c_[i-7])) * (d_[i] - 1 - c_[i-7] * pow(rho_red, c_[i-7])) - 
                    c_[i-7] * c_[i-7] * pow(rho_red, c_[i-7]));
        }

        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            // Sum
            s3 += n_[i] * pow(T_red, t_[i]) * exp_term * pow(rho_red, d_[i]) * (
                -2 * alpha_[i-34] + 
                4 * alpha_[i-34] * alpha_[i-34] * (rho_red - epsilon_[i-34]) * (rho_red - epsilon_[i-34]) -
                4 * d_[i] * alpha_[i-34] * (1 / rho_red) * (rho_red - epsilon_[i-34]) +
                d_[i] * (d_[i] - 1) * (1 / (rho_red * rho_red))
                );
        }

        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation dDeltaPowB_ddelta;
        Evaluation d2DeltaPowB_d2delta;
        Evaluation Psi;
        Evaluation dPsi_ddelta;
        Evaluation d2Psi_d2delta;
        for (int i = 39; i < 42; ++i) {
            // Delta function terms
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);
            dDeltaPowB_ddelta = dDistPowB_dRhoRed(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);
            d2DeltaPowB_d2delta = ddDistPowB_ddRhoRed(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);

            // Psi function terms
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_ddelta = dPsi_dRhoRed(T_red, rho_red, C_[i-39], D_[i-39]);
            d2Psi_d2delta = ddPsi_ddRhoRed(T_red, rho_red, C_[i-39], D_[i-39]);

            // Sum
            s4 += n_[i] * (
                pow(Delta, b_[i-39]) * (2 * dPsi_ddelta + rho_red * d2Psi_d2delta) +
                2 * dDeltaPowB_ddelta * (Psi + rho_red * dPsi_ddelta) +
                d2DeltaPowB_d2delta * rho_red * Psi
                );
        }

        // Return total sum
        Evaluation s = s1 + s2 + s3 + s4;
        return s;
    }

    /*!
    * \brief Derivative of the residual part of Helmholtz energy wrt. reciprocal reduced temperature.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation derivResHelmholtzWrtRecipRedTemp(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term 
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 += n_[i] * t_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i] - 1);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * t_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i] - 1) * exp(-pow(rho_red, c_[i-7]));
        }

        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            s3 += n_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i]) * exp_term *
                ((t_[i] / T_red) - 2 * beta_[i-34] * (T_red - gamma_[i-34]));
        }

        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation dDeltaPowB_dtau;
        Evaluation Psi;
        Evaluation dPsi_dtau;
        for (int i = 39; i < 42; ++i) {
            // Delta function terms
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);
            dDeltaPowB_dtau = dDistPowB_dTred(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);

            // Psi function terms
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_dtau = dPsi_dTred(T_red, rho_red, C_[i-39], D_[i-39]);

            // Sum
            s4 += n_[i] * rho_red * (dDeltaPowB_dtau * Psi + pow(Delta, b_[i-39]) * dPsi_dtau);
        }

        // Return total sum
        return s1 + s2 + s3 + s4;
    }

    /*!
    * \brief Second derivative of the residual part of Helmholtz energy wrt. reciprocal reduced temperature.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation secDerivResHelmholtzWrtRecipRedTemp(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term 
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 +=  n_[i] * t_[i] * (t_[i] - 1) * pow(rho_red, d_[i]) * pow(T_red, t_[i] - 2);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * t_[i] * (t_[i] - 1) * pow(rho_red, d_[i]) * pow(T_red, t_[i] - 2) * exp(-pow(rho_red, c_[i-7]));
        }

        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        Evaluation interm;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            // Intermediate calculation
            interm = (t_[i] / T_red) - 2 * beta_[i-34] * (T_red - gamma_[i-34]);

            // Sum
            s3 += n_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i]) * exp_term *
                (interm * interm - (t_[i] / (T_red * T_red)) - 2 * beta_[i-34]);
        }

        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation dDeltaPowB_dtau;
        Evaluation d2DeltaPowB_d2tau;
        Evaluation Psi;
        Evaluation dPsi_dtau;
        Evaluation d2Psi_d2tau;
        for (int i = 39; i < 42; ++i) {
            // Delta function terms
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);
            dDeltaPowB_dtau = dDistPowB_dTred(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);
            d2DeltaPowB_d2tau = ddDistPowB_ddTred(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);

            // Psi function terms
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_dtau = dPsi_dTred(T_red, rho_red, C_[i-39], D_[i-39]);
            d2Psi_d2tau = ddPsi_ddTred(T_red, rho_red, C_[i-39], D_[i-39]);

            // Sum
            s4 += n_[i] * rho_red * (
                d2DeltaPowB_d2tau * Psi + 
                2 * dDeltaPowB_dtau * dPsi_dtau + 
                pow(Delta, b_[i-39]) * d2Psi_d2tau
            );
        }

        // Return total sum
        return s1 + s2 + s3 + s4;
    }

    /*!
    * \brief Second derivative of the residual part of Helmholtz energy wrt. reciprocal reduced temperature and 
    * reduced density.
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    template <class Evaluation> 
    static Evaluation secDerivResHelmholtzWrtRecipRedTempAndRedRho(const Evaluation& T_red, const Evaluation& rho_red)
    {
        // First sum term 
        Evaluation s1 = 0.0;
        for (int i = 0; i < 7; ++i) {
            s1 += n_[i] * d_[i] * t_[i] * pow(rho_red, d_[i] - 1) * pow(T_red, t_[i] - 1);
        }

        // Second sum term
        Evaluation s2 = 0.0;
        for (int i = 7; i < 34; ++i) {
            s2 += n_[i] * exp(-pow(rho_red, c_[i-7])) * pow(rho_red, d_[i] - 1) * t_[i] * pow(T_red, t_[i] - 1) * 
                (d_[i] - c_[i-7] * pow(rho_red, c_[i-7]));
        }
        
        // Third sum term
        Evaluation s3 = 0.0;
        Evaluation exp_term;
        for (int i = 34; i < 39; ++i) {
            // Exponential term
            exp_term = expTerm_(T_red, rho_red, alpha_[i-34], beta_[i-34], gamma_[i-34], epsilon_[i-34]);

            // Sum
            s3 += n_[i] * pow(rho_red, d_[i]) * pow(T_red, t_[i]) * exp_term * 
                ((d_[i] / rho_red) - 2 * alpha_[i-34] * (rho_red - epsilon_[i-34])) * 
                ((t_[i] / T_red) - 2 * beta_[i-34] * (T_red - gamma_[i-34]));
        }
        
        // Fourth sum term
        Evaluation s4 = 0.0;
        Evaluation Delta;
        Evaluation dDeltaPowB_ddelta;
        Evaluation dDeltaPowB_dtau;
        Evaluation d2DeltaPowB_ddelta_dtau;
        Evaluation Psi;
        Evaluation dPsi_ddelta;
        Evaluation dPsi_dtau;
        Evaluation d2Psi_ddelta_dtau;
        for (int i = 39; i < 42; ++i) {
            // Delta function terms
            Delta = distFunc_(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], betap_[i-39]);
            dDeltaPowB_ddelta = dDistPowB_dRhoRed(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);
            dDeltaPowB_dtau = dDistPowB_dTred(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);
            d2DeltaPowB_ddelta_dtau = ddDistPowB_dRhoRed_dTred(T_red, rho_red, A_[i-39], B_[i-39], a_[i-39], b_[i-39], betap_[i-39]);

            // Psi function terms
            Psi = PsiFunc_(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_ddelta = dPsi_dRhoRed(T_red, rho_red, C_[i-39], D_[i-39]);
            dPsi_dtau = dPsi_dTred(T_red, rho_red, C_[i-39], D_[i-39]);
            d2Psi_ddelta_dtau = ddPsi_dRhoRed_dTred(T_red, rho_red, C_[i-39], D_[i-39]);

            // Sum
            s4 += n_[i] * (
                pow(Delta, b_[i-39]) * (dPsi_dtau + rho_red * d2Psi_ddelta_dtau) + 
                rho_red * dDeltaPowB_ddelta * dPsi_dtau +
                dDeltaPowB_dtau * (Psi + rho_red * dPsi_ddelta) + 
                d2DeltaPowB_ddelta_dtau * rho_red * Psi
            );
        }

        // Return total sum
        Evaluation s = s1 + s2 + s3 + s4;
        return s;
    }

    // ///////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS FOR RESIDUAL HELMHOLTZ AND ITS DERIVATIVES (TABLE 32)
    // ///////////////////////////////////////////////////////////////////////
    template <class Evaluation>
    static Evaluation thetaFunc_(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, const Scalar beta)
    /*!
    * \brief Helper function for the residual part of Helmholtz energy. Eq. (6.5) in Span & Wagner (1996).
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    {
        return (1 - T_red) + A * pow((rho_red - 1) * (rho_red - 1), 1 / (2 * beta));
    }
    
    template <class Evaluation>
    static Evaluation distFunc_(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, const Scalar B, 
        const Scalar a, const Scalar beta)
    /*!
    * \brief Helper function for the residual part of Helmholtz energy. Eq. (6.5) in Span & Wagner (1996).
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    {
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        return theta * theta + B * pow((rho_red - 1) * (rho_red - 1), a);
    }

    template <class Evaluation>
    static Evaluation dDist_dRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, const Scalar B, 
        const Scalar a, const Scalar beta)
    /*!
    * \brief Dervivative of distance function wrt to reduced density
    */
    {
        // Intermediate calcuations
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        Evaluation rho_red_term_squared = (rho_red - 1) * (rho_red - 1);

        Evaluation interm_1 = A * theta * (2 / beta) * pow(rho_red_term_squared, (1 / (2 * beta)) - 1);
        Evaluation interm_2 = 2 * B * a * pow(rho_red_term_squared, (a - 1));

        return (rho_red - 1) * (interm_1 + interm_2);
    }
    
    template <class Evaluation>
    static Evaluation ddDist_ddRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar beta)      
    /*!
    * \brief Second dervivative of distance function wrt to reduced density
    */
    {
        // Intermediate calcuations
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        Evaluation dDelta_ddelta = dDist_dRhoRed(T_red, rho_red, A, B, a, beta);
        Evaluation rho_red_term_squared = (rho_red - 1) * (rho_red - 1);

        Evaluation interm_1 = 4 * B * a * (a - 1) * pow(rho_red_term_squared, a - 2);
        Evaluation interm_2 = 2 * A * A * (1 / beta) * (1 / beta) * pow(rho_red_term_squared, (1 / beta) - 2);
        Evaluation interm_3 = A * theta * (4 / beta) * (1 / (2 * beta) - 1) * pow(rho_red_term_squared, (1 / (2 * beta)) - 2);

        return (1 / (rho_red - 1)) * dDelta_ddelta + rho_red_term_squared * (interm_1 + interm_2 + interm_3);
    }
    
    template <class Evaluation>
    static Evaluation dDistPowB_dRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar b, const Scalar beta)
    /*!
    * \brief Dervivative of distance function to the power of b wrt to reduced density
    */
    {
        // Intermediate calculations
        Evaluation Delta = distFunc_(T_red, rho_red, A, B, a, beta);
        Evaluation dDelta_ddelta = dDist_dRhoRed(T_red, rho_red, A, B, a, beta);

        return b * pow(Delta, b - 1) * dDelta_ddelta;
    }

    template <class Evaluation>
    static Evaluation ddDistPowB_ddRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar b, const Scalar beta)
    /*!
    * \brief Second dervivative of distance function to the power of b wrt to reduced density
    */
    {
        // Intermediate calculations
        Evaluation Delta = distFunc_(T_red, rho_red, A, B, a, beta);
        Evaluation dDelta_ddelta = dDist_dRhoRed(T_red, rho_red, A, B, a, beta);
        Evaluation d2Delta_d2delta = ddDist_ddRhoRed(T_red, rho_red, A, B, a, beta);

        return b * (pow(Delta, b - 1) * d2Delta_d2delta + (b - 1) * pow(Delta, b - 2) * dDelta_ddelta * dDelta_ddelta);
    }

    template <class Evaluation>
    static Evaluation dDistPowB_dTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar b, const Scalar beta)
    /*!
    * \brief Dervivative of distance function to the power of b wrt to reciprocal reduced temperature
    */
    {
        // Intermediate calculation
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        Evaluation Delta = distFunc_(T_red, rho_red, A, B, a, beta);

        return -2 * theta * b * pow(Delta, b - 1);
    }

    template <class Evaluation>
    static Evaluation ddDistPowB_ddTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar b, const Scalar beta)
    /*!
    * \brief Second dervivative of distance function to the power of b wrt to reciprocal reduced temperature
    */
    {
        // Intermediate calculation
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        Evaluation Delta = distFunc_(T_red, rho_red, A, B, a, beta);

        return 2 * b * pow(Delta, b - 1) + 4 * theta * theta * b * (b - 1) * pow(Delta, b - 2);
    }

    template <class Evaluation>
    static Evaluation ddDistPowB_dRhoRed_dTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar A, 
        const Scalar B, const Scalar a, const Scalar b, const Scalar beta)
    /*!
    * \brief Second dervivative of distance function to the power of b wrt to reciprocal reduced temperature 
    * and reduced density
    */
    {
        // Intermediate calculation
        Evaluation theta = thetaFunc_(T_red, rho_red, A, beta);
        Evaluation Delta = distFunc_(T_red, rho_red, A, B, a, beta);
        Evaluation dDelta_ddelta = dDist_dRhoRed(T_red, rho_red, A, B, a, beta);
        Evaluation rho_red_term_squared = (rho_red - 1) * (rho_red - 1);

        Evaluation interm_1 = A * b * (2 / beta) * pow(Delta, b - 1) * (rho_red - 1) * pow(rho_red_term_squared, (1 / ( 2 *beta)) - 1);
        Evaluation interm_2 = 2 * theta * b * (b - 1) * pow(Delta, b - 2) * dDelta_ddelta;

        return -interm_1 - interm_2;
    }


    template <class Evaluation>
    static Evaluation PsiFunc_(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, const Scalar D)
    /*!
    * \brief Helper function for the residual part of Helmholtz energy. Eq. (6.5) in Span & Wagner (1996).
    * 
    * \param T_red reciprocal reduced temperature [-]
    * \param rho_red reduced density [-]
    */
    {
        return exp(-C * (rho_red - 1) * (rho_red - 1) - D * (T_red - 1) * (T_red - 1));
    }

    template <class Evaluation>
    static Evaluation dPsi_dRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, const Scalar D)
    /*!
    * \brief Derivative of Psi function wrt reduced density
    */
    {
        // Intermediate calculation
        Evaluation Psi = PsiFunc_(T_red, rho_red, C, D);

        return -2 * C * (rho_red - 1) * Psi;
    }

    template <class Evaluation>
    static Evaluation ddPsi_ddRhoRed(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, const Scalar D)
    /*!
    * \brief Second derivative of Psi function wrt reduced density
    */
    {
        // Intermediate calculation
        Evaluation Psi = PsiFunc_(T_red, rho_red, C, D);

        return 2 * C * Psi * (2 * C * (rho_red - 1) * (rho_red - 1) - 1);
    }

    template <class Evaluation>
    static Evaluation dPsi_dTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, const Scalar D)
    /*!
    * \brief Derivative of Psi function wrt reciprocal reduced temperature
    */
    {
        // Intermediate calculation
        Evaluation Psi = PsiFunc_(T_red, rho_red, C, D);

        return -2 * D * (T_red - 1) * Psi;
    }

    template <class Evaluation>
    static Evaluation ddPsi_ddTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, const Scalar D)
    /*!
    * \brief Second derivative of Psi function wrt reciprocal reduced temperature
    */
    {
        // Intermediate calculation
        Evaluation Psi = PsiFunc_(T_red, rho_red, C, D);

        return 2 * D * Psi * (2 * D * (T_red - 1) * (T_red - 1) - 1);
    }

    template <class Evaluation>
    static Evaluation ddPsi_dRhoRed_dTred(const Evaluation& T_red, const Evaluation& rho_red, const Scalar C, 
        const Scalar D)
    /*!
    * \brief Second derivative of Psi function wrt reduced density and reciprocal reduced temperature
    */
    {
        // Intermediate calculation
        Evaluation Psi = PsiFunc_(T_red, rho_red, C, D);

        return 4 * C * D * (rho_red - 1) * (T_red - 1) * Psi;
    }

    template <class Evaluation>
    static Evaluation expTerm_(const Evaluation& T_red, const Evaluation& rho_red, const Scalar alpha, 
        const Scalar beta, const Scalar gamma, const Scalar epsilon)
    {
        return exp(-alpha * (rho_red - epsilon) * (rho_red - epsilon) - beta * (T_red - gamma) * (T_red - gamma));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*!
    * \brief Objective function for converting pressure to reduced density in root-finding algorithm.
    * 
    * \param rho_red reduced density [-]
    * \param pg gas phase pressure [Pa]
    * \param temperature temperature [K]
    */
    template <class Evaluation> 
    static Evaluation rootFindingObj_(const Evaluation& rho_red, const Evaluation& temperature, const Evaluation& pg)
    {
        // Temporary calculations
        Evaluation T_red = criticalTemperature() / temperature;  // reciprocal reduced temp.
        Evaluation p_MPa = pg / 1.0e6;  // Pa --> MPa
        Scalar R = IdealGas::R / (molarMass() * 1e3);  // J/(mol*K) --> kJ/(kg*K)
        Evaluation rho_cRT = criticalDensity() * R * temperature * 1e-3;  // MPa

        // Eq. (56) in Span et al. (2000)
        Evaluation dphir_dRedRho = derivResHelmholtzWrtRedRho(T_red, rho_red);
        Evaluation obj = rho_red * rho_cRT * (1 + rho_red * dphir_dRedRho) - p_MPa;
        return obj;
    }

    /*!
    * Derivative of vapor pressure wrt temperature. Needed in saturationTemperature().
    */
    template <class Evaluation>
    static Evaluation derVaporPressure(const Evaluation& T)
    {
        static constexpr Scalar a[4] =
            { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
        static constexpr Scalar t[4] =
            { 1.0, 1.5, 2.0, 4.0 };

        // Exponents needed for vapor pressure (Eq. 3.13 on p. 1524) and its derivative
        Evaluation exponent = 0;
        Evaluation exponent_der = 0;
        Scalar T_c = criticalTemperature();
        Evaluation Tred = T / T_c;
        for (int i = 0; i < 4; ++i) {
            exponent += a[i] * pow(1 - Tred, t[i]);
            exponent_der += a[i] * t[i] * pow(1 - Tred, t[i] - 1);
        }
        exponent *= 1.0/Tred;
        exponent_der *= 1.0/Tred;

        // Vapor pressure (could call vaporPressure() also, and delete exponent calculation above...)
        Evaluation ps = exp(exponent)*criticalPressure();

        return (-ps / T_c) * (exponent_der + exponent / Tred);
    }

    /*!
    * \brief Saturation temperature
    */
    template <class Evaluation>
    static Evaluation saturationTemperature(const Evaluation& p)
    {
        // Some properties for simple referencing
        Scalar p_c = criticalPressure();
        Scalar T_c = criticalTemperature();
        Scalar p_t = triplePressure();
        Scalar T_t = tripleTemperature();

        // Newton parameters
        int max_iter = 60;
        Scalar tol = 1e-7;

        // Initial guess
        Evaluation Ts = T_c + (p - p_c) * (T_t - T_c) / (p_t - p_c);

        // Newton loop
        Evaluation ps;
        Evaluation dps_dT;
        Evaluation dTs;
        for (int iter = 1; iter <= max_iter; ++iter) {
            // Vapor pressure and its derivative
            ps = vaporPressure(Ts);
            dps_dT = derVaporPressure(Ts);

            // Newton step
            dTs = (p - ps) / dps_dT;
            Ts += dTs;

            // Check tolerance
            if (abs(dTs) < tol) {
                return Ts;
            }
        }
        throw std::runtime_error("No saturation temperature could be found using Newton!");
    }

    template <class Evaluation>
    static Evaluation saturatedLiquidDensity(const Evaluation& T) {
        // Reduced temperature
        Evaluation Tred = T / criticalTemperature();

        // Parameters for calculation
        static constexpr Scalar a[4] = {1.9245108, -0.62385555, -0.32731127, 0.39245142};
        static constexpr Scalar t[4] = {0.34, 0.5, 1.6666666666666667, 1.8333333333333333};

        // Calculate vapor density in similar way as in vaporPressure
        Evaluation exponent = 0;
        for (int i = 0; i < 4; ++i) {
            exponent += a[i] * pow(1 - Tred, t[i]);
        }
        return criticalDensity() * exp(exponent);
    }

    template <class Evaluation>
    static Evaluation saturatedVaporDensity(const Evaluation& T) {
        // Reduced temperature
        Evaluation Tred = T / criticalTemperature();

        // Parameters for calculation (different for liquid and vapor)
        static constexpr Scalar a[5] = {-1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252};
        static constexpr Scalar t[5] = {0.340, 0.5, 1., 2.3333333333333335, 4.666666666666667};

        // Calculate vapor density in similar way as in vaporPressure
        Evaluation exponent = 0;
        for (int i = 0; i < 5; ++i) {
            exponent += a[i] * pow(1 - Tred, t[i]);
        }
        return criticalDensity() * exp(exponent);
    }
    
    /*!
    * \brief Bracket for reduced density optimization based on phase region 
    */
    template <class Evaluation>
    static std::pair<Evaluation, Evaluation> reducedDensityBracket_(const Evaluation& T, const Evaluation& p)
    {
        // Some properties for simple referencing
        Scalar p_c = criticalPressure();
        Scalar T_c = criticalTemperature();
        Scalar rho_c = criticalDensity();

        // Declare upper and lower bracket
        Evaluation rho_red_min;
        Evaluation rho_red_max;

        // 
        // Sub-critical region
        //
        if (p <= p_c && T <= T_c) {
            // Find the saturation temperature
            Evaluation T_s = saturationTemperature(p);
            
            // Sub-critical, liquid region
            if (T <= T_s) {
                // Base lower bracket on saturated liquid density
                rho_red_min = saturatedLiquidDensity(T) / rho_c;

                // Upper bracket set high
                rho_red_max = 2000.0 / rho_c;
            }
            // Sub-critical, vapor region
            else {
                // Base lower estimate on saturated vapor density
                rho_red_min = saturatedVaporDensity(T) / rho_c;

                // Set upper bracket to a low number and adjust later (code below)
                rho_red_max = 1e-5;
            }
        }
        // 
        // Liquid region
        // 
        else if (p > p_c && T <= T_c) {
            // Set lower bracket to a reasonably low value and upper to some high value
            rho_red_min = 600.0 / rho_c;
            rho_red_max = 2000.0 / rho_c;
        }
        // 
        // Vapor region
        // 
        else if (p <= p_c && T > T_c) {
            // Set lower bracket to (effectively) zero and upper to some high value
            rho_red_min = 1e-5;
            rho_red_max = 2000.0 / rho_c;
        }
        // 
        // Supercritical
        // 
        else {
            // Set lower bracket to a reasonably low value and upper to some high value
            rho_red_min = 25.0 / rho_c;
            rho_red_max = 2000.0 / rho_c;
        }

        // 
        // Finalize bracket 
        //
        // Declarations
        double factor = 1.1;
        int max_iter = 100;

        // Initial objective function calls
        Evaluation fmin = rootFindingObj_(rho_red_min, T, p);
        Evaluation fmax = rootFindingObj_(rho_red_max, T, p);

        // Loop until fmin and fmax have opposite signs
        for (int i = 0; i < max_iter; ++i) {
            // Check if sign is opposite and return bracket if true
            if (fmin * fmax < 0.0) {
                // Check order before returning
                if (rho_red_min > rho_red_max) {
                    Evaluation tmp = rho_red_max;
                    rho_red_max = rho_red_min;
                    rho_red_min = tmp;
                }
                return {rho_red_min, rho_red_max};
            }
            
            // Check fmin and fmax, and adjust bracket
            if (abs(fmin) < abs(fmax)) {
                rho_red_min += factor * (rho_red_min - rho_red_max);
                fmin = rootFindingObj_(rho_red_min, T, p);
            }
            else {
                rho_red_max += factor * (rho_red_max - rho_red_min);
                fmax = rootFindingObj_(rho_red_max, T, p);
            }
        }
        throw std::runtime_error("Failed to generate a bracket for reduced density calculation!");
    }
};

}

#endif