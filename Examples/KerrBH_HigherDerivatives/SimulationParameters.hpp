/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "C2EFT.hpp"
#include "CCZ4.hpp"

#ifdef USE_EBSYSTEM
#include "EBSystem.hpp"
typedef EBSystem System;
#elif USE_CSYSTEM
#include "CSystem.hpp"
typedef CSystem System;
#else
#error "Please define either USE_CSYSTEM or USE_EBSYSTEM"
#endif

// #include "MinkowskiPerturbed.hpp"
// #include "SchwarzschildIsotropic.hpp"
// #include "SchwarzschildKS.hpp"
#include "KerrBH.hpp"

// which one to use:
// typedef MinkowskiPerturbed InitialData;
// typedef SchwarzschildIsotropic InitialData;
// typedef SchwarzschildKS InitialData;
typedef KerrBH InitialData;

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("mass", id_params.mass);
        pp.load("spin", id_params.spin); // for Kerr only
        // pp.load("amplitude", id_params.amplitude);
        // pp.load("r0", id_params.r0);
        id_params.center = center;

        pp.load("epsilon", hd_params.epsilon);
        // pp.load("chi_threshold", hd_params.chi_threshold); // automatic now
        pp.load("weak_field_threshold", hd_params.weak_field_threshold);
        pp.load("weak_field_width", hd_params.weak_field_width);
        // pp.load("chi_ignore_threshold", hd_params.chi_ignore_threshold); //
        // automatic now

        pp.load("chi_damp_coeff", hd_params.chi_damp_coeff);
        pp.load("chi_damp_timescale", hd_params.chi_damp_timescale);
        pp.load("chi_threshold_percentage", hd_params.chi_threshold_percentage);
        CH_assert(hd_params.chi_threshold_percentage <
                  0.98); // just to make sure I don't do anything stupid

        if (pp.contains("chi_width"))
            pp.load("chi_width", hd_params.chi_width);
        else
        {
            // ensure excision is O(10^-4) just before the horizon
            hd_params.chi_width =
                (1. - hd_params.chi_threshold_percentage) / 4.;
        }
        hd_params.update_min_chi(0., id_params.spin);

        // this is such that the  'epsilon' in the EOM is replaced by
        // 'epsilon' when doing 'kappa / 2 * EM-tensor'
        G_Newton = 1.;
        hd_params.epsilon /= (G_Newton * 8. * M_PI);

        pp.load("tau", system_params.tau);
#ifdef USE_CSYSTEM
        pp.load("c_sigma", system_params.sigma);
#endif

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * id_params.mass);
#endif
    }

    double G_Newton;

    // Schwarzschild bh initial data
    InitialData::params_t id_params;

    C2EFT<System>::params_t hd_params;
    System::params_t system_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
