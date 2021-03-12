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
#include "MatterCCZ4RHSWithDiffusion.hpp"

#ifdef USE_EBSYSTEM
#include "EBSystem.hpp"
typedef EBSystem System;
#elif USE_CSYSTEM
#include "CSystem.hpp"
typedef CSystem System;
#else
#error "Please define either USE_CSYSTEM or USE_EBSYSTEM"
#endif

#include "BoostedBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> centerA, centerB;
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        FOR1(idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }

        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);
        pp.load("track_punctures", track_punctures, false);
        pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);

#ifdef USE_AHFINDER
        pp.load("AH_1_initial_guess", AH_1_initial_guess,
                0.5 * bh1_params.mass);
        pp.load("AH_2_initial_guess", AH_2_initial_guess,
                0.5 * bh2_params.mass);
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);
#endif

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
        double spin = 0.;
        hd_params.update_min_chi(0., spin);

        // this is such that the  'epsilon' in the EOM is replaced by
        // 'epsilon' when doing 'kappa / 2 * EM-tensor'
        G_Newton = 1.;
        hd_params.epsilon /= (G_Newton * 8. * M_PI);

        pp.load("tau", system_params.tau);
#ifdef USE_CSYSTEM
        pp.load("c_sigma", system_params.sigma);
#endif

        /////////////
        // Diffusion parameters
        pp.load("diffCFLFact", diffusion_params.diffCFLFact, 1e20);
        pp.load("lapidusCoeff", diffusion_params.lapidusCoeff, 0.001);
        pp.load("lapidusPower", diffusion_params.lapidusPower, 1.0);
        pp.load("diffusion_maximum", diffusion_params.diffusion_maximum, 1e20);
        diffusion_params.chiCutoff = hd_params.chi_threshold;
        diffusion_params.chiCutoff_width = hd_params.chi_width;
        /////////////
    }

    void check_params()
    {
        warn_parameter("massA", bh1_params.mass, bh1_params.mass >= 0,
                       "should be >= 0");
        warn_parameter("massB", bh2_params.mass, bh2_params.mass >= 0,
                       "should be >= 0");
        warn_array_parameter(
            "momentumA", bh1_params.momentum,
            std::sqrt(ArrayTools::norm2(bh1_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        warn_array_parameter(
            "momentumB", bh2_params.momentum,
            std::sqrt(ArrayTools::norm2(bh2_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        FOR1(idir)
        {
            std::string nameA = "centerA[" + std::to_string(idir) + "]";
            std::string nameB = "centerB[" + std::to_string(idir) + "]";
            double center_A_dir = bh1_params.center[idir];
            double center_B_dir = bh2_params.center[idir];
            warn_parameter(nameA, center_A_dir,
                           (center_A_dir >= 0.0) &&
                               (center_A_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
            warn_parameter(nameB, center_B_dir,
                           (center_B_dir >= 0.0) &&
                               (center_B_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
        }
        check_parameter("puncture_tracking_level", puncture_tracking_level,
                        (puncture_tracking_level >= 0) &&
                            (puncture_tracking_level <= max_level),
                        "must be between 0 and max_level (inclusive)");
    }

    double G_Newton;

    C2EFT<System>::params_t hd_params;
    System::params_t system_params;

    diffusion_params_t diffusion_params;

    // Initial data
    bool activate_extraction, track_punctures, calculate_constraint_norms;
    int puncture_tracking_level;

    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    double AH_2_initial_guess;
    bool AH_set_origins_to_punctures;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
