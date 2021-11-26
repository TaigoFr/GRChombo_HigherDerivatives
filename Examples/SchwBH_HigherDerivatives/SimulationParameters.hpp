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
#include "MovingPunctureGaugeEtaRadialDecay.hpp"

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
// #include "KerrBH.hpp"
#include "Schwarzschild_SolvedConstraints.hpp"

// which one to use:
// typedef MinkowskiPerturbed InitialData;
// typedef SchwarzschildIsotropic InitialData;
// typedef SchwarzschildKS InitialData;
// typedef KerrBH InitialData;
typedef Schwarzschild_SolvedConstraints InitialData;

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("mass", id_params.mass);
        // pp.load("spin", id_params.spin); // for Kerr only
        // pp.load("amplitude", id_params.amplitude);
        // pp.load("r0", id_params.r0);
        id_params.center = center;

        pp.load("epsilon", hd_params.epsilon);
        pout() << "Using epsilon = " << hd_params.epsilon << std::endl;

        // pp.load("chi_threshold", hd_params.chi_threshold); // automatic now
        pp.load("weak_field_threshold", hd_params.weak_field_threshold);
        pp.load("weak_field_width", hd_params.weak_field_width);
        // pp.load("chi_ignore_threshold", hd_params.chi_ignore_threshold); //
        // automatic now

        pp.load("chi_damp_coeff", hd_params.chi_damp_coeff);
        pp.load("chi_damp_timescale", hd_params.chi_damp_timescale);
        pp.load("chi_threshold_percentage", hd_params.chi_threshold_percentage);

        if (pp.contains("chi_width"))
            pp.load("chi_width", hd_params.chi_width);
        else
        {
            // ensure excision is O(10^-4) just before the horizon
            hd_params.chi_width =
                (1. - hd_params.chi_threshold_percentage) / 4.;
        }
        double time = 0., spin = 0.;
        // hd_params.update_min_chi(time, id_params.spin);
        hd_params.update_min_chi(time, spin);

        // this is such that the  'epsilon' in the EOM is replaced by
        // 'epsilon' when doing 'kappa / 2 * EM-tensor'
        G_Newton = 1.;
        hd_params.epsilon /= (G_Newton * 8. * M_PI);

        pp.load("tau", system_params.tau);
        pout() << "Using tau = " << system_params.tau << std::endl;

        pp.load("Box_transition", system_params.Box_transition);
        pout() << "Using Box_transition = " << system_params.Box_transition
               << std::endl;

        pp.load("use_tau_radial_decay", system_params.use_tau_radial_decay);
        pout() << "Using use_tau_radial_decay = "
               << system_params.use_tau_radial_decay << std::endl;
        if (system_params.use_tau_radial_decay)
        {
            pp.load("tau_asymptotic", system_params.tau_asymptotic);
            pout() << "Using tau_asymptotic = " << system_params.tau_asymptotic
                   << std::endl;
            pp.load("tau_decay_length", system_params.tau_decay_length);
            pout() << "Using tau_decay_length = "
                   << system_params.tau_decay_length << std::endl;
        }

#ifdef USE_CSYSTEM
        pp.load("c_sigma", system_params.sigma);
        pout() << "Using sigma = " << system_params.sigma << std::endl;

        pp.load("c_version", system_params.version);
        pout() << "Using C system version " << system_params.version
               << std::endl;
        CH_assert(system_params.version >= 1 && system_params.version <= 2);
        pp.load("advection_type", system_params.advection_type);
        pout() << "Using advection_type = " << system_params.advection_type
               << std::endl;

        if (system_params.version == 2)
        {
            pp.load("rescale_tau_by_lapse", system_params.rescale_tau_by_lapse);
            pp.load("rescale_sigma_by_lapse",
                    system_params.rescale_sigma_by_lapse);
            CH_assert(system_params.rescale_sigma_by_lapse >= 0 &&
                      system_params.rescale_sigma_by_lapse <= 2);

            pout() << "Using rescale_tau_by_lapse = "
                   << system_params.rescale_tau_by_lapse << std::endl;
            pout() << "Using rescale_sigma_by_lapse = "
                   << system_params.rescale_sigma_by_lapse << std::endl;

            if (system_params.advection_type == 1 ||
                system_params.advection_type == 2)
            {
                pp.load("advection_coeff", system_params.advection_coeff);
                pout() << "Using advection_coeff = "
                       << system_params.advection_coeff << std::endl;
            }
        }
#else
        pp.load("eb_version", system_params.version);
        pout() << "Using EB system version " << system_params.version
               << std::endl;
        CH_assert(system_params.version >= 1 && system_params.version <= 3);

        pp.load("rescale_tau_by_lapse", system_params.rescale_tau_by_lapse);
        pp.load("rescale_sigma_by_lapse", system_params.rescale_sigma_by_lapse);
        CH_assert(system_params.rescale_sigma_by_lapse >= 0 &&
                  system_params.rescale_sigma_by_lapse <= 2);
        pout() << "Using rescale_tau_by_lapse = "
               << system_params.rescale_tau_by_lapse << std::endl;
        pout() << "Using rescale_sigma_by_lapse = "
               << system_params.rescale_sigma_by_lapse << std::endl;

        pp.load("use_last_index_raised", system_params.use_last_index_raised);
        pout() << "Using use_last_index_raised = "
               << system_params.use_last_index_raised << std::endl;

        if (system_params.version == 2 || system_params.version == 3)
        {
            pp.load("advection_type", system_params.advection_type);
            pp.load("eb_sigma", system_params.sigma);
            pout() << "Using advection_type = " << system_params.advection_type
                   << std::endl;
            pout() << "Using sigma = " << system_params.sigma << std::endl;

            if (system_params.advection_type == 1 ||
                system_params.advection_type == 2)
            {
                pp.load("advection_coeff", system_params.advection_coeff);
                pout() << "Using advection_coeff = "
                       << system_params.advection_coeff << std::endl;
            }
        }
#endif

        /////////////
        // Diffusion parameters
        pp.load("diffCFLFact", diffusion_params.diffCFLFact, 1e20);
        pp.load("lapidusCoeff", diffusion_params.lapidusCoeff, 0.0);
        pp.load("lapidusPower", diffusion_params.lapidusPower, 1.0);
        diffusion_params.chi_damp_coeff = hd_params.chi_damp_coeff;
        diffusion_params.chi_damp_timescale = hd_params.chi_damp_timescale;
        pp.load("diffusion_chi_threshold_percentage",
                diffusion_params.chi_threshold_percentage,
                hd_params.chi_threshold_percentage);
        if (pp.contains("diffusion_chi_width"))
            pp.load("diffusion_chi_width", diffusion_params.chi_width);
        else
        {
            // ensure excision is O(10^-4) just before the horizon
            diffusion_params.chi_width =
                (1. - diffusion_params.chi_threshold_percentage) / 4.;
        }
        diffusion_params.update_min_chi(time, spin);
        /////////////

        pp.load("activate_extraction", activate_extraction, false);

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * id_params.mass);
#endif

        ccz4_params_modifiedGauge.kappa1 = ccz4_params.kappa1;
        ccz4_params_modifiedGauge.kappa2 = ccz4_params.kappa2;
        ccz4_params_modifiedGauge.kappa3 = ccz4_params.kappa3;
        ccz4_params_modifiedGauge.covariantZ4 = ccz4_params.covariantZ4;
        ccz4_params_modifiedGauge.lapse_advec_coeff =
            ccz4_params.lapse_advec_coeff;
        ccz4_params_modifiedGauge.lapse_power = ccz4_params.lapse_power;
        ccz4_params_modifiedGauge.lapse_coeff = ccz4_params.lapse_coeff;
        ccz4_params_modifiedGauge.shift_Gamma_coeff =
            ccz4_params.shift_Gamma_coeff;
        ccz4_params_modifiedGauge.shift_advec_coeff =
            ccz4_params.shift_advec_coeff;
        ccz4_params_modifiedGauge.eta = ccz4_params.eta;
        pp.load("eta_sigmoid_decay",
                ccz4_params_modifiedGauge.eta_sigmoid_decay, 17.);
        pout() << "Using eta_sigmoid_decay = "
               << ccz4_params_modifiedGauge.eta_sigmoid_decay << std::endl;
        pp.load("eta_sigmoid_chi_threshold",
                ccz4_params_modifiedGauge.eta_sigmoid_chi_threshold, 0.92);
        pout() << "Using eta_sigmoid_chi_threshold = "
               << ccz4_params_modifiedGauge.eta_sigmoid_chi_threshold
               << std::endl;
        pp.load("eta_asymptotic", ccz4_params_modifiedGauge.eta_asymptotic,
                0.1);
        pout() << "Using eta_asymptotic = "
               << ccz4_params_modifiedGauge.eta_asymptotic << std::endl;
    }

    void check_params()
    {
        // check_parameter("spin", id_params.spin,
        //                 abs(id_params.spin) <= id_params.mass,
        //                 "must be between -mass and +mass");

        check_parameter("chi_threshold_percentage",
                        hd_params.chi_threshold_percentage,
                        hd_params.chi_threshold_percentage < 0.98,
                        "must be sufficiently below 1");

        // try to create InitialData -> will give error if it has to
        InitialData id(id_params, 1. /*dummy*/);
    }

    double G_Newton;
    bool activate_extraction;

    // Schwarzschild bh initial data
    InitialData::params_t id_params;

    C2EFT<System>::params_t hd_params;
    System::params_t system_params;

    CCZ4_params_t<MovingPunctureGaugeEtaRadialDecay::params_t>
        ccz4_params_modifiedGauge;

    diffusion_params_t diffusion_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
