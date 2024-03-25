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

// #include "MinkowskiPerturbed.hpp"
// #include "SchwarzschildIsotropic.hpp"
// #include "SchwarzschildKS.hpp"
// #include "KerrBH.hpp"
// #include "Schwarzschild_SolvedConstraints.hpp"
#include "BoostedSchwarzschild_SolvedConstraints.hpp"
#include "SingleBH.hpp"

// which one to use:
// typedef MinkowskiPerturbed InitialData;
// typedef SchwarzschildIsotropic InitialData;
// typedef SchwarzschildKS InitialData;
// typedef KerrBH InitialData;
// typedef Schwarzschild_SolvedConstraints InitialData;
// typedef BoostedSchwarzschild_SolvedConstraints InitialData;

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
        pp.load("use_initial_data_with_solved_constraints",
                use_initial_data_with_solved_constraints);

        // Initial data
        pp.load("mass", id_params_with_constraints.mass);
        id_params_no_constraints.mass = id_params_with_constraints.mass;
        // pp.load("spin", id_params.spin); // for Kerr only
        // pp.load("amplitude", id_params.amplitude);
        // pp.load("r0", id_params.r0);
        id_params_with_constraints.center = center;
        id_params_no_constraints.center = center;

        if (use_initial_data_with_solved_constraints)
        {
            pp.load("boost_velocity", id_params_with_constraints.boost_velocity,
                    {0.});
            pout() << "Using boost_velocity = ("
                   << id_params_with_constraints.boost_velocity[0] << ", "
                   << id_params_with_constraints.boost_velocity[1] << ", "
                   << id_params_with_constraints.boost_velocity[2] << ")"
                   << std::endl;
        }
        else
        {
            pp.load("momentum", id_params_no_constraints.momentum, {0.});
            pout() << "Using momentum = ("
                   << id_params_no_constraints.momentum[0] << ", "
                   << id_params_no_constraints.momentum[1] << ", "
                   << id_params_no_constraints.momentum[2] << ")" << std::endl;
        }

        // pp.load("epsilon", hd_params.epsilon);
        // pout() << "Using epsilon = " << hd_params.epsilon << std::endl;

        // not called 'epsilon_final' to leave it compatible with Ramiro's ID
        pp.load("epsilon", hd_params.epsilon_final);
        pp.load("time_epsilon_start", hd_params.time_epsilon_start);
        pp.load("time_epsilon_end", hd_params.time_epsilon_end);
        pout() << "Using epsilon_final = " << hd_params.epsilon_final
               << std::endl;
        CH_assert(hd_params.time_epsilon_start >= 0. &&
                  hd_params.time_epsilon_start <= hd_params.time_epsilon_end);

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
            hd_params.chi_width = (1. - hd_params.chi_threshold_percentage) /
                                  (4. * hd_params.chi_threshold_percentage);
        }
        double time = 0., spin = 0.;
        // hd_params.update_min_chi(time, id_params.spin);
        hd_params.update_min_chi(time, spin);

        // this is such that the  'epsilon' in the EOM is replaced by
        // 'epsilon' when doing 'kappa / 2 * EM-tensor'
        G_Newton = 1.;
        // hd_params.epsilon /= (G_Newton * 8. * M_PI);
        hd_params.epsilon_final /= (G_Newton * 8. * M_PI);
        hd_params.epsilon =
            hd_params.time_epsilon_end == 0. ? hd_params.epsilon_final : 0.;

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

        pp.load("use_sigma_radial_decay", system_params.use_sigma_radial_decay);
        pout() << "Using use_sigma_radial_decay = "
               << system_params.use_sigma_radial_decay << std::endl;
        if (system_params.use_sigma_radial_decay)
        {
            pp.load("sigma_asymptotic", system_params.sigma_asymptotic);
            pout() << "Using sigma_asymptotic = "
                   << system_params.sigma_asymptotic << std::endl;
            pp.load("sigma_decay_length", system_params.sigma_decay_length);
            pout() << "Using sigma_decay_length = "
                   << system_params.sigma_decay_length << std::endl;
        }

        pp.load("use_tau_chi_decay", system_params.use_tau_chi_decay);
        pout() << "Using use_tau_chi_decay = "
               << system_params.use_tau_chi_decay << std::endl;

        if (system_params.use_tau_chi_decay)
        {
            pp.load("tau_asymptotic", system_params.tau_asymptotic);
            pout() << "Using tau_asymptotic = " << system_params.tau_asymptotic
                   << std::endl;
            pp.load("tau_decay_length", system_params.tau_decay_length);
            pout() << "Using tau_decay_length = "
                   << system_params.tau_decay_length << std::endl;
            pp.load("tau_decay_width", system_params.tau_decay_width);
            pout() << "Using tau_decay_width = "
                   << system_params.tau_decay_width << std::endl;
        }

        pp.load("use_sigma_chi_decay", system_params.use_sigma_chi_decay);
        pout() << "Using use_sigma_chi_decay = "
               << system_params.use_sigma_chi_decay << std::endl;
        if (system_params.use_sigma_chi_decay)
        {
            pp.load("sigma_asymptotic", system_params.sigma_asymptotic);
            pout() << "Using sigma_asymptotic = "
                   << system_params.sigma_asymptotic << std::endl;
            pp.load("sigma_decay_length", system_params.sigma_decay_length);
            pout() << "Using sigma_decay_length = "
                   << system_params.sigma_decay_length << std::endl;
            pp.load("sigma_decay_width", system_params.sigma_decay_width);
            pout() << "Using sigma_decay_width = "
                   << system_params.sigma_decay_width << std::endl;
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
        CH_assert(system_params.version >= 1 && system_params.version <= 4);

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

        if (system_params.version == 2 || system_params.version == 3 || system_params.version == 4)
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
        pp.load("AH_initial_guess", AH_initial_guess,
                0.5 * id_params_with_constraints.mass);
#endif
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
        Schwarzschild_SolvedConstraints id(id_params_with_constraints,
                                           1. /*dummy*/);
        SingleBH id2(id_params_no_constraints, 1. /*dummy*/);
    }

    double G_Newton;
    bool activate_extraction;

    BoostedSchwarzschild_SolvedConstraints::params_t id_params_with_constraints;
    BoostedBH::params_t id_params_no_constraints;

    C2EFT<System>::params_t hd_params;
    System::params_t system_params;

    diffusion_params_t diffusion_params;

    bool use_initial_data_with_solved_constraints;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
