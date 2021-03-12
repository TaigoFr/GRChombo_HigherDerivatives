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

#include "MinkowskiPerturbed.hpp"
// #include "SchwarzschildIsotropic.hpp"
// #include "SchwarzschildKS.hpp"

// which one to use:
typedef MinkowskiPerturbed InitialData;
// typedef SchwarzschildIsotropic InitialData;
// typedef SchwarzschildKS InitialData;

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
        // Relaxation params
        pp.load("relaxtime", relaxtime);
        pp.load("relaxspeed", relaxspeed);

        // Initial data
        // pp.load("mass", id_params.mass);
        pp.load("amplitude", id_params.amplitude);
        pp.load("r0", id_params.r0);
        id_params.center = center;

        pp.load("energy_scale", energy_scale);

        // this is such that the ''epsilon' in the EOM is replaced by
        // 'energy_scale' when doing 'kappa / 2 * EM-tensor'
        G_Newton = energy_scale / (8. * M_PI);

        pp.load("sigma", hd_params.sigma);
        pp.load("tau", hd_params.tau);
    }

    double energy_scale, G_Newton;

    // Schwarzschild bh initial data
    InitialData::params_t id_params;

    C2EFT::params_t hd_params;

    // Relaxation params
    Real relaxtime, relaxspeed;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
