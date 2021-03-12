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
#include "CCZ4.hpp"
#include "SchBH_KerrSchild.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp)
    {
        // Initial Sch data
        pp.load("sch_mass", sch_params.mass);
        // pp.load("sch_spin", sch_params.spin);
        pp.load("sch_center", sch_params.center, center);
    }

    SchBH::params_t sch_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
