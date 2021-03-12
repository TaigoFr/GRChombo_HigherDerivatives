/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4_HPP_)
#error "This file should only be included through CCZ4.hpp"
#endif

#ifndef CCZ4_IMPL_HPP_
#define CCZ4_IMPL_HPP_

#define COVARIANTZ4
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

#include "GeometricQuantities.hpp"

inline CCZ4::CCZ4(params_t params, double dx, double sigma, int formulation,
                  double cosmological_constant)
    : m_params(params), m_sigma(sigma), m_formulation(formulation),
      m_cosmological_constant(cosmological_constant), m_deriv(dx)
{
    // A user who wants to use BSSN should also have damping paramters = 0
    if (m_formulation == USE_BSSN)
    {
        if ((m_params.kappa1 != 0.) || (params.kappa2 != 0.) ||
            (params.kappa3 != 0.))
        {
            MayDay::Error("BSSN formulation is selected - CCZ4 kappa values "
                          "should be set to zero in params");
        }
    }
    if (m_formulation > USE_BSSN)
        MayDay::Error("The requested formulation is not supported");
}

template <class data_t> void CCZ4::compute(Cell<data_t> current_cell) const
{
    CH_TIME("CCZ4::compute");

    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
    gq.set_advection(advec);
    gq.set_formulation(m_formulation, m_params);
    gq.set_cosmological_constant(m_cosmological_constant);

    Vars<data_t> rhs = gq.get_rhs_equations();

    m_deriv.add_dissipation(rhs, current_cell, m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

#endif /* CCZ4_IMPL_HPP_ */
