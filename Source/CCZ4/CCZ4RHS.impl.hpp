/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4RHS_HPP_)
#error "This file should only be included through CCZ4RHS.hpp"
#endif

#ifndef CCZ4RHS_IMPL_HPP_
#define CCZ4RHS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

#include "GeometricQuantities.hpp"

template <class gauge_t, class deriv_t>
inline CCZ4RHS<gauge_t, deriv_t>::CCZ4RHS(
    CCZ4_params_t<typename gauge_t::params_t> a_params, double a_dx,
    double a_sigma, int a_formulation, double a_cosmological_constant)
    : m_params(a_params), m_gauge(a_params), m_sigma(a_sigma),
      m_formulation(a_formulation),
      m_cosmological_constant(a_cosmological_constant), m_deriv(a_dx)
{
    // A user who wants to use BSSN should also have damping paramters = 0
    if (m_formulation == USE_BSSN)
    {
        if ((m_params.kappa1 != 0.) || (m_params.kappa2 != 0.) ||
            (m_params.kappa3 != 0.))
        {
            MayDay::Error("BSSN formulation is selected - CCZ4 kappa values "
                          "should be set to zero in params");
        }
    }
    if (m_formulation > USE_BSSN)
        MayDay::Error("The requested formulation is not supported");
}

template <class gauge_t, class deriv_t>
template <class data_t>
void CCZ4RHS<gauge_t, deriv_t>::compute(Cell<data_t> current_cell) const
{
    CH_TIME("CCZ4::compute");

    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    GeometricQuantities<data_t, Vars, Diff2Vars, gauge_t> gq(vars, d1, d2);
    gq.set_advection_and_gauge(advec, m_gauge);
    gq.set_formulation(m_formulation, m_params);
    gq.set_cosmological_constant(m_cosmological_constant);

    Vars<data_t> rhs = gq.get_rhs_equations();

    m_deriv.add_dissipation(rhs, current_cell, m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

#endif /* CCZ4RHS_IMPL_HPP_ */
