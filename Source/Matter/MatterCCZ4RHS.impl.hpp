/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCCZ4RHS_HPP_)
#error "This file should only be included through MatterCCZ4RHS.hpp"
#endif

#ifndef MATTERCCZ4RHS_IMPL_HPP_
#define MATTERCCZ4RHS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t, class gauge_t, class deriv_t>
MatterCCZ4RHS<matter_t, gauge_t, deriv_t>::MatterCCZ4RHS(
    matter_t a_matter, CCZ4_params_t<typename gauge_t::params_t> a_params,
    double a_dx, double a_sigma, int a_formulation, double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_G_Newton(a_G_Newton)
{
}

template <class matter_t, class gauge_t, class deriv_t>
template <class data_t>
void MatterCCZ4RHS<matter_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{
    CH_TIME("MatterCCZ4::compute");

    // copy data from chombo gridpoint into local variables
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        this->m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

    GeometricQuantities<data_t, Vars, Diff2Vars, gauge_t> gq(
        matter_vars, d1, d2, "MatterCCZ4RHS::compute");

    gq.set_advection(advec);
    gq.set_gauge(this->m_gauge);
    gq.set_formulation(this->m_formulation, this->m_params);
    gq.set_cosmological_constant(this->m_cosmological_constant);

    // add evolution of matter fields themselves
    // done before to avoid getting contributions from the EM-tensor
    Vars<data_t> matter_rhs;
    my_matter.add_matter_rhs(matter_rhs, gq);

    const auto emtensor = my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, m_G_Newton);

    gq.compute_rhs_equations(matter_rhs);

    // Add dissipation to all terms
    this->m_deriv.add_dissipation(matter_rhs, current_cell, this->m_sigma);

    // Write the matter_rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

#endif /* MATTERCCZ4_IMPL_HPP_ */
