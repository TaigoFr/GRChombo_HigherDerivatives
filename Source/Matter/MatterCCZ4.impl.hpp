/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCCZ4_HPP_)
#error "This file should only be included through MatterCCZ4.hpp"
#endif

#ifndef MATTERCCZ4_IMPL_HPP_
#define MATTERCCZ4_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterCCZ4<matter_t>::MatterCCZ4(const matter_t &a_matter, params_t params,
                                 double dx, double sigma, int formulation,
                                 double G_Newton)
    : CCZ4(params, dx, sigma, formulation, 0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void MatterCCZ4<matter_t>::compute(Cell<data_t> current_cell) const
{
    CH_TIME("MatterCCZ4::compute");

    // copy data from chombo gridpoint into local variables
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

    GeometricQuantities<data_t, Vars, Diff2Vars> gq(matter_vars, d1, d2);
    gq.set_advection(advec);
    gq.set_formulation(m_formulation, m_params);
    gq.set_cosmological_constant(m_cosmological_constant);

    const auto emtensor = my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, m_G_Newton);

    Vars<data_t> matter_rhs = gq.get_rhs_equations();

    // add evolution of matter fields themselves
    my_matter.add_matter_rhs(matter_rhs, gq);

    // Add dissipation to all terms
    m_deriv.add_dissipation(matter_rhs, current_cell, m_sigma);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

#endif /* MATTERCCZ4_IMPL_HPP_ */
