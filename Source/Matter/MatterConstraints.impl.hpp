/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCONSTRAINTS_HPP_)
#error "This file should only be included through MatterConstraints.hpp"
#endif

#ifndef MATTERCONSTRAINTS_IMPL_HPP_
#define MATTERCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterConstraints<matter_t>::MatterConstraints(const matter_t &a_matter,
                                               double dx, double G_Newton)
    : Constraints(dx, 0.0 /*No cosmological constant*/), my_matter(a_matter),
      m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void MatterConstraints<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<MatterMetricVars>();
    const auto d1 = m_deriv.template diff1<MatterMetricVars>(current_cell);
    const auto d2 = m_deriv.template diff2<MatterDiff2MetricVars>(current_cell);

    GeometricQuantities<data_t, MatterMetricVars, MatterDiff2MetricVars> gq(
        vars, d1, d2);
    gq.set_cosmological_constant(m_cosmological_constant);

    const auto emtensor = my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, m_G_Newton);

    Vars<data_t> out;
    out.Ham = gq.get_hamiltonian_constraint();
    out.Mom = gq.get_momentum_constraints();

    // Write the constraints into the output FArrayBox
    current_cell.store_vars(out);
}

#endif /* MATTERCONSTRAINTS_IMPL_HPP_ */
