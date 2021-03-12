/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

inline Constraints::Constraints(double dx,
                                double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_cosmological_constant(cosmological_constant)
{
}

template <class data_t>
void Constraints::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<MetricVars>();
    const auto d1 = m_deriv.template diff1<MetricVars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2MetricVars>(current_cell);

    GeometricQuantities<data_t, MetricVars, Diff2MetricVars> gq(vars, d1, d2);
    gq.set_cosmological_constant(m_cosmological_constant);

    Vars<data_t> out;
    out.Ham = gq.get_hamiltonian_constraint();
    out.Mom = gq.get_momentum_constraints();

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out);
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
