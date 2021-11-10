/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CHIRELAXATION_HPP_)
#error "This file should only be included through ChiRelaxation.hpp"
#endif

#ifndef CHIRELAXATION_IMPL_HPP_
#define CHIRELAXATION_IMPL_HPP_

template <class matter_t>
ChiRelaxation<matter_t>::ChiRelaxation(const matter_t &a_matter, double dx,
                                       double relax_speed, double G_Newton)
    : my_matter(a_matter), m_relax_speed(relax_speed), m_G_Newton(G_Newton),
      m_deriv(dx)
{
}

template <class matter_t>
template <class data_t>
void ChiRelaxation<matter_t>::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variable and calculate derivs
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
                                                    "ChiRelaxation::compute");

    gq.set_advection(advec);
    // never needed so far, but should be included
    // gq.set_cosmological_constant(m_cosmological_constant);

    const auto emtensor = my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, m_G_Newton);

    // work out RHS including advection
    Vars<data_t> rhs;
    // All components that are not explicitly set in rhs_equation are 0
    VarsTools::assign(rhs, 0.);

    rhs.chi = m_relax_speed * gq.get_hamiltonian_constraint() / vars.chi;

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(rhs);
}

#endif /* CHIRELAXATION_IMPL_HPP_ */
