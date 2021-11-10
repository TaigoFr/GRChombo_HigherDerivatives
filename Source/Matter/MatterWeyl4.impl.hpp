/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERWEYL4_HPP_)
#error "This file should only be included through MatterWeyl4.hpp"
#endif

#ifndef MATTERWEYL4_IMPL_HPP_
#define MATTERWEYL4_IMPL_HPP_

template <class matter_t>
template <class data_t>
void MatterWeyl4<matter_t>::compute(Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
                                                    "MatterWeyl4::compute");

    gq.set_formulation(m_formulation,
                       CCZ4::params_t() /*params don't matter here*/);

    const auto emtensor = m_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, m_G_Newton);

    // Compute the E and B fields
    EBFields_t<data_t> ebfields;
    ebfields.E = gq.get_weyl_electric_part();
    ebfields.B = gq.get_weyl_magnetic_part();

    // work out the Newman Penrose scalar
    const Coordinates<data_t> coords(current_cell, m_dx, m_center);
    NPScalar_t<data_t> out =
        compute_Weyl4(ebfields, vars, d1, d2, gq.get_h_UU(), coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_Weyl4_Re);
    current_cell.store_vars(out.Im, c_Weyl4_Im);
}

#endif /* MATTERWEYL4_IMPL_HPP_ */
