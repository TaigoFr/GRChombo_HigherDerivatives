/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEINITIALEB
#define COMPUTEINITIALEB

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "SystemEB.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeInitialEB
{

    // Use the variable definitions in MatterCCZ4
    template <class data_t>
    using Vars = typename MatterCCZ4<C2EFT<SystemEB>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4<C2EFT<SystemEB>>::template Diff2Vars<data_t>;

  public:
    ComputeInitialEB(double m_dx, int a_formulation,
                     const CCZ4::params_t &a_ccz4_params)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
        gq.set_formulation(m_formulation, m_ccz4_params);

        const auto &Eij = gq.get_weyl_electric_part();
        const auto &Bij = gq.get_weyl_magnetic_part();

        current_cell.store_vars(Eij[0][0], c_E11);
        current_cell.store_vars(Eij[0][1], c_E12);
        current_cell.store_vars(Eij[0][2], c_E13);
        current_cell.store_vars(Eij[1][1], c_E22);
        current_cell.store_vars(Eij[1][2], c_E23);
        current_cell.store_vars(Eij[2][2], c_E33);

        current_cell.store_vars(Bij[0][0], c_B11);
        current_cell.store_vars(Bij[0][1], c_B12);
        current_cell.store_vars(Bij[0][2], c_B13);
        current_cell.store_vars(Bij[1][1], c_B22);
        current_cell.store_vars(Bij[1][2], c_B23);
        current_cell.store_vars(Bij[2][2], c_B33);
    }

  protected:
    int m_formulation;
    const CCZ4::params_t &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
};

#endif /* COMPUTEINITIALEB */
