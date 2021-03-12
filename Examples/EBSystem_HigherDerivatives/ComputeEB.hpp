/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEEB
#define COMPUTEEB

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "EBSystem.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeEB
{

    // Use the variable definitions in MatterCCZ4
    template <class data_t>
    using Vars = typename MatterCCZ4<C2EFT<EBSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4<C2EFT<EBSystem>>::template Diff2Vars<data_t>;

  public:
    ComputeEB(double m_dx, int a_formulation,
              const CCZ4::params_t &a_ccz4_params, const Interval &E_comps,
              const Interval &B_comps)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx), m_E_comps(E_comps), m_B_comps(B_comps)
    {
        CH_assert(m_E_comps.size() == (GR_SPACEDIM * (GR_SPACEDIM + 1) / 2));
        CH_assert(m_B_comps.size() == (GR_SPACEDIM * (GR_SPACEDIM + 1) / 2));
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

        current_cell.store_vars(Eij[0][0], m_E_comps.begin() + 0);
        current_cell.store_vars(Eij[0][1], m_E_comps.begin() + 1);
        current_cell.store_vars(Eij[0][2], m_E_comps.begin() + 2);
        current_cell.store_vars(Eij[1][1], m_E_comps.begin() + 3);
        current_cell.store_vars(Eij[1][2], m_E_comps.begin() + 4);
        current_cell.store_vars(Eij[2][2], m_E_comps.begin() + 5);

        current_cell.store_vars(Bij[0][0], m_B_comps.begin() + 0);
        current_cell.store_vars(Bij[0][1], m_B_comps.begin() + 1);
        current_cell.store_vars(Bij[0][2], m_B_comps.begin() + 2);
        current_cell.store_vars(Bij[1][1], m_B_comps.begin() + 3);
        current_cell.store_vars(Bij[1][2], m_B_comps.begin() + 4);
        current_cell.store_vars(Bij[2][2], m_B_comps.begin() + 5);
    }

  protected:
    int m_formulation;
    const CCZ4::params_t &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
    const Interval m_E_comps, m_B_comps;
};

#endif /* COMPUTEEB */
