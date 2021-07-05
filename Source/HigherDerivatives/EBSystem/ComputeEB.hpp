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
#include "MatterCCZ4RHS.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeEB
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Diff2Vars<data_t>;

  public:
    ComputeEB(double m_dx, int a_formulation,
              const CCZ4_params_t<> &a_ccz4_params, const Interval &E_comps,
              const Interval &B_comps, bool compute_time_derivatives = false)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx), m_E_comps(E_comps), m_B_comps(B_comps),
          m_compute_time_derivatives(compute_time_derivatives)
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

        Tensor<2, data_t> Eij = gq.get_weyl_electric_part();
        Tensor<2, data_t> Bij = gq.get_weyl_magnetic_part();

        if (m_compute_time_derivatives)
        {
            const auto advec =
                m_deriv.template advection<Vars>(current_cell, vars.shift);

            Tensor<2, data_t> dt_Eij = gq.compute_dt_weyl_electric_part(
                d1.Eij, d1.Bij, vars.Eij, vars.Bij, advec.Eij, advec.Bij);
            Tensor<2, data_t> dt_Bij = gq.compute_dt_weyl_magnetic_part(
                d1.Eij, d1.Bij, vars.Eij, vars.Bij, advec.Eij, advec.Bij);

            // replace tensors
            Eij = dt_Eij;
            Bij = dt_Bij;
        }

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
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
    const Interval m_E_comps, m_B_comps;
    bool m_compute_time_derivatives;
};

#endif /* COMPUTEEB */
