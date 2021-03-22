/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CDIAGNOSTICS
#define CDIAGNOSTICS

#include "C2EFT.hpp"
#include "CSystem.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class CDiagnostics
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<CSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<CSystem>>::template Diff2Vars<data_t>;

  public:
    CDiagnostics(double m_dx, int a_formulation,
                 const CCZ4_params_t<> &a_ccz4_params, int a_C_comp = -1,
                 int a_C_diff_comp = -1)
        : m_formulation(a_formulation), m_ccz4_params(a_ccz4_params),
          m_deriv(m_dx), m_C_comp(a_C_comp), m_C_diff_comp(a_C_diff_comp)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
        gq.set_formulation(m_formulation, m_ccz4_params);

        data_t C = gq.get_kretschmann();

        if (m_C_diff_comp >= 0)
            current_cell.store_vars(abs(C - vars.C), m_C_diff_comp);
        if (m_C_comp >= 0)
            current_cell.store_vars(C, m_C_comp);
    }

  protected:
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;

    int m_C_comp, m_C_diff_comp;
};

#endif /* CDIAGNOSTICS */
