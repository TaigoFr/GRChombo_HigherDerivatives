/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEC
#define COMPUTEC

#include "C2EFT.hpp"
#include "CSystem.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

/////////////////////////////////////////////////////
// In this system, this is only for the initial data!
// (unlike 'EBSystem' where E_phys and B_phys must be stored)
/////////////////////////////////////////////////////

class ComputeC
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<CSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<CSystem>>::template Diff2Vars<data_t>;

  public:
    ComputeC(double m_dx, int a_formulation,
             const CCZ4_params_t<> &a_ccz4_params)
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

        data_t C = gq.get_kretschmann();

        current_cell.store_vars(C, c_C);
    }

  protected:
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
};

#endif /* COMPUTEC */
