/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEEBDIFF
#define COMPUTEEBDIFF

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "SystemEB.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class ComputeEBdiff
{

    // Use the variable definitions in MatterCCZ4
    template <class data_t>
    using Vars = typename MatterCCZ4<C2EFT<SystemEB>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4<C2EFT<SystemEB>>::template Diff2Vars<data_t>;

  public:
    ComputeEBdiff(double m_dx, int a_formulation,
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

        data_t E_diff_squared = 0.;
        data_t B_diff_squared = 0.;
        FOR(i, j)
        {
            E_diff_squared +=
                (vars.Eij[i][j] - Eij[i][j]) * (vars.Eij[i][j] - Eij[i][j]);
            B_diff_squared +=
                (vars.Bij[i][j] - Bij[i][j]) * (vars.Bij[i][j] - Bij[i][j]);
        }

        current_cell.store_vars(sqrt(E_diff_squared), c_E_diff);
        current_cell.store_vars(sqrt(B_diff_squared), c_B_diff);
    }

  protected:
    int m_formulation;
    const CCZ4::params_t &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
};

#endif /* COMPUTEEBDIFF */
