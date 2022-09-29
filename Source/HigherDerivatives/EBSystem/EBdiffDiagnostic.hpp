/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EBDIFFDIAGNOSTIC
#define EBDIFFDIAGNOSTIC

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "EBSystem.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ASSUMES ComputeEB WAS CALLED
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!

template <class gauge_t = MovingPunctureGauge> class EBdiffDiagnostic
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>, gauge_t>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>,
                               gauge_t>::template Diff2Vars<data_t>;

  public:
    EBdiffDiagnostic(
        double m_dx, int a_formulation,
        const CCZ4_params_t<typename gauge_t::params_t> &a_ccz4_params,
        bool use_last_index_raised)
        : m_deriv(m_dx), m_formulation(a_formulation),
          m_ccz4_params(a_ccz4_params),
          m_use_last_index_raised(use_last_index_raised)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        GeometricQuantities<data_t, Vars, Diff2Vars, gauge_t> gq(
            vars, d1, d2, "EBdiffDiagnostic::compute");

        gq.set_formulation(m_formulation, m_ccz4_params);

        // get the LL E and B, but raise to LU if using the raised system
        Tensor<2, data_t> Eij = gq.get_weyl_electric_part();
        Tensor<2, data_t> Bij = gq.get_weyl_magnetic_part();

        if (m_use_last_index_raised)
        {
            const auto &h_UU = gq.get_metric_UU_spatial();
            Eij = TensorAlgebra::compute_dot_product(Eij, h_UU);
            Bij = TensorAlgebra::compute_dot_product(Bij, h_UU);
        }

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
    FourthOrderDerivatives m_deriv;
    int m_formulation;
    const CCZ4_params_t<typename gauge_t::params_t> &m_ccz4_params;
    bool m_use_last_index_raised;
};

#endif /* EBDIFFDIAGNOSTIC */
