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

class EBdiffDiagnostic
{

    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<EBSystem>>::template Diff2Vars<data_t>;

  public:
    EBdiffDiagnostic(double m_dx) : m_deriv(m_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        data_t E_diff_squared = 0.;
        data_t B_diff_squared = 0.;
        FOR(i, j)
        {
            E_diff_squared += (vars.Eij[i][j] - vars.Ephys[i][j]) *
                              (vars.Eij[i][j] - vars.Ephys[i][j]);
            B_diff_squared += (vars.Bij[i][j] - vars.Bphys[i][j]) *
                              (vars.Bij[i][j] - vars.Bphys[i][j]);
        }

        current_cell.store_vars(sqrt(E_diff_squared), c_E_diff);
        current_cell.store_vars(sqrt(B_diff_squared), c_B_diff);
    }

  protected:
    FourthOrderDerivatives m_deriv;
};

#endif /* EBDIFFDIAGNOSTIC */
