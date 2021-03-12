/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEAKFIELDCONDITIONDIAGNOSTIC
#define WEAKFIELDCONDITIONDIAGNOSTIC

#include "C2EFT.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

template <class System> class WeakFieldConditionDiagnostic
{
    // Use the variable definitions in MatterCCZ4
    template <class data_t>
    using Vars = typename MatterCCZ4<C2EFT<System>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4<C2EFT<System>>::template Diff2Vars<data_t>;

  public:
    WeakFieldConditionDiagnostic(const C2EFT<System> &a_matter, double m_dx,
                                 int a_formulation,
                                 const CCZ4::params_t &a_ccz4_params)
        : my_matter(a_matter), m_formulation(a_formulation),
          m_ccz4_params(a_ccz4_params), m_deriv(m_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec =
            m_deriv.template advection<Vars>(current_cell, vars.shift);

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
        gq.set_advection(advec);
        gq.set_formulation(m_formulation, m_ccz4_params);

        const auto emtensor = my_matter.compute_emtensor(gq);
        data_t weak_field = my_matter.weak_field_var(emtensor, gq);
        data_t weak_field_condition =
            my_matter.weak_field_condition(weak_field, gq);

        data_t weak_field_after_WFC = weak_field * (1. - weak_field_condition);

        current_cell.store_vars(weak_field, c_WeakFieldVar);
        current_cell.store_vars(weak_field_condition, c_WeakFieldCondition);
        current_cell.store_vars(weak_field_after_WFC, c_WeakFieldVar_after_WFC);
    }

  protected:
    const C2EFT<System> &my_matter; //!< The matter object, e.g. a scalar field.
    int m_formulation;
    const CCZ4::params_t &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
};

#endif /* WEAKFIELDCONDITIONDIAGNOSTIC */
