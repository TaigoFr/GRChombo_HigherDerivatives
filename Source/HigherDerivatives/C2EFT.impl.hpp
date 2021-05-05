/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(C2EFT_HPP_)
#error "This file should only be included through C2EFT.hpp"
#endif

#ifndef C2EFT_IMPL_HPP_
#define C2EFT_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class System>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
emtensor_t<data_t> C2EFT<System>::compute_emtensor(
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    emtensor_t<data_t> emtensor;

    Tensor<2, data_t, CH_SPACEDIM + 1> Tmn; // 4D

    compute_emtensor_4D(Tmn, gq);

    const auto &vars = gq.get_vars();
    const auto &n_U = gq.get_normal_U_ST();
    const auto &proj_LU = gq.get_projector_LU_ST();

    Tensor<1, data_t, CH_SPACEDIM + 1> Tmn_dot_normal_U =
        TensorAlgebra::compute_dot_product(Tmn, n_U);

    Tensor<1, data_t, CH_SPACEDIM + 1> Si_4D =
        TensorAlgebra::compute_dot_product(Tmn_dot_normal_U, proj_LU);

    Tensor<2, data_t, CH_SPACEDIM + 1> Sij_4D =
        TensorAlgebra::compute_dot_product(
            proj_LU, TensorAlgebra::compute_dot_product(Tmn, proj_LU), 1, 0);

    emtensor.rho = TensorAlgebra::compute_dot_product(Tmn_dot_normal_U, n_U);

    FOR(i)
    {
        emtensor.Si[i] = -Si_4D[i + 1];
        FOR(j) { emtensor.Sij[i][j] = Sij_4D[i + 1][j + 1]; }
    }

    // WEAK FIELD EXCISION
    // chi_ignore_threshold is just to make it faster, to avoid entering in
    // cells far from the BH
    if (m_apply_weak_field &&
        simd_compare_lt_any(vars.chi, m_params.chi_ignore_threshold))
    {
        data_t weak_field = weak_field_var(emtensor, gq);
        data_t weak_field_damp = 1. - weak_field_condition(weak_field, gq);
        emtensor.rho *= weak_field_damp;

        FOR(i)
        {
            emtensor.Si[i] *= weak_field_damp;
            FOR(j) { emtensor.Sij[i][j] *= weak_field_damp; }
        }
    }

    emtensor.S =
        TensorAlgebra::compute_trace(emtensor.Sij, gq.get_metric_UU_spatial());

    return emtensor;
}

template <class System>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void C2EFT<System>::compute_emtensor_4D(
    Tensor<2, data_t, CH_SPACEDIM + 1> &Tmn,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    const auto &g = gq.get_metric_ST();
    const auto &chris_ST = gq.get_chris_ST();

    data_t C;
    Tensor<1, data_t, CH_SPACEDIM + 1> d1_C_4D;
    Tensor<2, data_t, CH_SPACEDIM + 1> d2_C_4D;
    m_system.compute_C(C, d1_C_4D, d2_C_4D, gq);

    // compute Riemann
    // Can either be from evolved variables (if evolving
    // Eij or Bij for example), or otherwise simply use
    // gq.get_riemann_LLLU() and gq.get_riemann_LULU()
    Tensor<4, data_t, CH_SPACEDIM + 1> riemann_LLLU;
    Tensor<4, data_t, CH_SPACEDIM + 1> riemann_LULU;
    m_system.compute_Riemann(riemann_LLLU, riemann_LULU, gq);

    const Tensor<2, data_t, CH_SPACEDIM + 1> covd2_C =
        TensorAlgebra::covariant_derivative(d2_C_4D, d1_C_4D, chris_ST);

    FOR_ST(a, b)
    {
        Tmn[a][b] = -0.5 * C * C * g[a][b];
        FOR_ST(c, d)
        {
            Tmn[a][b] += 8. * riemann_LULU[a][c][b][d] * covd2_C[c][d];

            FOR_ST(e)
            {
                Tmn[a][b] += -4. * C * riemann_LULU[a][c][d][e] *
                             riemann_LLLU[b][c][e][d];
            }
        }
        Tmn[a][b] *= m_params.epsilon;
    }
}

// Adds in the RHS for the matter vars
template <class System>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t,
          template <typename> class rhs_vars_t>
void C2EFT<System>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    m_system.add_matter_rhs(total_rhs, gq);
}

template <class System>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
data_t C2EFT<System>::weak_field_var(
    const emtensor_t<data_t> &emtensor,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    // estimate of how big things are:
    data_t weak_field_var = emtensor.rho * emtensor.rho;
    FOR(i)
    {
        weak_field_var += emtensor.Si[i] * emtensor.Si[i];
        FOR(j) { weak_field_var += emtensor.Sij[i][j] * emtensor.Sij[i][j]; }
    }

    return sqrt(weak_field_var);
}

template <class System>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
data_t C2EFT<System>::weak_field_condition(
    const data_t &weak_field_var,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    const auto &vars = gq.get_vars();

    data_t condition_chi =
        sigmoid(vars.chi, m_params.chi_width, m_params.chi_threshold);

    data_t weak_field_condition =
        condition_chi * sigmoid(weak_field_var, -m_params.weak_field_width,
                                m_params.weak_field_threshold);

    return weak_field_condition;
}

template <class System>
template <class data_t, template <typename> class rhs_vars_t,
          template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t>
void C2EFT<System>::add_diffusion_terms(
    rhs_vars_t<data_t> &rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
    data_t diffCoeffSafe) const
{
    m_system.add_diffusion_terms(rhs, gq, diffCoeffSafe);
}

#endif /* C2EFT_IMPL_HPP_ */
