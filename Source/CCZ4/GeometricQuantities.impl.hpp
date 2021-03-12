/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(GEOMETRICQUANTITIES_HPP_)
#error "This file should only be included through GeometricQuantities.hpp"
#endif

#ifndef GEOMETRICQUANTITIES_IMPL_HPP_
#define GEOMETRICQUANTITIES_IMPL_HPP_

// helpers to avoid clutter in the code
#define template_GQ                                                            \
    template <class data_t, template <typename> class vars_t,                  \
              template <typename> class diff2_vars_t>
#define GeometricQuantities_t GeometricQuantities<data_t, vars_t, diff2_vars_t>

template_GQ GeometricQuantities_t::GeometricQuantities()
    : m_formulation(-1), m_ccz4_params(nullptr), m_16_pi_G_Newton(0.),
      m_cosmological_constant(0.)
{
    set_all_to_null();
}

template_GQ GeometricQuantities_t::GeometricQuantities(
    const Vars &a_vars, const Diff1Vars &a_d1_vars, const Diff2Vars &a_d2_vars)
    : GeometricQuantities()
{
    m_vars = &a_vars;
    m_d1_vars = &a_d1_vars;
    m_d2_vars = &a_d2_vars;
}

template_GQ void GeometricQuantities_t::set_all_to_null()
{
    m_vars = nullptr;
    m_d1_vars = nullptr;
    m_d2_vars = nullptr;
    m_advection = nullptr;
    m_em_tensor = nullptr;

    m_h_UU = nullptr;
    m_chris = nullptr;
    m_Z_U_conformal = nullptr;
    m_d1_chris_contracted = nullptr;
    m_covd_chi_conformal = nullptr;
    m_riemann_conformal_LLLL = nullptr;

    m_metric_spatial = nullptr;
    m_metric_UU_spatial = nullptr;
    m_shift_L = nullptr;
    m_extrinsic_curvature = nullptr;
    m_chris_spatial = nullptr;
    m_Z_U = nullptr;
    m_Z = nullptr;
    m_covd_Z = nullptr;
    m_d1_extrinsic_curvature = nullptr;
    m_covd_extrinsic_curvature = nullptr;
    m_levi_civita_spatial = nullptr;
    m_covd_lapse = nullptr;
    m_ricci = nullptr;
    m_ricci_1DZ = nullptr;
    m_ricci_2DZ = nullptr;
    m_weyl_magnetic_part = nullptr;
    m_riemann_spatial_LLLL = nullptr;
    m_gauss_codazzi = nullptr;
    m_codazzi_mainardi = nullptr;

    m_momentum_constraints = nullptr;
    m_weyl_electric_part = nullptr;
    m_lie_Z = nullptr;

    m_em_tensor_effective = nullptr;
    m_hamiltonian_constraint = nullptr;
    m_lie_derivatives = nullptr;
    m_lie_extrinsic_curvature = nullptr;
    m_eom_double_normal_projection = nullptr;

    m_rhs_equations = nullptr;
    // m_dt_chris_contracted = nullptr;
    // m_dt_chris_spatial_contracted = nullptr;

    m_metric_ST = nullptr;
    m_projector_LU_ST = nullptr;
    m_metric_UU_ST = nullptr;
    m_normal_U_ST = nullptr;
    m_normal_L_ST = nullptr;
    m_shift_ST = nullptr;
    m_levi_civita_spatial_ST = nullptr;
    m_levi_civita_ST = nullptr;
    m_Z_L_ST = nullptr;
    m_grad_normal_LL = nullptr;
    m_covd_Z_L_ST = nullptr;

    m_em_tensor_ST = nullptr;
    m_em_tensor_trace_ST = nullptr;
    m_weyl_tensor_LLLL = nullptr;
    m_weyl_squared = nullptr;

    m_em_tensor_effective_ST = nullptr;
    m_em_tensor_effective_trace_ST = nullptr;
    m_riemann_LLLL_ST = nullptr;
    m_riemann_LLLL_ST_v2 = nullptr;
    m_ricci_ST = nullptr;
    m_ricci_scalar_ST = nullptr;
    m_ricci_squared = nullptr;
    m_kretschmann = nullptr;
    m_riemann_squared = nullptr;

    // m_chris_ST = nullptr;
    // m_d1_Z_L_ST = nullptr;
}

template_GQ void GeometricQuantities_t::clean()
{
    if (m_h_UU != nullptr)
    {
        delete m_h_UU;
        m_h_UU = nullptr;
    }
    if (m_chris != nullptr)
    {
        delete m_chris;
        m_chris = nullptr;
    }
    if (m_Z_U_conformal != nullptr)
    {
        delete m_Z_U_conformal;
        m_Z_U_conformal = nullptr;
    }
    if (m_d1_chris_contracted != nullptr)
    {
        delete m_d1_chris_contracted;
        m_d1_chris_contracted = nullptr;
    }
    if (m_covd_chi_conformal != nullptr)
    {
        delete m_covd_chi_conformal;
        m_covd_chi_conformal = nullptr;
    }
    if (m_riemann_conformal_LLLL != nullptr)
    {
        delete m_riemann_conformal_LLLL;
        m_riemann_conformal_LLLL = nullptr;
    }

    if (m_metric_spatial != nullptr)
    {
        delete m_metric_spatial;
        m_metric_spatial = nullptr;
    }
    if (m_metric_UU_spatial != nullptr)
    {
        delete m_metric_UU_spatial;
        m_metric_UU_spatial = nullptr;
    }
    if (m_shift_L != nullptr)
    {
        delete m_shift_L;
        m_shift_L = nullptr;
    }
    if (m_extrinsic_curvature != nullptr)
    {
        delete m_extrinsic_curvature;
        m_extrinsic_curvature = nullptr;
    }
    if (m_chris_spatial != nullptr)
    {
        delete m_chris_spatial;
        m_chris_spatial = nullptr;
    }
    if (m_Z_U != nullptr)
    {
        delete m_Z_U;
        m_Z_U = nullptr;
    }
    if (m_Z != nullptr)
    {
        delete m_Z;
        m_Z = nullptr;
    }
    if (m_covd_Z != nullptr)
    {
        delete m_covd_Z;
        m_covd_Z = nullptr;
    }
    if (m_d1_extrinsic_curvature != nullptr)
    {
        delete m_d1_extrinsic_curvature;
        m_d1_extrinsic_curvature = nullptr;
    }
    if (m_covd_extrinsic_curvature != nullptr)
    {
        delete m_d1_extrinsic_curvature;
        m_d1_extrinsic_curvature = nullptr;
    }
    if (m_levi_civita_spatial != nullptr)
    {
        delete m_levi_civita_spatial;
        m_levi_civita_spatial = nullptr;
    }
    if (m_covd_lapse != nullptr)
    {
        delete m_covd_lapse;
        m_covd_lapse = nullptr;
    }
    if (m_ricci != nullptr)
    {
        delete m_ricci;
        m_ricci = nullptr;
    }
    if (m_ricci_1DZ != nullptr)
    {
        delete m_ricci_1DZ;
        m_ricci_1DZ = nullptr;
    }
    if (m_ricci_2DZ != nullptr)
    {
        delete m_ricci_2DZ;
        m_ricci_2DZ = nullptr;
    }
    if (m_weyl_magnetic_part != nullptr)
    {
        delete m_weyl_magnetic_part;
        m_weyl_magnetic_part = nullptr;
    }
    if (m_riemann_spatial_LLLL != nullptr)
    {
        delete m_riemann_spatial_LLLL;
        m_riemann_spatial_LLLL = nullptr;
    }
    if (m_gauss_codazzi != nullptr)
    {
        delete m_gauss_codazzi;
        m_gauss_codazzi = nullptr;
    }
    if (m_codazzi_mainardi != nullptr)
    {
        delete m_codazzi_mainardi;
        m_codazzi_mainardi = nullptr;
    }

    if (m_metric_ST != nullptr)
    {
        delete m_metric_ST;
        m_metric_ST = nullptr;
    }
    if (m_projector_LU_ST != nullptr)
    {
        delete m_projector_LU_ST;
        m_projector_LU_ST = nullptr;
    }
    if (m_metric_UU_ST != nullptr)
    {
        delete m_metric_UU_ST;
        m_metric_UU_ST = nullptr;
    }
    if (m_normal_U_ST != nullptr)
    {
        delete m_normal_U_ST;
        m_normal_U_ST = nullptr;
    }
    if (m_normal_L_ST != nullptr)
    {
        delete m_normal_L_ST;
        m_normal_L_ST = nullptr;
    }
    if (m_shift_ST != nullptr)
    {
        delete m_shift_ST;
        m_shift_ST = nullptr;
    }
    if (m_levi_civita_spatial_ST != nullptr)
    {
        delete m_levi_civita_spatial_ST;
        m_levi_civita_spatial_ST = nullptr;
    }
    if (m_levi_civita_ST != nullptr)
    {
        delete m_levi_civita_ST;
        m_levi_civita_ST = nullptr;
    }
    if (m_Z_L_ST != nullptr)
    {
        delete m_Z_L_ST;
        m_Z_L_ST = nullptr;
    }
    if (m_grad_normal_LL != nullptr)
    {
        delete m_grad_normal_LL;
        m_grad_normal_LL = nullptr;
    }
    if (m_covd_Z_L_ST != nullptr)
    {
        delete m_covd_Z_L_ST;
        m_covd_Z_L_ST = nullptr;
    }

    clean_em_tensor_dependent();
}

template_GQ void GeometricQuantities_t::clean_em_tensor_dependent()
{
    if (m_momentum_constraints != nullptr)
    {
        delete m_momentum_constraints;
        m_momentum_constraints = nullptr;
    }
    if (m_weyl_electric_part != nullptr)
    {
        delete m_weyl_electric_part;
        m_weyl_electric_part = nullptr;
    }
    if (m_lie_Z != nullptr)
    {
        delete m_lie_Z;
        m_lie_Z = nullptr;
    }

    if (m_em_tensor_ST != nullptr)
    {
        delete m_em_tensor_ST;
        m_em_tensor_ST = nullptr;
    }
    if (m_em_tensor_trace_ST != nullptr)
    {
        delete m_em_tensor_trace_ST;
        m_em_tensor_trace_ST = nullptr;
    }
    if (m_weyl_tensor_LLLL != nullptr)
    {
        delete m_weyl_tensor_LLLL;
        m_weyl_tensor_LLLL = nullptr;
    }
    if (m_weyl_squared != nullptr)
    {
        delete m_weyl_squared;
        m_weyl_squared = nullptr;
    }

    clean_eom_dependent();
}

template_GQ void GeometricQuantities_t::clean_eom_dependent()
{
    if (m_em_tensor_effective != nullptr)
    {
        delete m_em_tensor_effective;
        m_em_tensor_effective = nullptr;
    }
    if (m_hamiltonian_constraint != nullptr)
    {
        delete m_hamiltonian_constraint;
        m_hamiltonian_constraint = nullptr;
    }
    if (m_lie_derivatives != nullptr)
    {
        delete m_lie_derivatives;
        m_lie_derivatives = nullptr;
    }
    if (m_lie_extrinsic_curvature != nullptr)
    {
        delete m_lie_extrinsic_curvature;
        m_lie_extrinsic_curvature = nullptr;
    }
    if (m_eom_double_normal_projection != nullptr)
    {
        delete m_eom_double_normal_projection;
        m_eom_double_normal_projection = nullptr;
    }

    if (m_em_tensor_effective_ST != nullptr)
    {
        delete m_em_tensor_effective_ST;
        m_em_tensor_effective_ST = nullptr;
    }
    if (m_em_tensor_effective_trace_ST != nullptr)
    {
        delete m_em_tensor_effective_trace_ST;
        m_em_tensor_effective_trace_ST = nullptr;
    }
    if (m_riemann_LLLL_ST != nullptr)
    {
        delete m_riemann_LLLL_ST;
        m_riemann_LLLL_ST = nullptr;
    }
    if (m_riemann_LLLL_ST_v2 != nullptr)
    {
        delete m_riemann_LLLL_ST_v2;
        m_riemann_LLLL_ST_v2 = nullptr;
    }
    if (m_ricci_ST != nullptr)
    {
        delete m_ricci_ST;
        m_ricci_ST = nullptr;
    }
    if (m_ricci_scalar_ST != nullptr)
    {
        delete m_ricci_scalar_ST;
        m_ricci_scalar_ST = nullptr;
    }
    if (m_ricci_squared != nullptr)
    {
        delete m_ricci_squared;
        m_ricci_squared = nullptr;
    }
    if (m_kretschmann != nullptr)
    {
        delete m_kretschmann;
        m_kretschmann = nullptr;
    }
    if (m_riemann_squared != nullptr)
    {
        delete m_riemann_squared;
        m_riemann_squared = nullptr;
    }

    clean_advection_dependent();
}

template_GQ void GeometricQuantities_t::clean_advection_dependent()
{
    if (m_rhs_equations != nullptr)
    {
        delete m_rhs_equations;
        m_rhs_equations = nullptr;
    }
    /*
    if (m_dt_chris_contracted != nullptr)
    {
        delete m_dt_chris_contracted;
        m_dt_chris_contracted = nullptr;
    }
    if (m_dt_chris_spatial_contracted != nullptr)
    {
        delete m_dt_chris_spatial_contracted;
        m_dt_chris_spatial_contracted = nullptr;
    }

    if (m_chris_ST != nullptr)
    {
        delete m_chris_ST;
        m_chris_ST = nullptr;
    }
    if (m_d1_Z_L_ST != nullptr)
    {
        delete m_d1_Z_L_ST;
        m_d1_Z_L_ST = nullptr;
    }
    */
}

template_GQ GeometricQuantities_t::~GeometricQuantities() { clean(); }

//////////////////////////////////////////////////////////////////////////
///////////////////////////////    SETS    ///////////////////////////////
//////////////////////////////////////////////////////////////////////////

template_GQ void GeometricQuantities_t::set_vars(const Vars &a_vars)
{
    m_vars = &a_vars;
    clean();
}
template_GQ void GeometricQuantities_t::set_d1_vars(const Diff1Vars &a_d1_vars)
{
    m_d1_vars = &a_d1_vars;
    clean();
}
template_GQ void GeometricQuantities_t::set_d2_vars(const Diff2Vars &a_d2_vars)
{
    m_d2_vars = &a_d2_vars;
    clean();
}
template_GQ void GeometricQuantities_t::set_advection(const Vars &a_advection)
{
    m_advection = &a_advection;
    clean_advection_dependent();
}
template_GQ void
GeometricQuantities_t::set_formulation(int formulation,
                                       const CCZ4::params_t &a_ccz4_params)
{
    m_formulation = formulation;
    m_ccz4_params = &a_ccz4_params;
    clean_eom_dependent();
}
template_GQ void
GeometricQuantities_t::set_em_tensor(const emtensor_t<data_t> &a_em_tensor,
                                     double G_Newton)
{
    m_em_tensor = &a_em_tensor;
    CH_assert(G_Newton != 0.);
    m_16_pi_G_Newton = 16. * M_PI * G_Newton;
    clean_em_tensor_dependent();
}
template_GQ void
GeometricQuantities_t::set_cosmological_constant(double cosmological_constant)
{
    m_cosmological_constant = cosmological_constant;
    clean_eom_dependent();
}

template_GQ void GeometricQuantities_t::set_all_vars(const Vars &a_vars,
                                                     const Diff1Vars &a_d1_vars,
                                                     const Diff2Vars &a_d2_vars)
{
    m_vars = &a_vars;
    m_d1_vars = &a_d1_vars;
    m_d2_vars = &a_d2_vars;
    clean();
}

//////////////////////////////////////////////////////////////////////////
///////////////////////////////    GETS    ///////////////////////////////
//////////////////////////////////////////////////////////////////////////

template_GQ const double GeometricQuantities_t::get_formulation() const
{
    return m_formulation;
}
template_GQ const double
GeometricQuantities_t::get_cosmological_constant() const
{
    return m_cosmological_constant;
}
template_GQ const typename GeometricQuantities_t::Vars &
GeometricQuantities_t::get_vars() const
{
    CH_assert(m_vars != nullptr); // CH_assert of MayDay::Error ?
    return *m_vars;
}
template_GQ const typename GeometricQuantities_t::Diff1Vars &
GeometricQuantities_t::get_d1_vars() const
{
    CH_assert(m_d1_vars != nullptr);
    return *m_d1_vars;
}
template_GQ const typename GeometricQuantities_t::Diff2Vars &
GeometricQuantities_t::get_d2_vars() const
{
    CH_assert(m_d2_vars != nullptr);
    return *m_d2_vars;
}
template_GQ const typename GeometricQuantities_t::Vars &
GeometricQuantities_t::get_advection() const
{
    CH_assert(m_advection != nullptr);
    return *m_advection;
}
template_GQ const emtensor_t<data_t> &
GeometricQuantities_t::get_em_tensor() const
{
    CH_assert(m_em_tensor != nullptr);
    return *m_em_tensor;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_h_UU()
{
    if (m_h_UU == nullptr)
        compute_h_UU();
    return *m_h_UU;
}
template_GQ const chris_t<data_t> &GeometricQuantities_t::get_chris()
{
    if (m_chris == nullptr)
        compute_chris();
    return *m_chris;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Z_U_conformal()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_Z_U_conformal == nullptr)
        compute_Z_U_conformal();
    return *m_Z_U_conformal;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_d1_chris_contracted()
{
    if (m_d1_chris_contracted == nullptr)
        compute_d1_chris_contracted();
    return *m_d1_chris_contracted;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_covd_chi_conformal()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_covd_chi_conformal == nullptr)
        compute_covd_chi_conformal();
    return *m_covd_chi_conformal;
}
template_GQ const Tensor<4, data_t> &
GeometricQuantities_t::get_riemann_conformal_LLLL()
{
    if (m_riemann_conformal_LLLL == nullptr)
        compute_riemann_conformal_LLLL();
    return *m_riemann_conformal_LLLL;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_metric_spatial()
{
    if (m_metric_spatial == nullptr)
        compute_metric_spatial();
    return *m_metric_spatial;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_metric_UU_spatial()
{
    if (m_metric_UU_spatial == nullptr)
        compute_metric_UU_spatial();
    return *m_metric_UU_spatial;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_shift_L()
{
    if (m_shift_L == nullptr)
        compute_shift_L();
    return *m_shift_L;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_extrinsic_curvature()
{
    if (m_extrinsic_curvature == nullptr)
        compute_extrinsic_curvature();
    return *m_extrinsic_curvature;
}
template_GQ const Tensor<3, data_t> &GeometricQuantities_t::get_chris_spatial()
{
    if (m_chris_spatial == nullptr)
        compute_chris_spatial();
    return *m_chris_spatial;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Z_U()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_Z_U == nullptr)
        compute_Z_U();
    return *m_Z_U;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Z()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_Z == nullptr)
        compute_Z();
    return *m_Z;
}
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_covd_Z()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_covd_Z == nullptr)
        compute_covd_Z();
    return *m_covd_Z;
}
template_GQ const Tensor<3, data_t> &
GeometricQuantities_t::get_d1_extrinsic_curvature()
{
    if (m_d1_extrinsic_curvature == nullptr)
        compute_d1_extrinsic_curvature();
    return *m_d1_extrinsic_curvature;
}
template_GQ const Tensor<3, data_t> &
GeometricQuantities_t::get_covd_extrinsic_curvature()
{
    if (m_covd_extrinsic_curvature == nullptr)
        compute_covd_extrinsic_curvature();
    return *m_covd_extrinsic_curvature;
}
template_GQ const Tensor<3, data_t> &
GeometricQuantities_t::get_levi_civita_spatial()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_levi_civita_spatial == nullptr)
        compute_levi_civita_spatial();
    return *m_levi_civita_spatial;
}
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_covd_lapse()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_covd_lapse == nullptr)
        compute_covd_lapse();
    return *m_covd_lapse;
}
template_GQ const ricci_t<data_t> &GeometricQuantities_t::get_ricci_qDZ(int q)
{
    // either q==0 to get the pure Ricci, set BSSN or CCZ4 (for
    // example for q=2 one will get the Ricci with calculated Gammas
    // replaced by evolved Gammas, and extra Z terms for CCZ4)
    CH_assert(q == 0 || m_formulation >= 0);
    switch (q)
    {
    case 0:
        return get_ricci();
    case 1:
        return get_ricci_1DZ();
    case 2:
        return get_ricci_2DZ();
    default:
        MayDay::Error("Mmmm. Use 'compute_ricci_qDZ' for q != 0,1,2 .");
    }
}
template_GQ const ricci_t<data_t> &GeometricQuantities_t::get_ricci()
{
    if (m_ricci == nullptr)
        compute_ricci();
    return *m_ricci;
}
template_GQ const ricci_t<data_t> &GeometricQuantities_t::get_ricci_1DZ()
{
    CH_assert(m_formulation >= 0);
    if (m_ricci_1DZ == nullptr)
        compute_ricci_1DZ();
    return *m_ricci_1DZ;
}
template_GQ const ricci_t<data_t> &GeometricQuantities_t::get_ricci_2DZ()
{
    // will give the Ricci with calculated Gammas
    // replaced by evolved Gammas, and extra Z terms for CCZ4
    CH_assert(m_formulation >= 0);
    if (m_ricci_2DZ == nullptr)
        compute_ricci_2DZ();
    return *m_ricci_2DZ;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_weyl_magnetic_part()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_weyl_magnetic_part == nullptr)
        compute_weyl_magnetic_part();
    return *m_weyl_magnetic_part;
}
template_GQ const Tensor<4, data_t> &
GeometricQuantities_t::get_riemann_spatial_LLLL()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_riemann_spatial_LLLL == nullptr)
        compute_riemann_spatial_LLLL();
    return *m_riemann_spatial_LLLL;
}
template_GQ const Tensor<4, data_t> &GeometricQuantities_t::get_gauss_codazzi()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_gauss_codazzi == nullptr)
        compute_gauss_codazzi();
    return *m_gauss_codazzi;
}
template_GQ const Tensor<3, data_t> &
GeometricQuantities_t::get_codazzi_mainardi()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_codazzi_mainardi == nullptr)
        compute_codazzi_mainardi();
    return *m_codazzi_mainardi;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_momentum_constraints()
{
    if (m_momentum_constraints == nullptr)
        compute_momentum_constraints();
    return *m_momentum_constraints;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_weyl_electric_part()
{
    CH_assert(GR_SPACEDIM == 3);
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_weyl_electric_part == nullptr)
        compute_weyl_electric_part();
    return *m_weyl_electric_part;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_lie_Z()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_lie_Z == nullptr)
        compute_lie_Z();
    return *m_lie_Z;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const emtensor_t<data_t> &
GeometricQuantities_t::get_em_tensor_effective()
{
    // uses cosmological constant
    if (m_em_tensor_effective == nullptr)
        compute_em_tensor_effective();
    return *m_em_tensor_effective;
}
template_GQ const data_t &GeometricQuantities_t::get_hamiltonian_constraint()
{
    // uses cosmological constant
    if (m_hamiltonian_constraint == nullptr)
        compute_hamiltonian_constraint();
    return *m_hamiltonian_constraint;
}
template_GQ const vars_t<data_t> &GeometricQuantities_t::get_lie_derivatives()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_lie_derivatives == nullptr)
        compute_lie_derivatives();
    return *m_lie_derivatives;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_lie_extrinsic_curvature()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_lie_extrinsic_curvature == nullptr)
        compute_lie_extrinsic_curvature();
    return *m_lie_extrinsic_curvature;
}
template_GQ const Tensor<2, data_t> &
GeometricQuantities_t::get_eom_double_normal_projection()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_eom_double_normal_projection == nullptr)
        compute_eom_double_normal_projection();
    return *m_eom_double_normal_projection;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const vars_t<data_t> &GeometricQuantities_t::get_rhs_equations()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_rhs_equations == nullptr)
        compute_rhs_equations();
    return *m_rhs_equations;
}
/*
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_dt_chris_contracted()
{
    if (m_dt_chris_contracted == nullptr)
        compute_dt_chris_contracted();
    return *m_dt_chris_contracted;
}
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_dt_chris_spatial_contracted()
{
    if (m_dt_chris_spatial_contracted == nullptr)
        compute_dt_chris_spatial_contracted();
    return *m_dt_chris_spatial_contracted;
}
*/
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_metric_ST()
{
    if (m_metric_ST == nullptr)
        compute_metric_ST();
    return *m_metric_ST;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_projector_LU_ST()
{
    if (m_projector_LU_ST == nullptr)
        compute_projector_LU_ST();
    return *m_projector_LU_ST;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_metric_UU_ST()
{
    if (m_metric_UU_ST == nullptr)
        compute_metric_UU_ST();
    return *m_metric_UU_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_normal_U_ST()
{
    if (m_normal_U_ST == nullptr)
        compute_normal_U_ST();
    return *m_normal_U_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_normal_L_ST()
{
    if (m_normal_L_ST == nullptr)
        compute_normal_L_ST();
    return *m_normal_L_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_shift_ST()
{
    if (m_shift_ST == nullptr)
        compute_shift_ST();
    return *m_shift_ST;
}
template_GQ const Tensor<3, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_levi_civita_spatial_ST()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_levi_civita_spatial_ST == nullptr)
        compute_levi_civita_spatial_ST();
    return *m_levi_civita_spatial_ST;
}
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_levi_civita_ST()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_levi_civita_ST == nullptr)
        compute_levi_civita_ST();
    return *m_levi_civita_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_Z_L_ST()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_Z_L_ST == nullptr)
        compute_Z_L_ST();
    return *m_Z_L_ST;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_grad_normal_LL()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_grad_normal_LL == nullptr)
        compute_grad_normal_LL();
    return *m_grad_normal_LL;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_covd_Z_L_ST()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_covd_Z_L_ST == nullptr)
        compute_covd_Z_L_ST();
    return *m_covd_Z_L_ST;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_em_tensor_ST()
{
    if (m_em_tensor_ST == nullptr)
        compute_em_tensor_ST();
    return *m_em_tensor_ST;
}
template_GQ const data_t &GeometricQuantities_t::get_em_tensor_trace_ST()
{
    if (m_em_tensor_trace_ST == nullptr)
        compute_em_tensor_trace_ST();
    return *m_em_tensor_trace_ST;
}
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_weyl_tensor_LLLL()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_weyl_tensor_LLLL == nullptr)
        compute_weyl_tensor_LLLL();
    return *m_weyl_tensor_LLLL;
}
template_GQ const data_t &GeometricQuantities_t::get_weyl_squared()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_weyl_squared == nullptr)
        compute_weyl_squared();
    return *m_weyl_squared;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_em_tensor_effective_ST()
{
    if (m_em_tensor_effective_ST == nullptr)
        compute_em_tensor_effective_ST();
    return *m_em_tensor_effective_ST;
}
template_GQ const data_t &
GeometricQuantities_t::get_em_tensor_effective_trace_ST()
{
    if (m_em_tensor_effective_trace_ST == nullptr)
        compute_em_tensor_effective_trace_ST();
    return *m_em_tensor_effective_trace_ST;
}
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_riemann_LLLL_ST()
{
    CH_assert(m_formulation >= 0); // formulation is set
    CH_assert(GR_SPACEDIM == 3);
    if (m_riemann_LLLL_ST == nullptr)
        compute_riemann_LLLL_ST();
    return *m_riemann_LLLL_ST;
}
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_riemann_LLLL_ST_v2()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_riemann_LLLL_ST_v2 == nullptr)
        compute_riemann_LLLL_ST_v2();
    return *m_riemann_LLLL_ST_v2;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_ricci_ST()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_ricci_ST == nullptr)
        compute_ricci_ST();
    return *m_ricci_ST;
}
template_GQ const data_t &GeometricQuantities_t::get_ricci_scalar_ST()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_ricci_scalar_ST == nullptr)
        compute_ricci_scalar_ST();
    return *m_ricci_scalar_ST;
}
template_GQ const data_t &GeometricQuantities_t::get_ricci_squared()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_ricci_squared == nullptr)
        compute_ricci_squared();
    return *m_ricci_squared;
}
template_GQ const data_t &GeometricQuantities_t::get_kretschmann()
{
    CH_assert(m_formulation >= 0); // formulation is set
    CH_assert(GR_SPACEDIM == 3);
    if (m_kretschmann == nullptr)
        compute_kretschmann();
    return *m_kretschmann;
}
template_GQ const data_t &GeometricQuantities_t::get_riemann_squared()
{
    CH_assert(m_formulation >= 0); // formulation is set
    CH_assert(GR_SPACEDIM == 3);
    if (m_riemann_squared == nullptr)
        compute_riemann_squared();
    return *m_riemann_squared;
}
//////////////////////////////////////////////////////////////////////////
/*
template_GQ const Tensor<3, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_chris_ST()
{
    if (m_chris_ST == nullptr)
        compute_chris_ST();
    return *m_chris_ST;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_d1_Z_L_ST()
{
    CH_assert(m_formulation == CCZ4::USE_CCZ4);
    if (m_d1_Z_L_ST == nullptr)
        compute_d1_Z_L_ST();
    return *m_d1_Z_L_ST;
}
*/
//////////////////////////////////////////////////////////////////////////
/////////////////////////////    COMPUTES    /////////////////////////////
//////////////////////////////////////////////////////////////////////////

template_GQ void GeometricQuantities_t::compute_h_UU()
{
    if (m_h_UU != nullptr)
        delete m_h_UU;
    m_h_UU = new const Tensor<2, data_t>(
        TensorAlgebra::compute_inverse_sym(get_vars().h));
}
template_GQ void GeometricQuantities_t::compute_chris()
{
    if (m_chris != nullptr)
        delete m_chris;
    m_chris = new const chris_t<data_t>(
        TensorAlgebra::compute_christoffel(get_d1_vars().h, get_h_UU()));
}
template_GQ void GeometricQuantities_t::compute_Z_U_conformal()
{
    if (m_Z_U_conformal != nullptr)
        delete m_Z_U_conformal;

    const auto &vars = get_vars();
    const auto &chris = get_chris();

    Tensor<1, data_t> Z_U_conformal;
    FOR(i) { Z_U_conformal[i] = 0.5 * (vars.Gamma[i] - chris.contracted[i]); }

    m_Z_U_conformal = new const Tensor<1, data_t>(Z_U_conformal);
}
template_GQ void GeometricQuantities_t::compute_d1_chris_contracted()
{
    if (m_d1_chris_contracted != nullptr)
        delete m_d1_chris_contracted;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();

    Tensor<2, data_t> d1_chris_contracted = 0.0;
    FOR(i, j)
    {
        FOR(m, n, p)
        {
            data_t d1_terms = 0.0;
            FOR(q, r)
            {
                d1_terms += -h_UU[q][r] * (d1.h[n][q][j] * d1.h[m][p][r] +
                                           d1.h[m][n][j] * d1.h[p][q][r]);
            }
            d1_chris_contracted[i][j] +=
                h_UU[i][m] * h_UU[n][p] * (d2.h[m][n][j][p] + d1_terms);
        }
    }

    m_d1_chris_contracted = new const Tensor<2, data_t>(d1_chris_contracted);
}
template_GQ void GeometricQuantities_t::compute_covd_chi_conformal()
{
    if (m_covd_chi_conformal != nullptr)
        delete m_covd_chi_conformal;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &chris = get_chris();

    m_covd_chi_conformal = new const Tensor<2, data_t>(
        TensorAlgebra::covariant_derivative(d2.chi, d1.chi, chris.ULL));
}
template_GQ void GeometricQuantities_t::compute_riemann_conformal_LLLL()
{
    if (m_riemann_conformal_LLLL != nullptr)
        delete m_riemann_conformal_LLLL;

    Tensor<4, data_t> riemann_conformal_LLLL;

    const auto &chris = get_chris();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();

    FOR(i, j, k, l)
    {
        riemann_conformal_LLLL[i][j][k][l] =
            0.5 * (d2.h[i][l][j][k] + d2.h[j][i][l][k] - d2.h[j][l][i][k]) -
            0.5 * (d2.h[i][k][j][l] + d2.h[j][i][k][l] - d2.h[j][k][i][l]);
        FOR(m)
        {
            riemann_conformal_LLLL[i][j][k][l] +=
                chris.ULL[m][j][k] * d1.h[i][m][l] -
                chris.ULL[m][j][l] * d1.h[i][m][k] +
                chris.LLL[i][m][k] * chris.ULL[m][l][j] -
                chris.LLL[i][m][l] * chris.ULL[m][k][j];
        }
    }

    m_riemann_conformal_LLLL =
        new const Tensor<4, data_t>(riemann_conformal_LLLL);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_metric_spatial()
{
    if (m_metric_spatial != nullptr)
        delete m_metric_spatial;

    Tensor<2, data_t> metric_spatial;

    const auto &vars = get_vars();

    FOR(i, j) { metric_spatial[i][j] = vars.h[i][j] / vars.chi; }
    m_metric_spatial = new Tensor<2, data_t>(metric_spatial);
}
template_GQ void GeometricQuantities_t::compute_metric_UU_spatial()
{
    if (m_metric_UU_spatial != nullptr)
        delete m_metric_UU_spatial;

    Tensor<2, data_t> metric_UU_spatial;

    const auto &vars = get_vars();
    const auto &h_UU = get_h_UU();

    FOR(i, j) { metric_UU_spatial[i][j] = h_UU[i][j] * vars.chi; }
    m_metric_UU_spatial = new Tensor<2, data_t>(metric_UU_spatial);
}
template_GQ void GeometricQuantities_t::compute_shift_L()
{
    if (m_shift_L != nullptr)
        delete m_shift_L;

    Tensor<1, data_t> shift_L;

    const auto &vars = get_vars();
    const auto &metric_spatial = get_metric_spatial();
    m_shift_L = new const Tensor<1, data_t>(
        TensorAlgebra::lower_all(vars.shift, metric_spatial));
}
template_GQ void GeometricQuantities_t::compute_extrinsic_curvature()
{
    if (m_extrinsic_curvature != nullptr)
        delete m_extrinsic_curvature;

    Tensor<2, data_t> extrinsic_curvature;

    const auto &vars = get_vars();

    FOR(i, j)
    {
        extrinsic_curvature[i][j] =
            (vars.A[i][j] + vars.K * vars.h[i][j] / GR_SPACEDIM) / vars.chi;
    }

    m_extrinsic_curvature = new Tensor<2, data_t>(extrinsic_curvature);
}
template_GQ void GeometricQuantities_t::compute_chris_spatial()
{
    if (m_chris_spatial != nullptr)
        delete m_chris_spatial;

    const auto &chris = get_chris();
    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &h_UU = get_h_UU();

    m_chris_spatial =
        new const Tensor<3, data_t>(TensorAlgebra::compute_phys_chris(
            d1.chi, vars.chi, vars.h, h_UU, chris.ULL));
}
template_GQ void GeometricQuantities_t::compute_Z_U()
{
    if (m_Z_U != nullptr)
        delete m_Z_U;

    const auto &vars = get_vars();

    Tensor<1, data_t> Z_U = get_Z_U_conformal();
    FOR(i) { Z_U[i] *= vars.chi; }

    m_Z_U = new const Tensor<1, data_t>(Z_U);
}
template_GQ void GeometricQuantities_t::compute_Z()
{
    if (m_Z != nullptr)
        delete m_Z;

    const auto &vars = get_vars();
    const auto &Z_U_conformal = get_Z_U_conformal();

    m_Z = new const Tensor<1, data_t>(
        TensorAlgebra::compute_dot_product(Z_U_conformal, vars.h));
}
template_GQ void GeometricQuantities_t::compute_covd_Z()
{
    if (m_covd_Z != nullptr)
        delete m_covd_Z;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &Z = get_Z();
    const auto &Z_U_conformal = get_Z_U_conformal();
    const auto &d1_chris_contracted = get_d1_chris_contracted();
    const auto &chris = get_chris();

    Tensor<2, data_t> d_Z_U_conformal; // \partial_j Z^i
    FOR(i, j)
    {
        d_Z_U_conformal[i][j] =
            0.5 * (d1.Gamma[i][j] - d1_chris_contracted[i][j]);
    }

    const Tensor<2, data_t> covd_Z_U_conformal =
        TensorAlgebra::covariant_derivative(d_Z_U_conformal, Z_U_conformal,
                                            chris.ULL, {true});

    const Tensor<2, data_t> covd_Z_U_conformal_dot_h =
        TensorAlgebra::compute_dot_product(vars.h, covd_Z_U_conformal, 0, 0);

    const data_t Z_U_dot_chi =
        TensorAlgebra::compute_dot_product(Z_U_conformal, d1.chi);

    Tensor<2, data_t> covd_Z;
    FOR(i, j)
    {
        covd_Z[i][j] =
            covd_Z_U_conformal_dot_h[i][j] +
            (Z[i] * d1.chi[j] + Z[j] * d1.chi[i] - vars.h[i][j] * Z_U_dot_chi) /
                (2. * vars.chi);
    }

    m_covd_Z = new const Tensor<2, data_t>(covd_Z);
}
template_GQ void GeometricQuantities_t::compute_d1_extrinsic_curvature()
{
    if (m_d1_extrinsic_curvature != nullptr)
        delete m_d1_extrinsic_curvature;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &K_tensor = get_extrinsic_curvature();

    Tensor<3, data_t> d1_K_tensor;
    FOR(i, j, k)
    {
        d1_K_tensor[i][j][k] = (d1.A[i][j][k] - d1.chi[k] * K_tensor[i][j] +
                                d1.h[i][j][k] * vars.K / GR_SPACEDIM +
                                vars.h[i][j] * d1.K[k] / GR_SPACEDIM) /
                               vars.chi;
    }

    m_d1_extrinsic_curvature = new const Tensor<3, data_t>(d1_K_tensor);
}
template_GQ void GeometricQuantities_t::compute_covd_extrinsic_curvature()
{
    if (m_covd_extrinsic_curvature != nullptr)
        delete m_covd_extrinsic_curvature;

    const auto &chris_spatial = get_chris_spatial();
    const auto &K_tensor = get_extrinsic_curvature();
    const auto &d1_K_tensor = get_d1_extrinsic_curvature();

    m_covd_extrinsic_curvature =
        new const Tensor<3, data_t>(TensorAlgebra::covariant_derivative(
            d1_K_tensor, K_tensor, chris_spatial));
}
template_GQ void GeometricQuantities_t::compute_levi_civita_spatial()
{
    if (m_levi_civita_spatial != nullptr)
        delete m_levi_civita_spatial;

    Tensor<3, data_t> levi_civita_spatial = {0.};

    const auto &levi_civita_ST = get_levi_civita_ST();
    const auto &n_U = get_normal_U_ST();

    FOR(i, j, k)
    {
        FOR_ST(l)
        {
            levi_civita_spatial[i][j][k] +=
                n_U[l] * levi_civita_ST[l][i + 1][j + 1][k + 1];
        }
    }

    m_levi_civita_spatial = new Tensor<3, data_t>(levi_civita_spatial);
}
template_GQ void GeometricQuantities_t::compute_covd_lapse()
{
    if (m_covd_lapse != nullptr)
        delete m_covd_lapse;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &chris_spatial = get_chris_spatial();

    m_covd_lapse = new const Tensor<2, data_t>(
        TensorAlgebra::covariant_derivative(d2.lapse, d1.lapse, chris_spatial));
}
template_GQ ricci_t<data_t> GeometricQuantities_t::compute_ricci_qDZ(int q)
{
    CH_assert(q == 0 || m_formulation >= 0);
    // either q==0 to get the pure Ricci, set BSSN or CCZ4 (for
    // example for q=2 one will get the Ricci with calculated Gammas
    // replaced by evolved Gammas, and extra Z terms for CCZ4)

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();
    const auto &chris = get_chris();

    ricci_t<data_t> ricci_qDZ;

    const Tensor<2, data_t> covDtilde2chi = get_covd_chi_conformal();
    const data_t boxtildechi =
        TensorAlgebra::compute_trace(covDtilde2chi, h_UU);
    const data_t dchi_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.chi, d1.chi, h_UU);

    Tensor<1, data_t> chris_q;
    Tensor<1, Tensor<1, data_t>> d1_chris_q;
    if (q != 2)
    {
        const auto &d1_chris_contracted = get_d1_chris_contracted();
        double q_over_2 = q / 2.;
        double one_m_q_over_2 = 1. - q_over_2;
        FOR(i)
        {
            chris_q[i] =
                q_over_2 * vars.Gamma[i] + one_m_q_over_2 * chris.contracted[i];
            FOR(j)
            {
                d1_chris_q[i][j] = q_over_2 * d1.Gamma[i][j] +
                                   one_m_q_over_2 * d1_chris_contracted[i][j];
            }
        }
    }
    else
    {
        chris_q = vars.Gamma;
        d1_chris_q = d1.Gamma;
    }

    Tensor<1, data_t> Z_U_q = {0.};
    if (q != 0 && m_formulation == CCZ4::USE_CCZ4) // only if in CCZ4
        TensorAlgebra::hard_copy(Z_U_q, get_Z_U_conformal());

    Tensor<3, data_t> chris_LLU = {0.};
    FOR4(i, j, k, l) { chris_LLU[i][j][k] += h_UU[k][l] * chris.LLL[i][j][l]; }

    FOR(i, j)
    {
        data_t ricci_tilde = 0;
        FOR(k)
        {
            ricci_tilde += 0.5 * (vars.h[k][i] * d1_chris_q[k][j] +
                                  vars.h[k][j] * d1_chris_q[k][i]);
            ricci_tilde +=
                0.5 * chris_q[k] * (chris.LLL[i][j][k] + chris.LLL[j][i][k]);
            // 0.5 * chris_q[k] * d1.h[i][j][k];
            FOR(l)
            {
                ricci_tilde += -0.5 * h_UU[k][l] * d2.h[i][j][k][l] +
                               (chris.ULL[k][l][i] * chris_LLU[j][k][l] +
                                chris.ULL[k][l][j] * chris_LLU[i][k][l] +
                                chris.ULL[k][i][l] * chris_LLU[k][j][l]);
                // ricci_tilde -= 0.5 * h_UU[k][l] * d2.h[i][j][k][l];
                // FOR(m)
                // {
                //     ricci_tilde +=
                //         h_UU[l][m] * (chris.ULL[k][l][i] * chris.LLL[j][k][m]
                //         +
                //                       chris.ULL[k][l][j] * chris.LLL[i][k][m]
                //                       + chris.ULL[k][i][m] *
                //                       chris.LLL[k][l][j]);
                // }
            }
        }

        const data_t ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covDtilde2chi[i][j] +
                   vars.h[i][j] * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1.chi[i] * d1.chi[j] +
                    GR_SPACEDIM * vars.h[i][j] * dchi_dot_dchi) /
                       (2 * vars.chi));

        data_t z_terms = 0.;
        if (q != 0 && m_formulation == CCZ4::USE_CCZ4)
        {
            FOR(k)
            {
                z_terms += Z_U_q[k] * (vars.h[i][k] * d1.chi[j] +
                                       vars.h[j][k] * d1.chi[i] -
                                       vars.h[i][j] * d1.chi[k]);
            }
        }

        ricci_qDZ.LL[i][j] =
            ricci_tilde + (ricci_chi + q / 2. * z_terms) / vars.chi;
    }

    ricci_qDZ.scalar =
        vars.chi * TensorAlgebra::compute_trace(ricci_qDZ.LL, h_UU);

    return ricci_qDZ;
}
template_GQ void GeometricQuantities_t::compute_ricci()
{
    if (m_ricci != nullptr)
        delete m_ricci;

    m_ricci = new const ricci_t<data_t>(compute_ricci_qDZ(0));
}
template_GQ void GeometricQuantities_t::compute_ricci_1DZ()
{
    if (m_ricci_1DZ != nullptr)
        delete m_ricci_1DZ;

    m_ricci_1DZ = new const ricci_t<data_t>(compute_ricci_qDZ(1));
}
template_GQ void GeometricQuantities_t::compute_ricci_2DZ()
{
    if (m_ricci_2DZ != nullptr)
        delete m_ricci_2DZ;

    m_ricci_2DZ = new const ricci_t<data_t>(compute_ricci_qDZ(2));
}
template_GQ void GeometricQuantities_t::compute_weyl_magnetic_part()
{
    if (m_weyl_magnetic_part != nullptr)
        delete m_weyl_magnetic_part;

    Tensor<3, data_t> levi_civita_spatial_LUU = {0.};
    const auto &levi_civita_spatial = get_levi_civita_spatial();
    const auto &h_UU_spatial = get_metric_UU_spatial();
    FOR(i, j, k)
    {
        FOR(m, n)
        {
            levi_civita_spatial_LUU[i][j][k] += levi_civita_spatial[i][m][n] *
                                                h_UU_spatial[m][j] *
                                                h_UU_spatial[n][k];
        }
    }

    const auto &covD_K_tensor = get_covd_extrinsic_curvature();
    Tensor<2, data_t> weyl_magnetic_part = {0.};
    FOR4(i, j, k, l)
    {
        weyl_magnetic_part[i][j] +=
            levi_civita_spatial_LUU[i][k][l] * covD_K_tensor[j][l][k];
    }

    TensorAlgebra::make_symmetric(weyl_magnetic_part);

    m_weyl_magnetic_part = new const Tensor<2, data_t>(weyl_magnetic_part);
}
template_GQ void GeometricQuantities_t::compute_riemann_spatial_LLLL()
{
    if (m_riemann_spatial_LLLL != nullptr)
        delete m_riemann_spatial_LLLL;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &h_UU = get_h_UU();
    const auto &riemann_conformal_LLLL = get_riemann_conformal_LLLL();
    const Tensor<2, data_t> covd_chi = get_covd_chi_conformal();

    const data_t &chi = vars.chi;
    const data_t chi3 = chi * chi * chi;
    const Tensor<2, data_t> &h = vars.h;
    const data_t dchi_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.chi, d1.chi, h_UU);

    Tensor<4, data_t> riemann_spatial_LLLL;

    FOR(i, j, k, l)
    {
        riemann_spatial_LLLL[i][j][k][l] =
            riemann_conformal_LLLL[i][j][k][l] / chi +
            ((-h[k][i] * h[j][l] + h[k][j] * h[i][l]) * dchi_dot_dchi +
             (-h[k][i] * d1.chi[j] * d1.chi[l] +
              h[k][j] * d1.chi[i] * d1.chi[l] +
              h[l][i] * d1.chi[j] * d1.chi[k] -
              h[l][j] * d1.chi[i] * d1.chi[k]) +
             2. * chi *
                 (-h[l][i] * covd_chi[k][j] + h[l][j] * covd_chi[k][i] +
                  h[k][i] * covd_chi[l][j] - h[k][j] * covd_chi[l][i])) /
                (4. * chi3);
    }

    m_riemann_spatial_LLLL = new const Tensor<4, data_t>(riemann_spatial_LLLL);
}
template_GQ void GeometricQuantities_t::compute_gauss_codazzi()
{
    if (m_gauss_codazzi != nullptr)
        delete m_gauss_codazzi;

    Tensor<4, data_t> gauss_codazzi;

    const auto &riemann_spatial_LLLL = get_riemann_spatial_LLLL();
    const auto &Kij = get_extrinsic_curvature();

    FOR(i, j, k, l)
    {
        gauss_codazzi[i][j][k][l] = riemann_spatial_LLLL[i][j][k][l] +
                                    Kij[i][k] * Kij[j][l] -
                                    Kij[i][l] * Kij[j][k];
    }

    m_gauss_codazzi = new const Tensor<4, data_t>(gauss_codazzi);
}
template_GQ void GeometricQuantities_t::compute_codazzi_mainardi()
{
    if (m_codazzi_mainardi != nullptr)
        delete m_codazzi_mainardi;

    Tensor<3, data_t> codazzi_mainardi;

    const auto &covd_Kij = get_covd_extrinsic_curvature();

    FOR(i, j, k)
    {
        codazzi_mainardi[i][j][k] = covd_Kij[i][k][j] - covd_Kij[j][k][i];
    }

    m_codazzi_mainardi = new const Tensor<3, data_t>(codazzi_mainardi);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_momentum_constraints()
{
    if (m_momentum_constraints != nullptr)
        delete m_momentum_constraints;

    const auto &d1 = get_d1_vars();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const auto &covD_K_tensor = get_covd_extrinsic_curvature();

    const Tensor<1, data_t> trace_covD_K_tensor =
        TensorAlgebra::compute_trace(covD_K_tensor, metric_UU_spatial, 0, 2);

    Tensor<1, data_t> M;
    FOR(i) { M[i] = trace_covD_K_tensor[i] - d1.K[i]; }

    if (m_em_tensor != nullptr) // if EM-tensor is set
    {
        const auto &em_tensor = get_em_tensor();
        FOR(i) { M[i] += -m_16_pi_G_Newton / 2. * em_tensor.Si[i]; }
    }

    m_momentum_constraints = new const Tensor<1, data_t>(M);
}
template_GQ void GeometricQuantities_t::compute_weyl_electric_part()
{
    if (m_weyl_electric_part != nullptr)
        delete m_weyl_electric_part;

    const auto &ricci =
        get_ricci_qDZ((m_formulation == CCZ4::USE_CCZ4 ? 1 : 0));
    const auto &Kij = get_extrinsic_curvature();
    const auto &h_UU_spatial = get_metric_UU_spatial();
    const auto &vars = get_vars();

    const Tensor<2, data_t> Kim_Kmj =
        TensorAlgebra::compute_dot_product(Kij, Kij, h_UU_spatial);

    data_t K_minus_Theta = vars.K;
    if (m_formulation == CCZ4::USE_CCZ4)
        K_minus_Theta -= vars.Theta;

    Tensor<2, data_t> weyl_electric_part;
    FOR(i, j)
    {
        weyl_electric_part[i][j] =
            ricci.LL[i][j] + K_minus_Theta * Kij[i][j] - Kim_Kmj[i][j];
    }

    // if EM-tensor is set (cosmological constant part is trace-free)
    if (m_em_tensor != nullptr)
    {
        const auto &em_tensor = get_em_tensor();
        FOR(i, j)
        {
            weyl_electric_part[i][j] +=
                -m_16_pi_G_Newton / 4. * em_tensor.Sij[i][j];
        }
    }

    const auto &h_UU = get_h_UU();
    TensorAlgebra::make_trace_free(weyl_electric_part, vars.h, h_UU);

    m_weyl_electric_part = new const Tensor<2, data_t>(weyl_electric_part);
}
template_GQ void GeometricQuantities_t::compute_lie_Z()
{
    if (m_lie_Z != nullptr)
        delete m_lie_Z;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &M = get_momentum_constraints();
    const auto &Z = get_Z();
    const auto &Kij_dot_Z = TensorAlgebra::compute_dot_product(
        get_extrinsic_curvature(), get_Z_U());

    data_t kappa1 = m_ccz4_params->kappa1;
    if (m_ccz4_params->covariantZ4)
        kappa1 /= vars.lapse;

    Tensor<1, data_t> lie_Z;
    FOR(i)
    {
        lie_Z[i] = M[i] - vars.Theta * d1.lapse[i] / vars.lapse + d1.Theta[i] -
                   kappa1 * Z[i] - 2. * Kij_dot_Z[i];
    }

    m_lie_Z = new const Tensor<1, data_t>(lie_Z);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_em_tensor_effective()
{
    if (m_em_tensor_effective != nullptr)
        delete m_em_tensor_effective;

    // 0 if EM-tensor not set -> just use cosmological constant
    emtensor_t<data_t> em_tensor_effective;

    if (m_em_tensor != nullptr)
        em_tensor_effective = get_em_tensor();
    else
    {
        em_tensor_effective.Sij = 0.;
        em_tensor_effective.Si = 0.;
        em_tensor_effective.S = 0.;
        em_tensor_effective.rho = 0.;
    }

    if (m_cosmological_constant != 0.)
    {

        const auto &vars = get_vars();
        const double two_cosmological_constant_over_16piG =
            2. * m_cosmological_constant / m_16_pi_G_Newton;

        FOR(i, j)
        {
            em_tensor_effective.Sij[i][j] -=
                two_cosmological_constant_over_16piG * vars.h[i][j];
        }
        em_tensor_effective.rho += two_cosmological_constant_over_16piG;
        em_tensor_effective.S -= 3. * two_cosmological_constant_over_16piG;
    }

    m_em_tensor_effective = new const emtensor_t<data_t>(em_tensor_effective);
}
template_GQ void GeometricQuantities_t::compute_hamiltonian_constraint()
{
    if (m_hamiltonian_constraint != nullptr)
        delete m_hamiltonian_constraint;

    const auto &vars = get_vars();
    const auto &h_UU = get_h_UU();
    const auto &ricci = get_ricci();

    const auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    const data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

    data_t H = ricci.scalar +
               (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;

    // if EM-tensor or cosmological constant are set
    if (m_em_tensor != nullptr || m_cosmological_constant != 0.)
        H -= m_16_pi_G_Newton * get_em_tensor_effective().rho;

    m_hamiltonian_constraint = new const data_t(H);
}
template_GQ void GeometricQuantities_t::compute_lie_derivatives()
{
    if (m_lie_derivatives != nullptr)
        delete m_lie_derivatives;

    Vars LIE;
    const auto &vars = get_vars();

    ///////////////////////////
    // RHS chi
    LIE.chi = (2.0 / GR_SPACEDIM) * vars.chi * vars.K;

    ///////////////////////////
    // RHS hij
    FOR(i, j) { LIE.h[i][j] = -2.0 * vars.A[i][j]; }

    ///////////////////////////
    // RHS Aij
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();
    const auto &chris_spatial = get_chris_spatial();

    const emtensor_t<data_t> *em_tensor_effective = nullptr;
    // if EM-tensor or cosmological constant are set
    if (m_em_tensor != nullptr || m_cosmological_constant != 0.)
        em_tensor_effective = &get_em_tensor_effective();

    const ricci_t<data_t> ricci = get_ricci_qDZ(2);
    const Tensor<2, data_t> covd2lapse = get_covd_lapse();

    const Tensor<2, data_t> A_LU =
        TensorAlgebra::compute_dot_product(vars.A, h_UU);
    const Tensor<2, data_t> Aik_Ajk =
        TensorAlgebra::compute_dot_product(vars.A, A_LU);
    data_t K_minus_2Theta = vars.K;
    if (m_formulation == CCZ4::USE_CCZ4)
        K_minus_2Theta -= 2. * vars.Theta;

    Tensor<2, data_t> LIE_A_TF_part;
    FOR(i, j)
    {
        LIE_A_TF_part[i][j] =
            vars.chi * (-covd2lapse[i][j] / vars.lapse + ricci.LL[i][j]);
    }
    // Matter terms
    if (m_em_tensor != nullptr) // if EM-tensor is set (cosmological constant
                                // does not matter for Aij RHS)
    {
        FOR(i, j)
        {
            LIE_A_TF_part[i][j] +=
                -vars.chi * m_16_pi_G_Newton / 2. * m_em_tensor->Sij[i][j];
        }
    }
    TensorAlgebra::make_trace_free(LIE_A_TF_part, vars.h, h_UU);

    FOR(i, j)
    {
        LIE.A[i][j] = LIE_A_TF_part[i][j] + vars.A[i][j] * K_minus_2Theta -
                      2. * Aik_Ajk[i][j];
    }

    ///////////////////////////
    ///////////////////////////
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const auto &chris = get_chris();

    const data_t tr_covd2lapse =
        TensorAlgebra::compute_trace(covd2lapse, metric_UU_spatial);
    const data_t tr_A2 = TensorAlgebra::compute_trace(Aik_Ajk, h_UU);

    data_t kappa1 = m_ccz4_params->kappa1;
    if (m_ccz4_params->covariantZ4)
        kappa1 /= vars.lapse;

    if (m_formulation == CCZ4::USE_BSSN)
    {
        ///////////////////////////
        // RHS K

        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        LIE.K = -tr_covd2lapse / vars.lapse +
                (tr_A2 + vars.K * vars.K / GR_SPACEDIM);

        // if EM-tensor or cosmological constant are set
        if (em_tensor_effective != nullptr)
            LIE.K += m_16_pi_G_Newton / (2. * (GR_SPACEDIM - 1.)) *
                     (em_tensor_effective->S +
                      (GR_SPACEDIM - 2.) * em_tensor_effective->rho);

        ///////////////////////////
        // RHS Theta
        LIE.Theta = 0.;
    }
    else
    {
        ///////////////////////////
        // RHS K
        LIE.K = -tr_covd2lapse / vars.lapse + ricci.scalar +
                vars.K * K_minus_2Theta -
                2. * GR_SPACEDIM / (GR_SPACEDIM - 1.) * kappa1 *
                    (1. + m_ccz4_params->kappa2) * vars.Theta;

        // if EM-tensor or cosmological constant are set
        if (em_tensor_effective != nullptr)
            LIE.K += m_16_pi_G_Newton / (2. * (GR_SPACEDIM - 1.)) *
                     (em_tensor_effective->S -
                      GR_SPACEDIM * em_tensor_effective->rho);

        ///////////////////////////
        // RHS Theta
        const auto &Z_U = get_Z_U();
        const data_t Z_U_dot_dlapse =
            TensorAlgebra::compute_dot_product(Z_U, d1.lapse);

        LIE.Theta =
            0.5 * (ricci.scalar - 2. * vars.K * vars.Theta +
                   (GR_SPACEDIM - 1.) / (double)GR_SPACEDIM * vars.K * vars.K -
                   tr_A2 - 2. / vars.lapse * Z_U_dot_dlapse -
                   2. * kappa1 * (2. + m_ccz4_params->kappa2) * vars.Theta);

        // if EM-tensor or cosmological constant are set
        if (em_tensor_effective != nullptr)
            LIE.Theta += -0.5 * m_16_pi_G_Newton * em_tensor_effective->rho;
    }

    ///////////////////////////
    // RHS Gamma

    const Tensor<2, data_t> A_UU =
        TensorAlgebra::compute_dot_product(A_LU, h_UU, 0, 0);
    const Tensor<1, data_t> A_dot_dlapse =
        TensorAlgebra::compute_dot_product(A_UU, d1.lapse);
    // \Gamma^i_^{jk} A^{jk}. - Note the abuse of the compute trace function.
    const Tensor<1, data_t> chris_dot_A =
        TensorAlgebra::compute_trace(chris.ULL, A_UU);
    const Tensor<1, data_t> A_dot_dchi =
        TensorAlgebra::compute_dot_product(A_UU, d1.chi);
    const Tensor<1, data_t> d1K_U =
        TensorAlgebra::compute_dot_product(d1.K, h_UU);
    const Tensor<1, data_t> d2shift_trace_U =
        TensorAlgebra::compute_dot_product(
            TensorAlgebra::compute_trace(d2.shift, 0, 2), h_UU);
    const Tensor<1, data_t> d2shift_laplacian =
        TensorAlgebra::compute_trace(d2.shift, h_UU);
    const Tensor<1, data_t> d1Theta_U =
        TensorAlgebra::compute_dot_product(d1.Theta, h_UU);
    const Tensor<1, data_t> d1lapse_U =
        TensorAlgebra::compute_dot_product(d1.lapse, h_UU);
    const data_t div_shift = TensorAlgebra::compute_trace(d1.shift);

    FOR(i)
    {
        LIE.Gamma[i] =
            -2. / vars.lapse * A_dot_dlapse[i] + 2. * chris_dot_A[i] -
            GR_SPACEDIM * A_dot_dchi[i] / vars.chi -
            2. * (GR_SPACEDIM - 1.) / (double)GR_SPACEDIM * d1K_U[i] +
            ((GR_SPACEDIM - 2.) / (double)GR_SPACEDIM * d2shift_trace_U[i] +
             d2shift_laplacian[i]) /
                vars.lapse;

        // if EM-tensor is set (cosmological
        // constant does not matter for Gamma^i RHS)
        if (m_em_tensor != nullptr)
        {
            const Tensor<1, data_t> Si_U = TensorAlgebra::compute_dot_product(
                em_tensor_effective->Si, h_UU);
            LIE.Gamma[i] += -m_16_pi_G_Newton * Si_U[i];
        }

        if (m_formulation == CCZ4::USE_CCZ4)
        {
            const auto &Z_U_conformal = get_Z_U_conformal();
            const Tensor<1, data_t> Z_U_conformal_dot_dshift =
                TensorAlgebra::compute_dot_product(Z_U_conformal, d1.shift, 0,
                                                   1);

            LIE.Gamma[i] +=
                2. * (d1Theta_U[i] - vars.Theta * d1lapse_U[i] / vars.lapse -
                      2. / GR_SPACEDIM * vars.K * Z_U_conformal[i]) -
                2. * kappa1 * Z_U_conformal[i];
            if (m_ccz4_params->kappa3 != 1.)
                LIE.Gamma[i] +=
                    2. * (m_ccz4_params->kappa3 - 1.) / vars.lapse *
                    (2. / GR_SPACEDIM * Z_U_conformal[i] * div_shift -
                     Z_U_conformal_dot_dshift[i]);
        }
    }

    LIE.lapse = 0.;
    LIE.shift = 0.;
    LIE.B = 0.;

    m_lie_derivatives = new const Vars(LIE);
}
template_GQ void GeometricQuantities_t::compute_lie_extrinsic_curvature()
{
    if (m_lie_extrinsic_curvature != nullptr)
        delete m_lie_extrinsic_curvature;

    Tensor<2, data_t> lie_extrinsic_curvature;

    const auto &vars = get_vars();
    const auto &LIE = get_lie_derivatives();
    const auto &Kij = get_extrinsic_curvature();

    FOR(i, j)
    {
        lie_extrinsic_curvature[i][j] =
            -LIE.chi / vars.chi * Kij[i][j] +
            (LIE.A[i][j] +
             (LIE.K * vars.h[i][j] + vars.K * LIE.h[i][j]) / GR_SPACEDIM) /
                vars.chi;
    }

    m_lie_extrinsic_curvature =
        new const Tensor<2, data_t>(lie_extrinsic_curvature);
}
template_GQ void GeometricQuantities_t::compute_eom_double_normal_projection()
{
    if (m_eom_double_normal_projection != nullptr)
        delete m_eom_double_normal_projection;

    Tensor<2, data_t> eom_double_normal_projection;

    const auto &vars = get_vars();
    const auto &LIE = get_lie_derivatives();
    const auto &Kij = get_extrinsic_curvature();
    const auto &h_UU_spatial = get_metric_UU_spatial();
    const auto &covd_lapse = get_covd_lapse();
    const auto &lie_extrinsic_curvature = get_lie_extrinsic_curvature();

    const Tensor<2, data_t> Kim_Kmj =
        TensorAlgebra::compute_dot_product(Kij, Kij, h_UU_spatial);

    FOR(i, j)
    {
        eom_double_normal_projection[i][j] = lie_extrinsic_curvature[i][j] +
                                             Kim_Kmj[i][j] +
                                             covd_lapse[i][j] / vars.lapse;
    }

    m_eom_double_normal_projection =
        new const Tensor<2, data_t>(eom_double_normal_projection);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_rhs_equations()
{
    if (m_rhs_equations != nullptr)
        delete m_rhs_equations;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &advec = get_advection();
    const auto &chris = get_chris();
    const auto &LIE = get_lie_derivatives();

    const data_t div_shift = TensorAlgebra::compute_trace(d1.shift);

    vars_t<data_t> rhs;
    rhs.chi = vars.lapse * LIE.chi +
              TensorAlgebra::lie_derivative(advec.chi, vars.chi, div_shift,
                                            -2. / GR_SPACEDIM);
    rhs.h = TensorAlgebra::lie_derivative(advec.h, vars.h, d1.shift, vars.shift,
                                          div_shift, -2. / GR_SPACEDIM);
    FOR(i, j) { rhs.h[i][j] += LIE.h[i][j] * vars.lapse; }

    rhs.K = vars.lapse * LIE.K +
            TensorAlgebra::lie_derivative(advec.K, vars.K, div_shift, 0.);

    rhs.A = TensorAlgebra::lie_derivative(advec.A, vars.A, d1.shift, vars.shift,
                                          div_shift, -2. / GR_SPACEDIM);
    FOR(i, j) { rhs.A[i][j] += LIE.A[i][j] * vars.lapse; }

    rhs.Theta =
        vars.lapse * LIE.Theta +
        TensorAlgebra::lie_derivative(advec.Theta, vars.Theta, div_shift, 0.);

    // Use the calculated christoffels in the lie derivative terms, not the
    // evolved ones, for BSSN (as according to the old code for BSSN and
    // according with Alcubierre page 87)
    // rhs.Gamma = TensorAlgebra::lie_derivative(advec.Gamma, vars.Gamma,
    // d1.shift, vars.shift, divshift,
    // 2. / GR_SPACEDIM, {true});
    rhs.Gamma = TensorAlgebra::lie_derivative(
        advec.Gamma,
        (m_formulation == CCZ4::USE_CCZ4 ? vars.Gamma : chris.contracted),
        d1.shift, vars.shift, div_shift, 2. / GR_SPACEDIM, {true});
    FOR(i) { rhs.Gamma[i] += LIE.Gamma[i] * vars.lapse; }

    // Gauge evolution equations
    rhs.lapse = m_ccz4_params->lapse_advec_coeff * advec.lapse -
                m_ccz4_params->lapse_coeff *
                    pow(vars.lapse, m_ccz4_params->lapse_power) *
                    (vars.K - 2 * vars.Theta);

    FOR(i)
    {
        rhs.shift[i] = m_ccz4_params->shift_advec_coeff * advec.shift[i] +
                       m_ccz4_params->shift_Gamma_coeff * vars.B[i];
        rhs.B[i] =
            rhs.Gamma[i] +
            m_ccz4_params->shift_advec_coeff * (advec.B[i] - advec.Gamma[i]) -
            m_ccz4_params->eta * vars.B[i];
    }
    // rhs.lapse = 0.;
    // rhs.shift = 0.;
    // rhs.B = 0.;

    m_rhs_equations = new const Vars(rhs);
}
/*
template_GQ void GeometricQuantities_t::compute_dt_chris_contracted()
{
    if (m_dt_chris_contracted != nullptr)
        delete m_dt_chris_contracted;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();
    const auto &rhs = get_rhs_equations();

    Tensor<3, data_t> d1_dt_h;
    FOR(i, j, l)
    {
        d1_dt_h[i][j][l] =
            -2. * (d1.lapse[l] * vars.A[i][j] + vars.lapse * d1.A[i][j][l]);
        FOR(k)
        {
            d1_dt_h[i][j][l] += (d1.shift[k][l] * d1.h[i][j][k] +
                                 vars.shift[k] * d2.h[i][j][k][l]) +
                                (d1.h[i][k][l] * d1.shift[k][j] +
                                 vars.h[i][k] * d2.shift[k][j][l]) +
                                (d1.h[k][j][l] * d1.shift[k][i] +
                                 vars.h[k][j] * d2.shift[k][i][l]) -
                                2. / 3. *
                                    (d1.h[i][j][l] * d1.shift[k][k] +
                                     vars.h[i][j] * d2.shift[k][k][l]);
        }
    }

    Tensor<1, data_t> dt_chris_contracted = 0.;
    FOR(i)
    {
        FOR(m, n, p)
        {
            data_t d1_terms = 0.;
            FOR(q, r)
            {
                d1_terms += -h_UU[q][r] * (rhs.h[n][q] * d1.h[m][p][r] +
                                           rhs.h[m][n] * d1.h[p][q][r]);
            }
            dt_chris_contracted[i] +=
                h_UU[i][m] * h_UU[n][p] * (d1_dt_h[m][n][p] + d1_terms);
        }
    }

    m_dt_chris_contracted = new const Tensor<1, data_t>(dt_chris_contracted);
}
template_GQ void GeometricQuantities_t::compute_dt_chris_spatial_contracted()
{
    if (m_dt_chris_spatial_contracted != nullptr)
        delete m_dt_chris_spatial_contracted;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();
    const auto &rhs = get_rhs_equations();
    const auto &dt_chris_contracted = get_dt_chris_contracted();
    const auto &chris = get_chris();

    Tensor<1, data_t> dtdichi;
    Tensor<1, data_t> dichi_U = 0.;
    FOR(i)
    {
        dtdichi[i] = (2. / GR_SPACEDIM *
                      (d1.chi[i] * vars.K * vars.lapse +
                       vars.chi * d1.K[i] * vars.lapse +
                       vars.chi * vars.K * d1.lapse[i]));
        FOR(j)
        {
            dtdichi[i] +=
                d1.shift[j][i] * d1.chi[j] + vars.shift[j] * d2.chi[j][i] -
                2. / GR_SPACEDIM *
                    (d1.chi[i] * d1.shift[j][j] + vars.chi * d2.shift[j][j][i]);

            dichi_U[i] += h_UU[i][j] * d1.chi[j];
        }
    }

    Tensor<1, data_t> dt_chris_spatial_contracted;
    FOR(i)
    {
        dt_chris_spatial_contracted[i] =
            rhs.chi * chris.contracted[i] + vars.chi * dt_chris_contracted[i];
        FOR(j)
        {
            dt_chris_spatial_contracted[i] +=
                (GR_SPACEDIM - 2.) / 2. * h_UU[i][j] * dtdichi[j];
            FOR(k)
            {
                dt_chris_spatial_contracted[i] +=
                    (GR_SPACEDIM - 2.) / 2. *
                    (-h_UU[i][k] * rhs.h[k][j] * dichi_U[j]);
            }
        }
    }

    m_dt_chris_spatial_contracted =
        new const Tensor<1, data_t>(dt_chris_spatial_contracted);
}
*/
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_metric_ST()
{
    if (m_metric_ST != nullptr)
        delete m_metric_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> g;

    const auto &vars = get_vars();
    const auto &metric_spatial = get_metric_spatial();
    const auto &shift_L = get_shift_L();

    const data_t shift2 =
        TensorAlgebra::compute_dot_product(vars.shift, shift_L);

    g[0][0] = -vars.lapse * vars.lapse + shift2;

    FOR(i)
    {
        g[0][i + 1] = shift_L[i];
        g[i + 1][0] = g[0][i + 1];
        FOR(j) { g[i + 1][j + 1] = metric_spatial[i][j]; }
    }
    m_metric_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(g);
}
template_GQ void GeometricQuantities_t::compute_projector_LU_ST()
{
    if (m_projector_LU_ST != nullptr)
        delete m_projector_LU_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> g_LU = 0.;

    const auto &vars = get_vars();

    FOR(i)
    {
        g_LU[0][i + 1] = vars.shift[i];
        g_LU[i + 1][i + 1] = 1.;
    }
    m_projector_LU_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(g_LU);
}
template_GQ void GeometricQuantities_t::compute_metric_UU_ST()
{
    if (m_metric_UU_ST != nullptr)
        delete m_metric_UU_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> g_UU;

    const auto &vars = get_vars();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const data_t lapse_squared = vars.lapse * vars.lapse;

    g_UU[0][0] = -1. / lapse_squared;

    FOR(i)
    {
        g_UU[0][i + 1] = vars.shift[i] / lapse_squared;
        g_UU[i + 1][0] = g_UU[0][i + 1];
        FOR(j)
        {
            g_UU[i + 1][j + 1] = metric_UU_spatial[i][j] -
                                 vars.shift[i] * vars.shift[j] / lapse_squared;
        }
    }
    m_metric_UU_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(g_UU);
}
template_GQ void GeometricQuantities_t::compute_normal_U_ST()
{
    if (m_normal_U_ST != nullptr)
        delete m_normal_U_ST;

    Tensor<1, data_t, CH_SPACEDIM + 1> normal_U_ST;

    const auto &vars = get_vars();

    normal_U_ST[0] = 1. / vars.lapse;

    FOR(i) { normal_U_ST[i + 1] = -vars.shift[i] / vars.lapse; }
    m_normal_U_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(normal_U_ST);
}
template_GQ void GeometricQuantities_t::compute_normal_L_ST()
{
    if (m_normal_L_ST != nullptr)
        delete m_normal_L_ST;

    Tensor<1, data_t, CH_SPACEDIM + 1> normal_L_ST = {0.};

    const auto &vars = get_vars();

    normal_L_ST[0] = -vars.lapse;

    m_normal_L_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(normal_L_ST);
}
template_GQ void GeometricQuantities_t::compute_shift_ST()
{
    if (m_shift_ST != nullptr)
        delete m_shift_ST;

    Tensor<1, data_t, CH_SPACEDIM + 1> shift_ST = {0.};

    const auto &vars = get_vars();
    FOR(i) { shift_ST[i + 1] = vars.shift[i]; }

    m_shift_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(shift_ST);
}
template_GQ void GeometricQuantities_t::compute_levi_civita_spatial_ST()
{
    if (m_levi_civita_spatial_ST != nullptr)
        delete m_levi_civita_spatial_ST;

    Tensor<3, data_t, CH_SPACEDIM + 1> levi_civita_spatial_ST = {0.};

    const auto &levi_civita_ST = get_levi_civita_ST();
    const auto &n_U = get_normal_U_ST();

    FOR_ST(i, j, k)
    {
        FOR_ST(l)
        {
            levi_civita_spatial_ST[i][j][k] +=
                n_U[l] * levi_civita_ST[l][i][j][k];
        }
    }

    m_levi_civita_spatial_ST =
        new Tensor<3, data_t, CH_SPACEDIM + 1>(levi_civita_spatial_ST);
}
template_GQ void GeometricQuantities_t::compute_levi_civita_ST()
{
    if (m_levi_civita_ST != nullptr)
        delete m_levi_civita_ST;

    Tensor<4, data_t, CH_SPACEDIM + 1> levi_civita_ST;

    const auto &vars = get_vars();

    const auto epsilon4_symbol = TensorAlgebra::epsilon4D();
    const data_t sqrt_g_det = vars.lapse * pow(vars.chi, -1.5);

    FOR_ST(i, j, k, l)
    {
        levi_civita_ST[i][j][k][l] = sqrt_g_det * epsilon4_symbol[i][j][k][l];
    }

    m_levi_civita_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>(levi_civita_ST);
}
template_GQ void GeometricQuantities_t::compute_Z_L_ST()
{
    if (m_Z_L_ST != nullptr)
        delete m_Z_L_ST;

    Tensor<1, data_t, CH_SPACEDIM + 1> Z_L_ST;

    const auto &vars = get_vars();
    const auto &Z = get_Z();

    Z_L_ST[0] = -vars.lapse * vars.Theta +
                TensorAlgebra::compute_dot_product(vars.shift, Z);
    FOR(i) { Z_L_ST[i + 1] = Z[i]; }

    m_Z_L_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(Z_L_ST);
}
template_GQ void GeometricQuantities_t::compute_grad_normal_LL()
{
    if (m_grad_normal_LL != nullptr)
        delete m_grad_normal_LL;

    Tensor<2, data_t, CH_SPACEDIM + 1> grad_normal_LL = 0.;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &Kij = get_extrinsic_curvature();
    const auto &shift_ST = get_shift_ST();
    const auto &normal_L_ST = get_normal_L_ST();
    const auto Kij_ST = TensorAlgebra::make_spatial_tensor_ST(Kij, shift_ST);
    const auto d1_lapse_ST =
        TensorAlgebra::make_spatial_tensor_ST(d1.lapse, shift_ST);

    FOR_ST(m, n)
    {
        grad_normal_LL[m][n] =
            -normal_L_ST[n] * d1_lapse_ST[m] / vars.lapse - Kij_ST[n][m];
    }

    m_grad_normal_LL = new Tensor<2, data_t, CH_SPACEDIM + 1>(grad_normal_LL);
}
template_GQ void GeometricQuantities_t::compute_covd_Z_L_ST()
{
    if (m_covd_Z_L_ST != nullptr)
        delete m_covd_Z_L_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> covd_Z_L_ST = 0.;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &LIE = get_lie_derivatives();
    const auto &lie_Z = get_lie_Z();
    const auto &Z_U = get_Z_U();
    const auto &Kij = get_extrinsic_curvature();
    const auto &chris_spatial = get_chris_spatial();
    const auto &covd_Z = get_covd_Z();
    const auto &shift_ST = get_shift_ST();
    const auto &normal_L_ST = get_normal_L_ST();

    // Call the projections of Z by Q (Theta), Qi, Qj, and Qij
    const data_t Q =
        -LIE.Theta -
        TensorAlgebra::compute_dot_product(Z_U, d1.lapse) / vars.lapse;

    const Tensor<1, data_t> Kij_dot_Z =
        TensorAlgebra::compute_dot_product(Kij, Z_U);

    Tensor<1, data_t> Qi;
    FOR(i)
    {
        Qi[i] = Kij_dot_Z[i] + vars.Theta * d1.lapse[i] / vars.lapse + lie_Z[i];
    }

    Tensor<1, data_t> Qj;
    FOR(i) { Qj[i] = Kij_dot_Z[i] - d1.Theta[i]; }

    Tensor<2, data_t> Qij;
    FOR(i, j) { Qij[i][j] = -Kij[i][j] * vars.Theta + covd_Z[i][j]; }

    const auto Qij_ST = TensorAlgebra::make_spatial_tensor_ST(Qij, shift_ST);
    const auto Qi_ST = TensorAlgebra::make_spatial_tensor_ST(Qi, shift_ST);
    const auto Qj_ST = TensorAlgebra::make_spatial_tensor_ST(Qj, shift_ST);

    FOR_ST(m, n)
    {
        covd_Z_L_ST[m][n] = Q * normal_L_ST[m] * normal_L_ST[n] -
                            normal_L_ST[m] * Qj_ST[n] -
                            normal_L_ST[n] * Qi_ST[m] + Qij_ST[m][n];
    }

    m_covd_Z_L_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(covd_Z_L_ST);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_em_tensor_ST()
{
    if (m_em_tensor_ST != nullptr)
        delete m_em_tensor_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> em_tensor_ST;

    const auto &em_tensor = get_em_tensor();
    const auto &normal_L_ST = get_normal_L_ST();
    const auto &shift_ST = get_shift_ST();

    const auto Smn_ST =
        TensorAlgebra::make_spatial_tensor_ST(em_tensor.Sij, shift_ST);
    const auto Sm_ST =
        TensorAlgebra::make_spatial_tensor_ST(em_tensor.Si, shift_ST);

    FOR_ST(m, n)
    {
        em_tensor_ST[m][n] = em_tensor.rho * normal_L_ST[m] * normal_L_ST[n] +
                             normal_L_ST[m] * Sm_ST[n] +
                             normal_L_ST[n] * Sm_ST[m] + Smn_ST[m][n];
    }
    m_em_tensor_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(em_tensor_ST);
}
template_GQ void GeometricQuantities_t::compute_em_tensor_trace_ST()
{
    if (m_em_tensor_trace_ST != nullptr)
        delete m_em_tensor_trace_ST;

    const auto &em_tensor = get_em_tensor();
    m_em_tensor_trace_ST = new data_t(em_tensor.S - em_tensor.rho);
}
template_GQ void GeometricQuantities_t::compute_weyl_tensor_LLLL()
{
    if (m_weyl_tensor_LLLL != nullptr)
        delete m_weyl_tensor_LLLL;

    const auto &g = get_metric_ST();
    const auto &g_UU = get_metric_UU_ST();
    const auto &n_L = get_normal_L_ST();
    const auto &levi_civita_spatial_ST = get_levi_civita_spatial_ST();
    const auto &shift_ST = get_shift_ST();
    const auto &Eij = get_weyl_electric_part();
    const auto &Bij = get_weyl_magnetic_part();

    // l_LL
    Tensor<2, data_t, 4> l_LL; // l[a][b] = g[a][b] + 2n[a]n[b]
    FOR_ST(a, b) l_LL[a][b] = g[a][b] + 2. * n_L[a] * n_L[b];

    Tensor<3, data_t, CH_SPACEDIM + 1> epsilon3_ULL = {0.};
    FOR_ST(a, b, c, d)
    {
        epsilon3_ULL[a][b][c] += g_UU[a][d] * levi_civita_spatial_ST[d][b][c];
    }

    const auto E_LL = TensorAlgebra::make_spatial_tensor_ST(Eij, shift_ST);
    const auto B_LL = TensorAlgebra::make_spatial_tensor_ST(Bij, shift_ST);

    // Weyl
    Tensor<4, data_t, CH_SPACEDIM + 1> weyl;

    FOR_ST(m, n, r, s)
    {
        weyl[m][n][r][s] = (l_LL[m][r] * E_LL[s][n] - l_LL[m][s] * E_LL[r][n]) -
                           (l_LL[n][r] * E_LL[s][m] - l_LL[n][s] * E_LL[r][m]);
        FOR_ST(b)
        weyl[m][n][r][s] += -(n_L[r] * B_LL[s][b] * epsilon3_ULL[b][m][n] -
                              n_L[s] * B_LL[r][b] * epsilon3_ULL[b][m][n]) -
                            (n_L[m] * B_LL[n][b] * epsilon3_ULL[b][r][s] -
                             n_L[n] * B_LL[m][b] * epsilon3_ULL[b][r][s]);
    }

    m_weyl_tensor_LLLL = new Tensor<4, data_t, CH_SPACEDIM + 1>(weyl);
}
template_GQ void GeometricQuantities_t::compute_weyl_squared()
{
    if (m_weyl_squared != nullptr)
        delete m_weyl_squared;

    const auto &metric_UU = get_metric_UU_spatial();
    const auto &Eij = get_weyl_electric_part();
    const auto &Bij = get_weyl_magnetic_part();

    data_t weyl_squared = 0.;
    FOR(i, j, k, l)
    {
        weyl_squared += 8. * metric_UU[i][k] * metric_UU[j][l] *
                        (Eij[i][j] * Eij[k][l] - Bij[i][j] * Bij[k][l]);
    }

    m_weyl_squared = new data_t(weyl_squared);
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_em_tensor_effective_ST()
{
    if (m_em_tensor_effective_ST != nullptr)
        delete m_em_tensor_effective_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> em_tensor_effective_ST;

    if (m_em_tensor != nullptr)
        em_tensor_effective_ST = get_em_tensor_ST();
    else
        em_tensor_effective_ST = 0.;

    if (m_cosmological_constant != 0.)
    {
        const auto &metric_ST = get_metric_ST();
        const double two_cosmological_constant_over_16piG =
            2. * m_cosmological_constant / m_16_pi_G_Newton;

        FOR_ST(m, n)
        {
            em_tensor_effective_ST[m][n] +=
                -two_cosmological_constant_over_16piG * metric_ST[m][n];
        }
    }

    m_em_tensor_effective_ST =
        new Tensor<2, data_t, CH_SPACEDIM + 1>(em_tensor_effective_ST);
}
template_GQ void GeometricQuantities_t::compute_em_tensor_effective_trace_ST()
{
    if (m_em_tensor_effective_trace_ST != nullptr)
        delete m_em_tensor_effective_trace_ST;

    m_em_tensor_effective_trace_ST = new data_t(
        get_em_tensor_trace_ST() -
        2. * (GR_SPACEDIM + 1.) * m_cosmological_constant / m_16_pi_G_Newton);
}
template_GQ void GeometricQuantities_t::compute_riemann_LLLL_ST()
{
    if (m_riemann_LLLL_ST != nullptr)
        delete m_riemann_LLLL_ST;

    const auto &weyl = get_weyl_tensor_LLLL();
    const auto &ricci = get_ricci_ST();
    const auto &ricci_scalar = get_ricci_scalar_ST();
    const auto &g = get_metric_ST();

    Tensor<4, data_t, CH_SPACEDIM + 1> riemann;

    FOR_ST(a, b, m, n)
    {
        riemann[a][b][m][n] = weyl[a][b][m][n] +
                              (g[a][m] * ricci[n][b] - g[a][n] * ricci[m][b] -
                               g[b][m] * ricci[n][a] + g[b][n] * ricci[m][a]) /
                                  (GR_SPACEDIM - 1.) -
                              (g[a][m] * g[n][b] - g[a][n] * g[m][b]) *
                                  ricci_scalar /
                                  (GR_SPACEDIM * (GR_SPACEDIM - 1.));
    }

    m_riemann_LLLL_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>(riemann);
}
template_GQ void GeometricQuantities_t::compute_riemann_LLLL_ST_v2()
{
    if (m_riemann_LLLL_ST_v2 != nullptr)
        delete m_riemann_LLLL_ST_v2;

    const auto &gauss_codazzi = get_gauss_codazzi();
    const auto &codazzi_mainardi = get_codazzi_mainardi();
    const auto &eom_double_normal_projection =
        get_eom_double_normal_projection();
    const auto &shift_ST = get_shift_ST();
    const auto &n_L = get_normal_L_ST();
    // const auto &projector_LU = get_projector_LU_ST();

    const auto GC_4D =
        TensorAlgebra::make_spatial_tensor_ST(gauss_codazzi, shift_ST);
    const auto CM_4D =
        TensorAlgebra::make_spatial_tensor_ST(codazzi_mainardi, shift_ST);
    const auto EOMnn_4D = TensorAlgebra::make_spatial_tensor_ST(
        eom_double_normal_projection, shift_ST);

    Tensor<4, data_t, CH_SPACEDIM + 1> riemann;
    FOR_ST(m, n, r, s)
    {
        riemann[m][n][r][s] =
            GC_4D[m][n][r][s] +
            (-n_L[s] * CM_4D[m][n][r] + n_L[r] * CM_4D[m][n][s] -
             n_L[n] * CM_4D[r][s][m] + n_L[m] * CM_4D[r][s][n]) +
            (n_L[s] * n_L[n] * EOMnn_4D[m][r] -
             n_L[r] * n_L[n] * EOMnn_4D[m][s] -
             n_L[s] * n_L[m] * EOMnn_4D[n][r] +
             n_L[r] * n_L[m] * EOMnn_4D[n][s]);
    }

    m_riemann_LLLL_ST_v2 = new Tensor<4, data_t, CH_SPACEDIM + 1>(riemann);
}
template_GQ void GeometricQuantities_t::compute_ricci_ST()
{
    if (m_ricci_ST != nullptr)
        delete m_ricci_ST;

    const auto &vars = get_vars();
    const auto &covd_Z_L_ST = get_covd_Z_L_ST();
    const auto &Z_L_ST = get_Z_L_ST();
    const auto &n_U = get_normal_U_ST();
    const auto &n_L = get_normal_L_ST();
    const auto &g = get_metric_ST();
    const auto &Tmn = get_em_tensor_effective_ST();
    const auto &T = get_em_tensor_effective_trace_ST();

    data_t kappa1 = m_ccz4_params->kappa1;
    if (m_ccz4_params->covariantZ4)
        kappa1 /= vars.lapse;

    const data_t Z_dot_n = TensorAlgebra::compute_dot_product(Z_L_ST, n_U);

    Tensor<2, data_t, CH_SPACEDIM + 1> ricci;
    FOR_ST(m, n)
    {
        ricci[m][n] = m_16_pi_G_Newton / 2. *
                      (Tmn[m][n] - g[m][n] * T / (GR_SPACEDIM - 1.));
        if (m_formulation == CCZ4::USE_CCZ4)
        {
            ricci[m][n] +=
                -covd_Z_L_ST[m][n] - covd_Z_L_ST[n][m] +
                kappa1 * (n_L[m] * Z_L_ST[n] + n_L[n] * Z_L_ST[m] -
                          2. / (GR_SPACEDIM - 1.) *
                              (1. + m_ccz4_params->kappa2) * g[m][n] * Z_dot_n);
        }
    }
    m_ricci_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(ricci);
}
template_GQ void GeometricQuantities_t::compute_ricci_scalar_ST()
{
    if (m_ricci_scalar_ST != nullptr)
        delete m_ricci_scalar_ST;

    const auto &ricci = get_ricci_ST();
    const auto &g_UU = get_metric_UU_ST();

    m_ricci_scalar_ST = new data_t(TensorAlgebra::compute_trace(ricci, g_UU));
}
template_GQ void GeometricQuantities_t::compute_ricci_squared()
{
    if (m_ricci_squared != nullptr)
        delete m_ricci_squared;

    const auto &g_UU = get_metric_UU_ST();
    const auto &ricci = get_ricci_ST();

    data_t ricci_squared = 0.;
    FOR_ST(a, b, c, d)
    {
        ricci_squared += g_UU[a][c] * g_UU[b][d] * ricci[a][b] * ricci[c][d];
    }

    m_ricci_squared = new data_t(ricci_squared);
}
template_GQ void GeometricQuantities_t::compute_kretschmann()
{
    if (m_kretschmann != nullptr)
        delete m_kretschmann;

    const auto &ricci_squared = get_ricci_squared();
    const auto &ricci_scalar = get_ricci_scalar_ST();
    const auto &weyl_squared = get_weyl_squared();

    m_kretschmann = new data_t(
        weyl_squared + 4. / (GR_SPACEDIM - 1.) * ricci_squared -
        2. / (GR_SPACEDIM * (GR_SPACEDIM - 1.)) * ricci_scalar * ricci_scalar);
}
template_GQ void GeometricQuantities_t::compute_riemann_squared()
{
    if (m_riemann_squared != nullptr)
        delete m_riemann_squared;

    const auto &riemann = get_riemann_LLLL_ST();
    const auto &g_UU = get_metric_UU_ST();

    data_t riemann_squared = 0.;
    FOR_ST(a, b, c, d)
    {
        FOR_ST(m, n, r, s)
        {
            riemann_squared += riemann[a][b][c][d] * riemann[m][n][r][s] *
                               g_UU[a][m] * g_UU[b][n] * g_UU[c][r] *
                               g_UU[d][s];
        }
    }

    m_riemann_squared = new data_t();
}
//////////////////////////////////////////////////////////////////////////
/*
template_GQ void GeometricQuantities_t::compute_chris_ST()
{
    if (m_chris_ST != nullptr)
        delete m_chris_ST;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &advec = get_advection();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const auto &chris_spatial = get_chris_spatial();
    const auto &Kij = get_extrinsic_curvature();
    const auto &rhs = get_rhs_equations();

    const Tensor<2, data_t> covd_shift = TensorAlgebra::covariant_derivative(
        d1.shift, vars.shift, chris_spatial, {true});
    const Tensor<1, data_t> shift_dot_Kij =
        TensorAlgebra::compute_dot_product(vars.shift, Kij);
    const data_t shift_dot_shift_dot_Kij =
        TensorAlgebra::compute_dot_product(vars.shift, shift_dot_Kij);

    Tensor<3, data_t, CH_SPACEDIM + 1> chris_ST;

    chris_ST[0][0][0] =
        (rhs.lapse + advec.lapse - shift_dot_shift_dot_Kij) / vars.lapse;

    Tensor<1, data_t> aux_term;
    FOR(i)
    {
        chris_ST[0][0][i + 1] = (d1.lapse[i] - shift_dot_Kij[i]) / vars.lapse;
        chris_ST[0][i + 1][0] = chris_ST[0][0][i + 1];

        FOR(j) { chris_ST[0][i + 1][j + 1] = -Kij[i][j] / vars.lapse; }

        aux_term[i] = d1.lapse[i] - 2. * shift_dot_Kij[i];
    }
    const Tensor<1, data_t> aux_term_U =
        TensorAlgebra::compute_dot_product(aux_term, metric_UU_spatial);

    const Tensor<1, data_t> shift_dot_covd_shift =
        TensorAlgebra::compute_dot_product(vars.shift, covd_shift);
    const Tensor<2, data_t> Kij_UL =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Kij, 0, 0);

    FOR(i)
    {
        chris_ST[i + 1][0][0] = vars.lapse * aux_term_U[i] -
                                vars.shift[i] * chris_ST[0][0][0] +
                                rhs.shift[i] + shift_dot_covd_shift[i];

        FOR(j)
        {
            chris_ST[i + 1][0][j + 1] = -vars.shift[i] * chris_ST[0][0][j + 1] -
                                        vars.lapse * Kij_UL[i][j] +
                                        covd_shift[i][j];
            chris_ST[i + 1][j + 1][0] = chris_ST[i + 1][0][j + 1];

            FOR(k)
            {
                chris_ST[i + 1][j + 1][k + 1] =
                    chris_spatial[i][j][k] +
                    vars.shift[i] * Kij[j][k] / vars.lapse;
            }
        }
    }

    m_chris_ST = new Tensor<3, data_t, CH_SPACEDIM + 1>(chris_ST);
}
template_GQ void GeometricQuantities_t::compute_d1_Z_L_ST()
{
    if (m_d1_Z_L_ST != nullptr)
        delete m_d1_Z_L_ST;

    Tensor<2, data_t, CH_SPACEDIM + 1> d1_Z_L_ST;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &Z = get_Z();
    const auto &chris = get_chris();
    const auto &d1_chris_contracted = get_d1_chris_contracted();
    const auto &dt_chris_contracted = get_dt_chris_contracted();
    const auto &rhs = get_rhs_equations();

    Tensor<1, data_t> dt_Z = 0.;
    Tensor<2, data_t> d1_Z = 0.;
    FOR(i, j)
    {
        dt_Z[i] += 0.5 * rhs.h[i][j] * (vars.Gamma[j] - chris.contracted[j]) +
                   0.5 * vars.h[i][j] * (rhs.Gamma[j] - dt_chris_contracted[j]);
        FOR(k)
        {
            d1_Z[i][k] +=
                0.5 * d1.h[i][j][k] * (vars.Gamma[j] - chris.contracted[j]) +
                0.5 * vars.h[i][j] *
                    (d1.Gamma[j][k] - d1_chris_contracted[j][k]);
        }
    }

    d1_Z_L_ST[0][0] = -vars.lapse * rhs.Theta - rhs.lapse * vars.Theta +
                      TensorAlgebra::compute_dot_product(rhs.shift, Z) +
                      TensorAlgebra::compute_dot_product(vars.shift, dt_Z);

    const auto &d1_shift_dot_Z =
        TensorAlgebra::compute_dot_product(Z, d1.shift, 0, 0);
    const auto &shift_dot_d1_Z =
        TensorAlgebra::compute_dot_product(vars.shift, d1_Z, 0, 0);

    FOR(i)
    {
        d1_Z_L_ST[i + 1][0] = dt_Z[i];
        d1_Z_L_ST[0][i + 1] = -vars.lapse * d1.Theta[i] -
                              d1.lapse[i] * vars.Theta + d1_shift_dot_Z[i] +
                              shift_dot_d1_Z[i];

        FOR(j) { d1_Z_L_ST[i + 1][j + 1] = d1_Z[i][j]; }
    }

    m_d1_Z_L_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>(d1_Z_L_ST);
}
*/

#endif /* GEOMETRICQUANTITIES_IMPL_HPP_ */
