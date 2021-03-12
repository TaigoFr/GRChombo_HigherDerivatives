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
              template <typename> class diff2_vars_t, class gauge_t>
#define GeometricQuantities_t                                                  \
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t>

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
    m_gauge = nullptr;
    m_em_tensor = nullptr;

    m_h_UU = nullptr;
    m_chris = nullptr;
    m_Z_U_conformal = nullptr;
    m_d1_chris_contracted = nullptr;
    m_covd_chi_conformal = nullptr;
    m_riemann_conformal_LLLL = nullptr;
    m_A_LU = nullptr;
    m_A_UU = nullptr;
    m_tr_A2 = nullptr;
    m_div_shift = nullptr;
    m_Gamma_L = nullptr;

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
    m_levi_civita_spatial_LUU = nullptr;
    m_covd_lapse = nullptr;
    m_ricci = nullptr;
    m_ricci_1DZ = nullptr;
    m_ricci_2DZ = nullptr;
    m_weyl_magnetic_part = nullptr;
    m_riemann_spatial_LLLL = nullptr;
    m_gauss_codazzi = nullptr;
    m_codazzi_mainardi = nullptr;
    m_Gamma_spatial = nullptr;
    m_Gamma_L_spatial = nullptr;
    m_acceleration_spatial = nullptr;

    m_momentum_constraints = nullptr;
    m_weyl_electric_part = nullptr;
    m_lie_Z = nullptr;

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
    m_acceleration_ST = nullptr;

    m_em_tensor_ST = nullptr;
    m_em_tensor_trace_ST = nullptr;
    m_weyl_tensor_LLLL = nullptr;
    m_weyl_squared = nullptr;

    m_riemann_LLLL_ST = nullptr;
    m_riemann_LLLL_ST_v2 = nullptr;
    m_ricci_ST = nullptr;
    m_ricci_scalar_ST = nullptr;
    m_ricci_squared = nullptr;
    m_kretschmann = nullptr;
    m_riemann_squared = nullptr;
    m_riemann_LLLU_ST = nullptr;
    m_riemann_LULU_ST = nullptr;

    m_chris_ST = nullptr;
    m_Gamma_ST = nullptr;
    m_Gamma_L_ST = nullptr;
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
    if (m_A_LU != nullptr)
    {
        delete m_A_LU;
        m_A_LU = nullptr;
    }
    if (m_A_UU != nullptr)
    {
        delete m_A_UU;
        m_A_UU = nullptr;
    }
    if (m_tr_A2 != nullptr)
    {
        delete m_tr_A2;
        m_tr_A2 = nullptr;
    }
    if (m_div_shift != nullptr)
    {
        delete m_div_shift;
        m_div_shift = nullptr;
    }
    if (m_Gamma_L != nullptr)
    {
        delete m_Gamma_L;
        m_Gamma_L = nullptr;
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
        delete m_covd_extrinsic_curvature;
        m_covd_extrinsic_curvature = nullptr;
    }
    if (m_levi_civita_spatial != nullptr)
    {
        delete m_levi_civita_spatial;
        m_levi_civita_spatial = nullptr;
    }
    if (m_levi_civita_spatial_LUU != nullptr)
    {
        delete m_levi_civita_spatial_LUU;
        m_levi_civita_spatial_LUU = nullptr;
    }
    if (m_covd_lapse != nullptr)
    {
        delete m_covd_lapse;
        m_covd_lapse = nullptr;
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
    if (m_Gamma_spatial != nullptr)
    {
        delete m_Gamma_spatial;
        m_Gamma_spatial = nullptr;
    }
    if (m_Gamma_L_spatial != nullptr)
    {
        delete m_Gamma_L_spatial;
        m_Gamma_L_spatial = nullptr;
    }
    if (m_acceleration_spatial != nullptr)
    {
        delete m_acceleration_spatial;
        m_acceleration_spatial = nullptr;
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
    if (m_acceleration_ST != nullptr)
    {
        delete m_acceleration_ST;
        m_acceleration_ST = nullptr;
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
    if (m_ricci != nullptr) // depe
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
    if (m_weyl_electric_part != nullptr)
    {
        delete m_weyl_electric_part;
        m_weyl_electric_part = nullptr;
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
    if (m_riemann_LLLU_ST != nullptr)
    {
        delete m_riemann_LLLU_ST;
        m_riemann_LLLU_ST = nullptr;
    }
    if (m_riemann_LULU_ST != nullptr)
    {
        delete m_riemann_LULU_ST;
        m_riemann_LULU_ST = nullptr;
    }

    clean_advection_and_gauge_dependent();
}

template_GQ void GeometricQuantities_t::clean_advection_and_gauge_dependent()
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
    */

    if (m_chris_ST != nullptr)
    {
        delete m_chris_ST;
        m_chris_ST = nullptr;
    }
    if (m_Gamma_ST != nullptr)
    {
        delete m_Gamma_ST;
        m_Gamma_ST = nullptr;
    }
    if (m_Gamma_L_ST != nullptr)
    {
        delete m_Gamma_L_ST;
        m_Gamma_L_ST = nullptr;
    }
    /*
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
template_GQ void
GeometricQuantities_t::set_advection_and_gauge(const Vars &a_advection,
                                               const gauge_t &a_gauge)
{
    m_advection = &a_advection;
    m_gauge = &a_gauge;
    clean_advection_and_gauge_dependent();
}
template_GQ void
GeometricQuantities_t::set_formulation(int formulation,
                                       const CCZ4_params_t<> &a_ccz4_params)
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

template_GQ int GeometricQuantities_t::get_formulation() const
{
    return m_formulation;
}
template_GQ const CCZ4_params_t<> &
GeometricQuantities_t::get_formulation_params() const
{
    return *m_ccz4_params;
}
template_GQ double GeometricQuantities_t::get_cosmological_constant() const
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
template_GQ const gauge_t &GeometricQuantities_t::get_gauge() const
{
    CH_assert(m_gauge != nullptr);
    return *m_gauge;
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
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
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
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_A_LU()
{
    if (m_A_LU == nullptr)
        compute_A_LU();
    return *m_A_LU;
}
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_A_UU()
{
    if (m_A_UU == nullptr)
        compute_A_UU();
    return *m_A_UU;
}
template_GQ const data_t &GeometricQuantities_t::get_tr_A2()
{
    if (m_tr_A2 == nullptr)
        compute_tr_A2();
    return *m_tr_A2;
}
template_GQ const data_t &GeometricQuantities_t::get_div_shift()
{
    if (m_div_shift == nullptr)
        compute_div_shift();
    return *m_div_shift;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Gamma_L()
{
    if (m_Gamma_L == nullptr)
        compute_Gamma_L();
    return *m_Gamma_L;
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
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
    if (m_Z_U == nullptr)
        compute_Z_U();
    return *m_Z_U;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Z()
{
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
    if (m_Z == nullptr)
        compute_Z();
    return *m_Z;
}
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_covd_Z()
{
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
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
template_GQ const Tensor<3, data_t> &
GeometricQuantities_t::get_levi_civita_spatial_LUU()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_levi_civita_spatial_LUU == nullptr)
        compute_levi_civita_spatial_LUU();
    return *m_levi_civita_spatial_LUU;
}
template_GQ const Tensor<2, data_t> &GeometricQuantities_t::get_covd_lapse()
{
    CH_assert(GR_SPACEDIM == 3);
    if (m_covd_lapse == nullptr)
        compute_covd_lapse();
    return *m_covd_lapse;
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
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_Gamma_spatial()
{
    if (m_Gamma_spatial == nullptr)
        compute_Gamma_spatial();
    return *m_Gamma_spatial;
}
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_Gamma_L_spatial()
{
    if (m_Gamma_L_spatial == nullptr)
        compute_Gamma_L_spatial();
    return *m_Gamma_L_spatial;
}
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_acceleration_spatial()
{
    if (m_acceleration_spatial == nullptr)
        compute_acceleration_spatial();
    return *m_acceleration_spatial;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<1, data_t> &
GeometricQuantities_t::get_momentum_constraints()
{
    if (m_momentum_constraints == nullptr)
        compute_momentum_constraints();
    return *m_momentum_constraints;
}
template_GQ const Tensor<1, data_t> &GeometricQuantities_t::get_lie_Z()
{
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_lie_Z == nullptr)
        compute_lie_Z();
    return *m_lie_Z;
}
//////////////////////////////////////////////////////////////////////////
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
GeometricQuantities_t::get_weyl_electric_part()
{
    CH_assert(GR_SPACEDIM == 3);
    CH_assert(m_formulation >= 0); // formulation is set
    if (m_weyl_electric_part == nullptr)
        compute_weyl_electric_part();
    return *m_weyl_electric_part;
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
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
    if (m_Z_L_ST == nullptr)
        compute_Z_L_ST();
    return *m_Z_L_ST;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_grad_normal_LL()
{
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
    if (m_grad_normal_LL == nullptr)
        compute_grad_normal_LL();
    return *m_grad_normal_LL;
}
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_covd_Z_L_ST()
{
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
    if (m_covd_Z_L_ST == nullptr)
        compute_covd_Z_L_ST();
    return *m_covd_Z_L_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_acceleration_ST()
{
    if (m_acceleration_ST == nullptr)
        compute_acceleration_ST();
    return *m_acceleration_ST;
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
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_riemann_LLLU_ST()
{
    CH_assert(m_formulation >= 0); // formulation is set
    CH_assert(GR_SPACEDIM == 3);
    if (m_riemann_LLLU_ST == nullptr)
        compute_riemann_LLLU_ST();
    return *m_riemann_LLLU_ST;
}
template_GQ const Tensor<4, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_riemann_LULU_ST()
{
    CH_assert(m_formulation >= 0); // formulation is set
    CH_assert(GR_SPACEDIM == 3);
    if (m_riemann_LULU_ST == nullptr)
        compute_riemann_LULU_ST();
    return *m_riemann_LULU_ST;
}
//////////////////////////////////////////////////////////////////////////
template_GQ const Tensor<3, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_chris_ST()
{
    if (m_chris_ST == nullptr)
        compute_chris_ST();
    return *m_chris_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_Gamma_ST()
{
    if (m_Gamma_ST == nullptr)
        compute_Gamma_ST();
    return *m_Gamma_ST;
}
template_GQ const Tensor<1, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_Gamma_L_ST()
{
    if (m_Gamma_L_ST == nullptr)
        compute_Gamma_L_ST();
    return *m_Gamma_L_ST;
}
/*
template_GQ const Tensor<2, data_t, CH_SPACEDIM + 1> &
GeometricQuantities_t::get_d1_Z_L_ST()
{
    CH_assert(m_formulation == CCZ4RHS<>::USE_CCZ4);
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
    m_h_UU =
        new Tensor<2, data_t>(TensorAlgebra::compute_inverse_sym(get_vars().h));
}
template_GQ void GeometricQuantities_t::compute_chris()
{
    if (m_chris != nullptr)
        delete m_chris;
    m_chris = new chris_t<data_t>(
        TensorAlgebra::compute_christoffel(get_d1_vars().h, get_h_UU()));
}
template_GQ void GeometricQuantities_t::compute_Z_U_conformal()
{
    if (m_Z_U_conformal != nullptr)
        delete m_Z_U_conformal;

    const auto &vars = get_vars();
    const auto &chris = get_chris();

    m_Z_U_conformal = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_Z_U_conformal)[i] = 0.5 * (vars.Gamma[i] - chris.contracted[i]);
    }
}
template_GQ void GeometricQuantities_t::compute_d1_chris_contracted()
{
    if (m_d1_chris_contracted != nullptr)
        delete m_d1_chris_contracted;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();

    m_d1_chris_contracted = new Tensor<2, data_t>({0.});
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
            (*m_d1_chris_contracted)[i][j] +=
                h_UU[i][m] * h_UU[n][p] * (d2.h[m][n][j][p] + d1_terms);
        }
    }
}
template_GQ void GeometricQuantities_t::compute_covd_chi_conformal()
{
    if (m_covd_chi_conformal != nullptr)
        delete m_covd_chi_conformal;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &chris = get_chris();

    m_covd_chi_conformal = new Tensor<2, data_t>(
        TensorAlgebra::covariant_derivative(d2.chi, d1.chi, chris.ULL));
}
template_GQ void GeometricQuantities_t::compute_riemann_conformal_LLLL()
{
    if (m_riemann_conformal_LLLL != nullptr)
        delete m_riemann_conformal_LLLL;

    m_riemann_conformal_LLLL = new Tensor<4, data_t>;

    const auto &chris = get_chris();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();

    FOR(i, j, k, l)
    {
        (*m_riemann_conformal_LLLL)[i][j][k][l] =
            0.5 * (d2.h[i][l][j][k] + d2.h[j][i][l][k] - d2.h[j][l][i][k]) -
            0.5 * (d2.h[i][k][j][l] + d2.h[j][i][k][l] - d2.h[j][k][i][l]);
        FOR(m)
        {
            (*m_riemann_conformal_LLLL)[i][j][k][l] +=
                chris.ULL[m][j][k] * d1.h[i][m][l] -
                chris.ULL[m][j][l] * d1.h[i][m][k] +
                chris.LLL[i][m][k] * chris.ULL[m][l][j] -
                chris.LLL[i][m][l] * chris.ULL[m][k][j];
        }
    }
}
template_GQ void GeometricQuantities_t::compute_A_LU()
{
    if (m_A_LU != nullptr)
        delete m_A_LU;

    m_A_LU = new Tensor<2, data_t>(
        TensorAlgebra::compute_dot_product(get_vars().A, get_h_UU()));
}
template_GQ void GeometricQuantities_t::compute_A_UU()
{
    if (m_A_UU != nullptr)
        delete m_A_UU;

    m_A_UU = new Tensor<2, data_t>(
        TensorAlgebra::compute_dot_product(get_A_LU(), get_h_UU(), 0, 0));
}
template_GQ void GeometricQuantities_t::compute_tr_A2()
{
    if (m_tr_A2 != nullptr)
        delete m_tr_A2;

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    m_tr_A2 =
        new data_t(TensorAlgebra::compute_trace(get_vars().A, get_A_UU()));
}
template_GQ void GeometricQuantities_t::compute_div_shift()
{
    if (m_div_shift != nullptr)
        delete m_div_shift;

    m_div_shift = new data_t(TensorAlgebra::compute_trace(get_d1_vars().shift));
}
template_GQ void GeometricQuantities_t::compute_Gamma_L()
{
    if (m_Gamma_L != nullptr)
        delete m_Gamma_L;
    m_Gamma_L = new Tensor<1, data_t>(TensorAlgebra::compute_dot_product(
        get_chris().contracted, get_vars().h));
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_metric_spatial()
{
    if (m_metric_spatial != nullptr)
        delete m_metric_spatial;

    const auto &vars = get_vars();

    m_metric_spatial = new Tensor<2, data_t>;
    FOR(i, j) { (*m_metric_spatial)[i][j] = vars.h[i][j] / vars.chi; }
}
template_GQ void GeometricQuantities_t::compute_metric_UU_spatial()
{
    if (m_metric_UU_spatial != nullptr)
        delete m_metric_UU_spatial;

    const auto &vars = get_vars();
    const auto &h_UU = get_h_UU();

    m_metric_UU_spatial = new Tensor<2, data_t>;
    FOR(i, j) { (*m_metric_UU_spatial)[i][j] = h_UU[i][j] * vars.chi; }
}
template_GQ void GeometricQuantities_t::compute_shift_L()
{
    if (m_shift_L != nullptr)
        delete m_shift_L;

    const auto &vars = get_vars();
    const auto &metric_spatial = get_metric_spatial();
    m_shift_L = new Tensor<1, data_t>(
        TensorAlgebra::lower_all(vars.shift, metric_spatial));
}
template_GQ void GeometricQuantities_t::compute_extrinsic_curvature()
{
    if (m_extrinsic_curvature != nullptr)
        delete m_extrinsic_curvature;

    const auto &vars = get_vars();

    m_extrinsic_curvature = new Tensor<2, data_t>;
    FOR(i, j)
    {
        (*m_extrinsic_curvature)[i][j] =
            (vars.A[i][j] + vars.K * vars.h[i][j] / GR_SPACEDIM) / vars.chi;
    }
}
template_GQ void GeometricQuantities_t::compute_chris_spatial()
{
    if (m_chris_spatial != nullptr)
        delete m_chris_spatial;

    const auto &chris = get_chris();
    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &h_UU = get_h_UU();

    m_chris_spatial = new Tensor<3, data_t>(TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris.ULL));
}
template_GQ void GeometricQuantities_t::compute_Z_U()
{
    if (m_Z_U != nullptr)
        delete m_Z_U;

    const auto &vars = get_vars();

    m_Z_U = new Tensor<1, data_t>(get_Z_U_conformal());
    FOR(i) { (*m_Z_U)[i] *= vars.chi; }
}
template_GQ void GeometricQuantities_t::compute_Z()
{
    if (m_Z != nullptr)
        delete m_Z;

    const auto &vars = get_vars();
    const auto &Z_U_conformal = get_Z_U_conformal();

    m_Z = new Tensor<1, data_t>(
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

    m_covd_Z = new Tensor<2, data_t>;
    FOR(i, j)
    {
        (*m_covd_Z)[i][j] =
            covd_Z_U_conformal_dot_h[i][j] +
            (Z[i] * d1.chi[j] + Z[j] * d1.chi[i] - vars.h[i][j] * Z_U_dot_chi) /
                (2. * vars.chi);
    }
}
template_GQ void GeometricQuantities_t::compute_d1_extrinsic_curvature()
{
    if (m_d1_extrinsic_curvature != nullptr)
        delete m_d1_extrinsic_curvature;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &K_tensor = get_extrinsic_curvature();

    m_d1_extrinsic_curvature = new Tensor<3, data_t>;
    FOR(i, j, k)
    {
        (*m_d1_extrinsic_curvature)[i][j][k] =
            (d1.A[i][j][k] - d1.chi[k] * K_tensor[i][j] +
             d1.h[i][j][k] * vars.K / GR_SPACEDIM +
             vars.h[i][j] * d1.K[k] / GR_SPACEDIM) /
            vars.chi;
    }
}
template_GQ void GeometricQuantities_t::compute_covd_extrinsic_curvature()
{
    if (m_covd_extrinsic_curvature != nullptr)
        delete m_covd_extrinsic_curvature;

    const auto &chris_spatial = get_chris_spatial();
    const auto &K_tensor = get_extrinsic_curvature();
    const auto &d1_K_tensor = get_d1_extrinsic_curvature();

    m_covd_extrinsic_curvature =
        new Tensor<3, data_t>(TensorAlgebra::covariant_derivative(
            d1_K_tensor, K_tensor, chris_spatial));
}
template_GQ void GeometricQuantities_t::compute_levi_civita_spatial()
{
    if (m_levi_civita_spatial != nullptr)
        delete m_levi_civita_spatial;

    const auto &levi_civita_ST = get_levi_civita_ST();
    const auto &n_U = get_normal_U_ST();

    m_levi_civita_spatial = new Tensor<3, data_t>({0.});
    FOR(i, j, k)
    {
        FOR_ST(l)
        {
            (*m_levi_civita_spatial)[i][j][k] +=
                n_U[l] * levi_civita_ST[l][i + 1][j + 1][k + 1];
        }
    }
}
template_GQ void GeometricQuantities_t::compute_levi_civita_spatial_LUU()
{
    if (m_levi_civita_spatial_LUU != nullptr)
        delete m_levi_civita_spatial_LUU;

    const auto &levi_civita_LLL = get_levi_civita_spatial();
    const auto &metric_UU = get_metric_UU_spatial();

    m_levi_civita_spatial_LUU = new Tensor<3, data_t>({0.});
    FOR3(i, j, k)
    {
        FOR2(m, n)
        {
            (*m_levi_civita_spatial_LUU)[i][j][k] +=
                levi_civita_LLL[i][m][n] * metric_UU[m][j] * metric_UU[n][k];
        }
    }
}
template_GQ void GeometricQuantities_t::compute_covd_lapse()
{
    if (m_covd_lapse != nullptr)
        delete m_covd_lapse;

    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &chris_spatial = get_chris_spatial();

    m_covd_lapse = new Tensor<2, data_t>(
        TensorAlgebra::covariant_derivative(d2.lapse, d1.lapse, chris_spatial));
}
template_GQ void GeometricQuantities_t::compute_weyl_magnetic_part()
{
    if (m_weyl_magnetic_part != nullptr)
        delete m_weyl_magnetic_part;

    const auto &levi_civita_spatial_LUU = get_levi_civita_spatial_LUU();
    const auto &covD_K_tensor = get_covd_extrinsic_curvature();

    m_weyl_magnetic_part = new Tensor<2, data_t>({0.});
    FOR4(i, j, k, l)
    {
        (*m_weyl_magnetic_part)[i][j] +=
            levi_civita_spatial_LUU[i][k][l] * covD_K_tensor[j][l][k];
    }

    TensorAlgebra::make_symmetric(*m_weyl_magnetic_part);
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

    m_riemann_spatial_LLLL = new Tensor<4, data_t>;
    FOR(i, j, k, l)
    {
        (*m_riemann_spatial_LLLL)[i][j][k][l] =
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
}
template_GQ void GeometricQuantities_t::compute_gauss_codazzi()
{
    if (m_gauss_codazzi != nullptr)
        delete m_gauss_codazzi;

    const auto &riemann_spatial_LLLL = get_riemann_spatial_LLLL();
    const auto &Kij = get_extrinsic_curvature();

    m_gauss_codazzi = new Tensor<4, data_t>;
    FOR(i, j, k, l)
    {
        (*m_gauss_codazzi)[i][j][k][l] = riemann_spatial_LLLL[i][j][k][l] +
                                         Kij[i][k] * Kij[j][l] -
                                         Kij[i][l] * Kij[j][k];
    }
}
template_GQ void GeometricQuantities_t::compute_codazzi_mainardi()
{
    if (m_codazzi_mainardi != nullptr)
        delete m_codazzi_mainardi;

    const auto &covd_Kij = get_covd_extrinsic_curvature();

    m_codazzi_mainardi = new Tensor<3, data_t>;
    FOR(i, j, k)
    {
        (*m_codazzi_mainardi)[i][j][k] = covd_Kij[i][k][j] - covd_Kij[j][k][i];
    }
}
template_GQ void GeometricQuantities_t::compute_Gamma_spatial()
{
    if (m_Gamma_spatial != nullptr)
        delete m_Gamma_spatial;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &h_UU = get_h_UU();
    const auto &chris = get_chris();

    m_Gamma_spatial = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_Gamma_spatial)[i] = chris.contracted[i] * vars.chi;
        FOR(j)
        {
            (*m_Gamma_spatial)[i] +=
                (GR_SPACEDIM - 2.) / 2. * h_UU[i][j] * d1.chi[j];
        }
    }
}
template_GQ void GeometricQuantities_t::compute_Gamma_L_spatial()
{
    if (m_Gamma_L_spatial != nullptr)
        delete m_Gamma_L_spatial;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &Gamma_L = get_Gamma_L();

    m_Gamma_L_spatial = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_Gamma_L_spatial)[i] =
            Gamma_L[i] + (GR_SPACEDIM - 2.) / 2. * d1.chi[i] / vars.chi;
    }
}
template_GQ void GeometricQuantities_t::compute_acceleration_spatial()
{
    if (m_acceleration_spatial != nullptr)
        delete m_acceleration_spatial;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();

    m_acceleration_spatial = new Tensor<1, data_t>;
    FOR(i) { (*m_acceleration_spatial)[i] = d1.lapse[i] / vars.lapse; }
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

    m_momentum_constraints = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_momentum_constraints)[i] =
            trace_covD_K_tensor[i] - d1.K[i] -
            (m_em_tensor == nullptr
                 ? 0.
                 : m_16_pi_G_Newton / 2. * m_em_tensor->Si[i]);
    }
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

    m_lie_Z = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_lie_Z)[i] = M[i] - vars.Theta * d1.lapse[i] / vars.lapse +
                        d1.Theta[i] - kappa1 * Z[i] - 2. * Kij_dot_Z[i];
    }
}
//////////////////////////////////////////////////////////////////////////
template_GQ ricci_t<data_t> GeometricQuantities_t::compute_ricci_qDZ(int q)
{
    CH_TIME("GeometricQuantities_t::compute_ricci_qDZ");

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
    if (q != 0 && m_formulation == CCZ4RHS<>::USE_CCZ4) // only if in CCZ4
        TensorAlgebra::hard_copy(Z_U_q, get_Z_U_conformal());

    /*
    const data_t boxtildechi =
        TensorAlgebra::compute_trace(covDtilde2chi, h_UU);
    const data_t dchi_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.chi, d1.chi, h_UU);
    */

    Tensor<3, data_t> chris_LLU = {0.};
    data_t boxtildechi = 0.;
    data_t dchi_dot_dchi = 0.;
    FOR(i, j)
    {
        boxtildechi += covDtilde2chi[i][j] * h_UU[i][j];
        dchi_dot_dchi += d1.chi[i] * d1.chi[j] * h_UU[i][j];
        FOR(k, l) { chris_LLU[i][j][k] += h_UU[k][l] * chris.LLL[i][j][l]; }
    }

    FOR(i, j)
    {
        data_t ricci_tilde = 0;
        data_t z_terms = 0.;
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

            if (q != 0 && m_formulation == CCZ4RHS<>::USE_CCZ4)
            {
                z_terms += Z_U_q[k] * (vars.h[i][k] * d1.chi[j] +
                                       vars.h[j][k] * d1.chi[i] -
                                       vars.h[i][j] * d1.chi[k]);
            }
        }

        const data_t ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covDtilde2chi[i][j] +
                   vars.h[i][j] * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1.chi[i] * d1.chi[j] +
                    GR_SPACEDIM * vars.h[i][j] * dchi_dot_dchi) /
                       (2 * vars.chi));

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

    m_ricci = new ricci_t<data_t>(compute_ricci_qDZ(0));
}
template_GQ void GeometricQuantities_t::compute_ricci_1DZ()
{
    if (m_ricci_1DZ != nullptr)
        delete m_ricci_1DZ;

    m_ricci_1DZ = new ricci_t<data_t>(compute_ricci_qDZ(1));
}
template_GQ void GeometricQuantities_t::compute_ricci_2DZ()
{
    if (m_ricci_2DZ != nullptr)
        delete m_ricci_2DZ;

    m_ricci_2DZ = new ricci_t<data_t>(compute_ricci_qDZ(2));
}
template_GQ void GeometricQuantities_t::compute_weyl_electric_part()
{
    if (m_weyl_electric_part != nullptr)
        delete m_weyl_electric_part;

    const auto &ricci =
        get_ricci_qDZ((m_formulation == CCZ4RHS<>::USE_CCZ4 ? 1 : 0));
    const auto &Kij = get_extrinsic_curvature();
    const auto &h_UU_spatial = get_metric_UU_spatial();
    const auto &vars = get_vars();

    const Tensor<2, data_t> Kim_Kmj =
        TensorAlgebra::compute_dot_product(Kij, Kij, h_UU_spatial);

    data_t K_minus_Theta = vars.K;
    if (m_formulation == CCZ4RHS<>::USE_CCZ4)
        K_minus_Theta -= vars.Theta;

    m_weyl_electric_part = new Tensor<2, data_t>;
    FOR(i, j)
    {
        (*m_weyl_electric_part)[i][j] =
            ricci.LL[i][j] + K_minus_Theta * Kij[i][j] - Kim_Kmj[i][j] -
            (m_em_tensor == nullptr
                 ? 0.
                 : m_16_pi_G_Newton / 4. * m_em_tensor->Sij[i][j]);
    }

    TensorAlgebra::make_trace_free(*m_weyl_electric_part, vars.h, get_h_UU());
}
template_GQ void GeometricQuantities_t::compute_hamiltonian_constraint()
{
    if (m_hamiltonian_constraint != nullptr)
        delete m_hamiltonian_constraint;

    const auto &vars = get_vars();
    const auto &ricci = get_ricci();
    const data_t &tr_A2 = get_tr_A2();

    m_hamiltonian_constraint = new data_t;
    (*m_hamiltonian_constraint) =
        ricci.scalar + (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM -
        tr_A2 - 2. * m_cosmological_constant -
        (m_em_tensor == nullptr ? 0. : m_16_pi_G_Newton * m_em_tensor->rho);
}
template_GQ void GeometricQuantities_t::compute_lie_derivatives()
{
    CH_TIME("GeometricQuantities::compute_lie_derivatives");

    if (m_lie_derivatives != nullptr)
        delete m_lie_derivatives;

    m_lie_derivatives = new Vars;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &h_UU = get_h_UU();
    const auto &A_LU = get_A_LU();
    const auto &chris = get_chris();
    const Tensor<2, data_t> &A_UU = get_A_UU();
    const data_t div_shift = get_div_shift();

    const ricci_t<data_t> ricci = get_ricci_qDZ(2);
    const Tensor<2, data_t> covd2lapse =
        get_covd_lapse(); // calculates chris_spatial -> not very efficient

    bool ccz4 = (m_formulation == CCZ4RHS<>::USE_CCZ4);
    bool k3 = (ccz4 && m_ccz4_params->kappa3 != 1.);
    bool vacuum = (m_em_tensor == nullptr);

    data_t kappa1 =
        m_ccz4_params->kappa1 / (m_ccz4_params->covariantZ4 ? vars.lapse : 1.);

    data_t K_minus_2Theta = vars.K - (ccz4 ? 2. * vars.Theta : 0.);

    Tensor<2, data_t> Aik_Ajk = 0.;
    data_t tr_covd2lapse = 0.;
    Tensor<2, data_t> LIE_A_TF_part;

    // LIE chi
    m_lie_derivatives->chi = (2.0 / GR_SPACEDIM) * vars.chi * vars.K;

    FOR(i, j)
    {
        // LIE hij
        m_lie_derivatives->h[i][j] = -2.0 * vars.A[i][j];

        FOR(k) { Aik_Ajk[i][j] += vars.A[i][k] * A_LU[j][k]; }

        tr_covd2lapse += covd2lapse[i][j] * h_UU[i][j] * vars.chi;

        LIE_A_TF_part[i][j] =
            vars.chi * (-covd2lapse[i][j] / vars.lapse + ricci.LL[i][j]) -
            (vacuum
                 ? 0.
                 : vars.chi * m_16_pi_G_Newton / 2. * m_em_tensor->Sij[i][j]);
    }

    TensorAlgebra::make_trace_free(LIE_A_TF_part, vars.h, h_UU);

    Tensor<1, data_t> Z_U_conformal = 0.;
    if (ccz4)
        Z_U_conformal = get_Z_U_conformal();

    data_t tr_A2 = 0.;
    FOR(i)
    {
        // LIE Gamma
        m_lie_derivatives->Gamma[i] =
            (-4. / GR_SPACEDIM * vars.K - 2. * kappa1 +
             (k3 ? 4. / GR_SPACEDIM * (m_ccz4_params->kappa3 - 1.) /
                       vars.lapse * div_shift
                 : 0.)) *
            Z_U_conformal[i];

        FOR(j)
        {
            // LIE Gamma
            m_lie_derivatives->Gamma[i] +=
                -(2. / vars.lapse * d1.lapse[j] +
                  GR_SPACEDIM * d1.chi[j] / vars.chi) *
                    A_UU[i][j] +
                2. *
                    (-(GR_SPACEDIM - 1.) / (double)GR_SPACEDIM * d1.K[j] +
                     d1.Theta[j] - vars.Theta * d1.lapse[j] / vars.lapse -
                     (vacuum ? 0.
                             : 0.5 * m_16_pi_G_Newton * m_em_tensor->Si[j])) *
                    h_UU[i][j] -
                (k3 ? 2. * (m_ccz4_params->kappa3 - 1.) / vars.lapse *
                          Z_U_conformal[j] * d1.shift[i][j]
                    : 0.);

            FOR(k)
            {
                // LIE Gamma
                m_lie_derivatives->Gamma[i] +=
                    2. * chris.ULL[i][j][k] * A_UU[j][k] +
                    ((GR_SPACEDIM - 2.) / (double)GR_SPACEDIM *
                         d2.shift[k][j][k] * h_UU[j][i] +
                     d2.shift[i][j][k] * h_UU[j][k]) /
                        vars.lapse;
            }

            // LIE Aij
            m_lie_derivatives->A[i][j] = LIE_A_TF_part[i][j] +
                                         vars.A[i][j] * K_minus_2Theta -
                                         2. * Aik_Ajk[i][j];

            tr_A2 += Aik_Ajk[i][j] * h_UU[i][j];
        }
    }

    if (ccz4)
    {
        // LIE K
        m_lie_derivatives->K =
            -tr_covd2lapse / vars.lapse + ricci.scalar +
            vars.K * K_minus_2Theta -
            2. * GR_SPACEDIM / (GR_SPACEDIM - 1.) * kappa1 *
                (1. + m_ccz4_params->kappa2) * vars.Theta -
            2. * GR_SPACEDIM / (GR_SPACEDIM - 1.) * m_cosmological_constant +
            (vacuum ? 0.
                    : m_16_pi_G_Newton / (2. * (GR_SPACEDIM - 1.)) *
                          (m_em_tensor->S - GR_SPACEDIM * m_em_tensor->rho));

        const auto &Z_U = get_Z_U();
        const data_t Z_U_dot_dlapse =
            TensorAlgebra::compute_dot_product(Z_U, d1.lapse);

        // LIE Theta
        m_lie_derivatives->Theta =
            0.5 * (ricci.scalar - 2. * vars.K * vars.Theta +
                   (GR_SPACEDIM - 1.) / (double)GR_SPACEDIM * vars.K * vars.K -
                   tr_A2 - 2. / vars.lapse * Z_U_dot_dlapse -
                   2. * kappa1 * (2. + m_ccz4_params->kappa2) * vars.Theta) -
            m_cosmological_constant -
            (vacuum ? 0. : 0.5 * m_16_pi_G_Newton * m_em_tensor->rho);
    }
    else
    {
        // LIE K
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        m_lie_derivatives->K =
            -tr_covd2lapse / vars.lapse +
            (tr_A2 + vars.K * vars.K / GR_SPACEDIM) -
            2. / (GR_SPACEDIM - 1.) * m_cosmological_constant +
            (vacuum ? 0.
                    : m_16_pi_G_Newton / (2. * (GR_SPACEDIM - 1.)) *
                          (m_em_tensor->S +
                           (GR_SPACEDIM - 2.) * m_em_tensor->rho));

        // LIE Theta
        m_lie_derivatives->Theta = 0.;
    }

    m_lie_derivatives->lapse = 0.;
    m_lie_derivatives->shift = 0.;
    m_lie_derivatives->B = 0.;
}
template_GQ void GeometricQuantities_t::compute_lie_extrinsic_curvature()
{
    if (m_lie_extrinsic_curvature != nullptr)
        delete m_lie_extrinsic_curvature;

    const auto &vars = get_vars();
    const auto &LIE = get_lie_derivatives();
    const auto &Kij = get_extrinsic_curvature();

    m_lie_extrinsic_curvature = new Tensor<2, data_t>;
    FOR(i, j)
    {
        (*m_lie_extrinsic_curvature)[i][j] =
            -LIE.chi / vars.chi * Kij[i][j] +
            (LIE.A[i][j] +
             (LIE.K * vars.h[i][j] + vars.K * LIE.h[i][j]) / GR_SPACEDIM) /
                vars.chi;
    }
}
template_GQ void GeometricQuantities_t::compute_eom_double_normal_projection()
{
    if (m_eom_double_normal_projection != nullptr)
        delete m_eom_double_normal_projection;

    const auto &vars = get_vars();
    const auto &Kij = get_extrinsic_curvature();
    const auto &h_UU_spatial = get_metric_UU_spatial();
    const auto &covd_lapse = get_covd_lapse();
    const auto &lie_extrinsic_curvature = get_lie_extrinsic_curvature();

    const Tensor<2, data_t> Kim_Kmj =
        TensorAlgebra::compute_dot_product(Kij, Kij, h_UU_spatial);

    m_eom_double_normal_projection = new Tensor<2, data_t>;
    FOR(i, j)
    {
        (*m_eom_double_normal_projection)[i][j] =
            lie_extrinsic_curvature[i][j] + Kim_Kmj[i][j] +
            covd_lapse[i][j] / vars.lapse;
    }
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_rhs_equations()
{
    CH_TIME("GeometricQuantities::compute_rhs_equations");

    if (m_rhs_equations != nullptr)
        delete m_rhs_equations;

    m_rhs_equations = new Vars;
    compute_rhs_equations(*m_rhs_equations);
}
template_GQ void GeometricQuantities_t::compute_rhs_equations(Vars &rhs)
{
    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &d2 = get_d2_vars();
    const auto &advec = get_advection();
    const auto &LIE = get_lie_derivatives();
    const auto &chris = get_chris();
    const auto &gauge = get_gauge();

    const data_t &div_shift = get_div_shift();

    {
        CH_TIME("GeometricQuantities::compute_rhs_equations::RHS");

        rhs.chi = vars.lapse * LIE.chi + advec.chi -
                  2. / GR_SPACEDIM * vars.chi * div_shift;
        rhs.K = vars.lapse * LIE.K + advec.K;
        rhs.Theta = vars.lapse * LIE.Theta + advec.Theta;

        FOR(i)
        {
            // Use the calculated christoffels in the lie derivative terms, not
            // the evolved ones, for BSSN (as according to the old code for BSSN
            // and according with Alcubierre page 87)
            rhs.Gamma[i] = advec.Gamma[i] + LIE.Gamma[i] * vars.lapse +
                           2. / GR_SPACEDIM *
                               (m_formulation == CCZ4RHS<>::USE_CCZ4
                                    ? vars.Gamma[i]
                                    : chris.contracted[i]) *
                               div_shift;

            FOR1(j)
            {
                // Use the calculated christoffels in the lie derivative terms,
                // not the evolved ones, for BSSN (as according to the old code
                // for BSSN and according with Alcubierre page 87)
                rhs.Gamma[i] -= (m_formulation == CCZ4RHS<>::USE_CCZ4
                                     ? vars.Gamma[j]
                                     : chris.contracted[j]) *
                                d1.shift[i][j];

                rhs.h[i][j] = advec.h[i][j] + LIE.h[i][j] * vars.lapse -
                              2. / GR_SPACEDIM * vars.h[i][j] * div_shift;
                rhs.A[i][j] = advec.A[i][j] + LIE.A[i][j] * vars.lapse -
                              2. / GR_SPACEDIM * vars.A[i][j] * div_shift;

                FOR1(k)
                {
                    rhs.h[i][j] += vars.h[i][k] * d1.shift[k][j] +
                                   vars.h[j][k] * d1.shift[k][i];
                    rhs.A[i][j] += vars.A[i][k] * d1.shift[k][j] +
                                   vars.A[j][k] * d1.shift[k][i];
                }
            }
        }

        // Gauge evolution equations
        gauge.rhs_gauge(rhs, vars, d1, d2, advec);

        /*
        // Alternative using TensorAlgebra::lie_derivative (is correct, but
        slower)

        rhs.chi =
            vars.lapse * LIE.chi +
            TensorAlgebra::lie_derivative(advec.chi, vars.chi, div_shift,
                                          -2. / GR_SPACEDIM);
        rhs.h =
            TensorAlgebra::lie_derivative(advec.h, vars.h, d1.shift, vars.shift,
                                          div_shift, -2. / GR_SPACEDIM);
        FOR(i, j) { rhs.h[i][j] += LIE.h[i][j] * vars.lapse; }

        rhs.K =
            vars.lapse * LIE.K +
            TensorAlgebra::lie_derivative(advec.K, vars.K, div_shift, 0.);

        rhs.A =
            TensorAlgebra::lie_derivative(advec.A, vars.A, d1.shift, vars.shift,
                                          div_shift, -2. / GR_SPACEDIM);
        FOR(i, j) { rhs.A[i][j] += LIE.A[i][j] * vars.lapse; }

        rhs.Theta = vars.lapse * LIE.Theta +
                                 TensorAlgebra::lie_derivative(
                                     advec.Theta, vars.Theta, div_shift, 0.);

        // Use the calculated christoffels in the lie derivative terms, not the
        // evolved ones, for BSSN (as according to the old code for BSSN and
        // according with Alcubierre page 87)
        // rhs.Gamma = TensorAlgebra::lie_derivative(advec.Gamma,
        // vars.Gamma, d1.shift, vars.shift, div_shift,
        // 2. / GR_SPACEDIM, {true});
        rhs.Gamma = TensorAlgebra::lie_derivative(
            advec.Gamma,
            (m_formulation == CCZ4RHS<>::USE_CCZ4 ? vars.Gamma :
        chris.contracted), d1.shift, vars.shift, div_shift, 2. / GR_SPACEDIM,
        {true}); FOR(i) { rhs.Gamma[i] += LIE.Gamma[i] * vars.lapse; }

        // Gauge evolution equations
        rhs.lapse =
            m_ccz4_params->lapse_advec_coeff * advec.lapse -
            m_ccz4_params->lapse_coeff *
                pow(vars.lapse, m_ccz4_params->lapse_power) *
                (vars.K - 2 * vars.Theta);

        FOR(i)
        {
            rhs.shift[i] =
                m_ccz4_params->shift_advec_coeff * advec.shift[i] +
                m_ccz4_params->shift_Gamma_coeff * vars.B[i];
            rhs.B[i] = rhs.Gamma[i] +
                                    m_ccz4_params->shift_advec_coeff *
                                        (advec.B[i] - advec.Gamma[i]) -
                                    m_ccz4_params->eta * vars.B[i];
        }
        // rhs.lapse = 0.;
        // rhs.shift = 0.;
        // rhs.B = 0.;

        */
    }
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

    m_dt_chris_contracted = new const Tensor<1, data_t>({0.});
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
            (*m_dt_chris_contracted)[i] +=
                h_UU[i][m] * h_UU[n][p] * (d1_dt_h[m][n][p] + d1_terms);
        }
    }

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

    m_dt_chris_spatial_contracted = new Tensor<1, data_t>;
    FOR(i)
    {
        (*m_dt_chris_spatial_contracted)[i] =
            rhs.chi * chris.contracted[i] + vars.chi * dt_chris_contracted[i];
        FOR(j)
        {
            (*m_dt_chris_spatial_contracted)[i] +=
                (GR_SPACEDIM - 2.) / 2. * h_UU[i][j] * dtdichi[j];
            FOR(k)
            {
                (*m_dt_chris_spatial_contracted)[i] +=
                    (GR_SPACEDIM - 2.) / 2. *
                    (-h_UU[i][k] * rhs.h[k][j] * dichi_U[j]);
            }
        }
    }

}
*/
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_metric_ST()
{
    if (m_metric_ST != nullptr)
        delete m_metric_ST;

    const auto &vars = get_vars();
    const auto &metric_spatial = get_metric_spatial();
    const auto &shift_L = get_shift_L();

    const data_t shift2 =
        TensorAlgebra::compute_dot_product(vars.shift, shift_L);

    m_metric_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    (*m_metric_ST)[0][0] = -vars.lapse * vars.lapse + shift2;

    FOR(i)
    {
        (*m_metric_ST)[0][i + 1] = shift_L[i];
        (*m_metric_ST)[i + 1][0] = (*m_metric_ST)[0][i + 1];
        FOR(j) { (*m_metric_ST)[i + 1][j + 1] = metric_spatial[i][j]; }
    }
}
template_GQ void GeometricQuantities_t::compute_projector_LU_ST()
{
    if (m_projector_LU_ST != nullptr)
        delete m_projector_LU_ST;

    const auto &vars = get_vars();

    m_projector_LU_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>({0.});
    FOR(i)
    {
        (*m_projector_LU_ST)[0][i + 1] = vars.shift[i];
        (*m_projector_LU_ST)[i + 1][i + 1] = 1.;
    }
}
template_GQ void GeometricQuantities_t::compute_metric_UU_ST()
{
    if (m_metric_UU_ST != nullptr)
        delete m_metric_UU_ST;

    const auto &vars = get_vars();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const data_t lapse_squared = vars.lapse * vars.lapse;

    m_metric_UU_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    (*m_metric_UU_ST)[0][0] = -1. / lapse_squared;

    FOR(i)
    {
        (*m_metric_UU_ST)[0][i + 1] = vars.shift[i] / lapse_squared;
        (*m_metric_UU_ST)[i + 1][0] = (*m_metric_UU_ST)[0][i + 1];
        FOR(j)
        {
            (*m_metric_UU_ST)[i + 1][j + 1] =
                metric_UU_spatial[i][j] -
                vars.shift[i] * vars.shift[j] / lapse_squared;
        }
    }
}
template_GQ void GeometricQuantities_t::compute_normal_U_ST()
{
    if (m_normal_U_ST != nullptr)
        delete m_normal_U_ST;

    const auto &vars = get_vars();

    m_normal_U_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>;
    (*m_normal_U_ST)[0] = 1. / vars.lapse;

    FOR(i) { (*m_normal_U_ST)[i + 1] = -vars.shift[i] / vars.lapse; }
}
template_GQ void GeometricQuantities_t::compute_normal_L_ST()
{
    if (m_normal_L_ST != nullptr)
        delete m_normal_L_ST;

    m_normal_L_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>({0.});
    (*m_normal_L_ST)[0] = -get_vars().lapse;
}
template_GQ void GeometricQuantities_t::compute_shift_ST()
{
    if (m_shift_ST != nullptr)
        delete m_shift_ST;

    const auto &vars = get_vars();

    m_shift_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>({0.});
    FOR(i) { (*m_shift_ST)[i + 1] = vars.shift[i]; }
}
template_GQ void GeometricQuantities_t::compute_levi_civita_spatial_ST()
{
    if (m_levi_civita_spatial_ST != nullptr)
        delete m_levi_civita_spatial_ST;

    const auto &levi_civita_ST = get_levi_civita_ST();
    const auto &n_U = get_normal_U_ST();

    m_levi_civita_spatial_ST = new Tensor<3, data_t, CH_SPACEDIM + 1>({0.});
    FOR_ST(i, j, k)
    {
        FOR_ST(l)
        {
            (*m_levi_civita_spatial_ST)[i][j][k] +=
                n_U[l] * levi_civita_ST[l][i][j][k];
        }
    }
}
template_GQ void GeometricQuantities_t::compute_levi_civita_ST()
{
    if (m_levi_civita_ST != nullptr)
        delete m_levi_civita_ST;

    Tensor<4, data_t, CH_SPACEDIM + 1> levi_civita_ST;

    const auto &vars = get_vars();

    const auto epsilon4_symbol = TensorAlgebra::epsilon4D();
    const data_t sqrt_g_det = vars.lapse * pow(vars.chi, -1.5);

    m_levi_civita_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>;
    FOR_ST(i, j, k, l)
    {
        (*m_levi_civita_ST)[i][j][k][l] =
            sqrt_g_det * epsilon4_symbol[i][j][k][l];
    }
}
template_GQ void GeometricQuantities_t::compute_Z_L_ST()
{
    if (m_Z_L_ST != nullptr)
        delete m_Z_L_ST;

    const auto &vars = get_vars();
    const auto &Z = get_Z();

    m_Z_L_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>;
    (*m_Z_L_ST)[0] = -vars.lapse * vars.Theta +
                     TensorAlgebra::compute_dot_product(vars.shift, Z);
    FOR(i) { (*m_Z_L_ST)[i + 1] = Z[i]; }
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

    m_grad_normal_LL = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    FOR_ST(m, n)
    {
        (*m_grad_normal_LL)[m][n] =
            -normal_L_ST[n] * d1_lapse_ST[m] / vars.lapse - Kij_ST[n][m];
    }
}
template_GQ void GeometricQuantities_t::compute_covd_Z_L_ST()
{
    if (m_covd_Z_L_ST != nullptr)
        delete m_covd_Z_L_ST;

    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();
    const auto &LIE = get_lie_derivatives();
    const auto &lie_Z = get_lie_Z();
    const auto &Z_U = get_Z_U();
    const auto &Kij = get_extrinsic_curvature();
    const auto &covd_Z = get_covd_Z();
    const auto &shift_ST = get_shift_ST();
    const auto &normal_L_ST = get_normal_L_ST();

    // Call the projections of Z by Q (Theta), Qi, Qj, and Qij

    data_t Z_U_dot_d1_lapse = 0.;
    Tensor<1, data_t> Kij_dot_Z = {0.};
    Tensor<1, data_t> Qi;
    Tensor<1, data_t> Qj;
    Tensor<2, data_t> Qij;
    FOR(i)
    {
        Z_U_dot_d1_lapse += Z_U[i] * d1.lapse[i];

        FOR(j)
        {
            Kij_dot_Z[i] += Kij[i][j] * Z_U[j];

            Qij[i][j] = -Kij[i][j] * vars.Theta + covd_Z[i][j];
        }

        Qi[i] = Kij_dot_Z[i] + vars.Theta * d1.lapse[i] / vars.lapse + lie_Z[i];
        Qj[i] = Kij_dot_Z[i] - d1.Theta[i];
    }

    const data_t Q = -LIE.Theta - Z_U_dot_d1_lapse / vars.lapse;

    const auto Qij_ST = TensorAlgebra::make_spatial_tensor_ST(Qij, shift_ST);
    const auto Qi_ST = TensorAlgebra::make_spatial_tensor_ST(Qi, shift_ST);
    const auto Qj_ST = TensorAlgebra::make_spatial_tensor_ST(Qj, shift_ST);

    m_covd_Z_L_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    FOR_ST(m, n)
    {
        (*m_covd_Z_L_ST)[m][n] = Q * normal_L_ST[m] * normal_L_ST[n] -
                                 normal_L_ST[m] * Qj_ST[n] -
                                 normal_L_ST[n] * Qi_ST[m] + Qij_ST[m][n];
    }
}
template_GQ void GeometricQuantities_t::compute_acceleration_ST()
{
    if (m_acceleration_ST != nullptr)
        delete m_acceleration_ST;

    const auto &vars = get_vars();
    const auto &acceleration = get_acceleration_spatial();

    m_acceleration_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>;
    (*m_acceleration_ST)[0] =
        TensorAlgebra::compute_dot_product(vars.shift, acceleration);
    FOR(i) { (*m_acceleration_ST)[i + 1] = acceleration[i]; }
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_em_tensor_ST()
{
    if (m_em_tensor_ST != nullptr)
        delete m_em_tensor_ST;

    const auto &em_tensor = get_em_tensor();
    const auto &normal_L_ST = get_normal_L_ST();
    const auto &shift_ST = get_shift_ST();

    const auto Smn_ST =
        TensorAlgebra::make_spatial_tensor_ST(em_tensor.Sij, shift_ST);
    const auto Sm_ST =
        TensorAlgebra::make_spatial_tensor_ST(em_tensor.Si, shift_ST);

    m_em_tensor_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    FOR_ST(m, n)
    {
        (*m_em_tensor_ST)[m][n] =
            em_tensor.rho * normal_L_ST[m] * normal_L_ST[n] +
            normal_L_ST[m] * Sm_ST[n] + normal_L_ST[n] * Sm_ST[m] +
            Smn_ST[m][n];
    }
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

    const auto &Eij = get_weyl_electric_part();
    const auto &Bij = get_weyl_magnetic_part();

    m_weyl_tensor_LLLL = new Tensor<4, data_t, CH_SPACEDIM + 1>(
        compute_weyl_tensor_LLLL(Eij, Bij));
}
template_GQ Tensor<4, data_t, CH_SPACEDIM + 1>
GeometricQuantities_t::compute_weyl_tensor_LLLL(const Tensor<2, data_t> &Eij,
                                                const Tensor<2, data_t> &Bij)
{

    const auto &g = get_metric_ST();
    const auto &g_UU = get_metric_UU_ST();
    const auto &n_L = get_normal_L_ST();
    const auto &levi_civita_spatial_ST = get_levi_civita_spatial_ST();
    const auto &shift_ST = get_shift_ST();

    Tensor<2, data_t, 4> l_LL; // l[a][b] = g[a][b] + 2n[a]n[b]
    Tensor<3, data_t, CH_SPACEDIM + 1> epsilon3_ULL = {0.};
    FOR_ST(a, b)
    {
        l_LL[a][b] = g[a][b] + 2. * n_L[a] * n_L[b];

        FOR_ST(c, d)
        {
            epsilon3_ULL[a][b][c] +=
                g_UU[a][d] * levi_civita_spatial_ST[d][b][c];
        }
    }

    const auto E_LL = TensorAlgebra::make_spatial_tensor_ST(Eij, shift_ST);
    const auto B_LL = TensorAlgebra::make_spatial_tensor_ST(Bij, shift_ST);

    Tensor<3, data_t, CH_SPACEDIM + 1> B_dot_epsilon = 0.;
    FOR_ST(a, b, c, d)
    {
        B_dot_epsilon[a][b][c] += B_LL[a][d] * epsilon3_ULL[d][b][c];
    }

    // Weyl
    Tensor<4, data_t, CH_SPACEDIM + 1> weyl_LLLL;
    FOR_ST(m, n, r, s)
    {
        weyl_LLLL[m][n][r][s] =
            (l_LL[m][r] * E_LL[s][n] - l_LL[m][s] * E_LL[r][n]) -
            (l_LL[n][r] * E_LL[s][m] - l_LL[n][s] * E_LL[r][m]) -
            (n_L[r] * B_dot_epsilon[s][m][n] -
             n_L[s] * B_dot_epsilon[r][m][n]) -
            (n_L[m] * B_dot_epsilon[n][r][s] - n_L[n] * B_dot_epsilon[m][r][s]);
    }

    return weyl_LLLL;
}
template_GQ void GeometricQuantities_t::compute_weyl_squared()
{
    if (m_weyl_squared != nullptr)
        delete m_weyl_squared;

    const auto &metric_UU = get_metric_UU_spatial();
    const auto &Eij = get_weyl_electric_part();
    const auto &Bij = get_weyl_magnetic_part();

    m_weyl_squared = new data_t(0.);
    FOR(i, j, k, l)
    {
        (*m_weyl_squared) += 8. * metric_UU[i][k] * metric_UU[j][l] *
                             (Eij[i][j] * Eij[k][l] - Bij[i][j] * Bij[k][l]);
    }
}
//////////////////////////////////////////////////////////////////////////
template_GQ void GeometricQuantities_t::compute_riemann_LLLL_ST()
{
    if (m_riemann_LLLL_ST != nullptr)
        delete m_riemann_LLLL_ST;

    const auto &weyl = get_weyl_tensor_LLLL();
    const auto &ricci = get_ricci_ST();
    const auto &ricci_scalar = get_ricci_scalar_ST();
    const auto &g = get_metric_ST();

    m_riemann_LLLL_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>;
    FOR_ST(a, b, m, n)
    {
        (*m_riemann_LLLL_ST)[a][b][m][n] =
            weyl[a][b][m][n] +
            (g[a][m] * ricci[n][b] - g[a][n] * ricci[m][b] -
             g[b][m] * ricci[n][a] + g[b][n] * ricci[m][a]) /
                (GR_SPACEDIM - 1.) -
            (g[a][m] * g[n][b] - g[a][n] * g[m][b]) * ricci_scalar /
                (GR_SPACEDIM * (GR_SPACEDIM - 1.));
    }
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

    const auto GC_4D =
        TensorAlgebra::make_spatial_tensor_ST(gauss_codazzi, shift_ST);
    const auto CM_4D =
        TensorAlgebra::make_spatial_tensor_ST(codazzi_mainardi, shift_ST);
    const auto EOMnn_4D = TensorAlgebra::make_spatial_tensor_ST(
        eom_double_normal_projection, shift_ST);

    m_riemann_LLLL_ST_v2 = new Tensor<4, data_t, CH_SPACEDIM + 1>;
    FOR_ST(m, n, r, s)
    {
        (*m_riemann_LLLL_ST_v2)[m][n][r][s] =
            GC_4D[m][n][r][s] +
            (-n_L[s] * CM_4D[m][n][r] + n_L[r] * CM_4D[m][n][s] -
             n_L[n] * CM_4D[r][s][m] + n_L[m] * CM_4D[r][s][n]) +
            (n_L[s] * n_L[n] * EOMnn_4D[m][r] -
             n_L[r] * n_L[n] * EOMnn_4D[m][s] -
             n_L[s] * n_L[m] * EOMnn_4D[n][r] +
             n_L[r] * n_L[m] * EOMnn_4D[n][s]);
    }
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

    if (m_em_tensor != nullptr)
    {
        // compute these tensors! To use below directly
        get_em_tensor_ST();
        get_em_tensor_trace_ST();
    }

    data_t kappa1 = m_ccz4_params->kappa1;
    if (m_ccz4_params->covariantZ4)
        kappa1 /= vars.lapse;

    data_t Z_dot_n = 0.;
    if (m_formulation == CCZ4RHS<>::USE_CCZ4)
        Z_dot_n = TensorAlgebra::compute_dot_product(Z_L_ST, n_U);

    m_ricci_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    FOR_ST(m, n)
    {
        (*m_ricci_ST)[m][n] =
            m_16_pi_G_Newton / 2. *
                (m_em_tensor == nullptr
                     ? 0.
                     : (*m_em_tensor_ST)[m][n] - g[m][n] *
                                                     (*m_em_tensor_trace_ST) /
                                                     (GR_SPACEDIM - 1.)) +
            m_cosmological_constant * g[m][n] * 2. / (GR_SPACEDIM - 1.);

        if (m_formulation == CCZ4RHS<>::USE_CCZ4)
        {
            (*m_ricci_ST)[m][n] +=
                -covd_Z_L_ST[m][n] - covd_Z_L_ST[n][m] +
                kappa1 * (n_L[m] * Z_L_ST[n] + n_L[n] * Z_L_ST[m] -
                          2. / (GR_SPACEDIM - 1.) *
                              (1. + m_ccz4_params->kappa2) * g[m][n] * Z_dot_n);
        }
    }
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

    m_ricci_squared = new data_t(0.);
    FOR_ST(a, b, c, d)
    {
        (*m_ricci_squared) +=
            g_UU[a][c] * g_UU[b][d] * ricci[a][b] * ricci[c][d];
    }
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

    const auto &riemann_LULU = get_riemann_LULU_ST();

    m_riemann_squared = new data_t(0.);
    FOR_ST(a, b, c, d)
    {
        (*m_riemann_squared) +=
            riemann_LULU[a][b][c][d] * riemann_LULU[b][a][d][c];
    }
}
template_GQ void GeometricQuantities_t::compute_riemann_LLLU_ST()
{
    if (m_riemann_LLLU_ST != nullptr)
        delete m_riemann_LLLU_ST;

    const auto &g_UU = get_metric_UU_ST();
    const auto &riemann_LLLL = get_riemann_LLLL_ST();

    m_riemann_LLLU_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>({0.});
    FOR_ST(a, b, c, d, e)
    {
        (*m_riemann_LLLU_ST)[a][b][c][d] +=
            riemann_LLLL[a][b][c][e] * g_UU[e][d];
    }
}
template_GQ void GeometricQuantities_t::compute_riemann_LULU_ST()
{
    if (m_riemann_LULU_ST != nullptr)
        delete m_riemann_LULU_ST;

    const auto &g_UU = get_metric_UU_ST();
    const auto &riemann_LLLU = get_riemann_LLLU_ST();

    m_riemann_LULU_ST = new Tensor<4, data_t, CH_SPACEDIM + 1>({0.});
    FOR_ST(a, b, c, d, e)
    {
        (*m_riemann_LULU_ST)[a][b][c][d] +=
            riemann_LLLU[a][e][c][d] * g_UU[e][b];
    }
}
//////////////////////////////////////////////////////////////////////////

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

    /*
        const Tensor<2, data_t> covd_shift =
       TensorAlgebra::covariant_derivative( d1.shift, vars.shift, chris_spatial,
       {true}); const Tensor<1, data_t> shift_dot_Kij =
            TensorAlgebra::compute_dot_product(vars.shift, Kij);
        const data_t shift_dot_shift_dot_Kij =
            TensorAlgebra::compute_dot_product(vars.shift, shift_dot_Kij);
    const Tensor<1, data_t> shift_dot_covd_shift =
        TensorAlgebra::compute_dot_product(vars.shift, covd_shift);
    const Tensor<2, data_t> Kij_UL =
        TensorAlgebra::compute_dot_product(metric_UU_spatial, Kij, 0, 0);
    */

    Tensor<2, data_t> covd_shift;
    Tensor<1, data_t> shift_dot_Kij = 0.;
    Tensor<2, data_t> Kij_UL = 0.;

    FOR(i, j)
    {
        covd_shift[i][j] = d1.shift[i][j];
        FOR(k)
        {
            covd_shift[i][j] += chris_spatial[i][j][k] * vars.shift[k];

            Kij_UL[i][j] += metric_UU_spatial[i][k] * Kij[k][j];
        }

        shift_dot_Kij[i] += vars.shift[j] * Kij[i][j];
    }

    data_t shift_dot_shift_dot_Kij = 0.;
    Tensor<1, data_t> shift_dot_covd_shift = 0.;

    FOR(i)
    {
        shift_dot_shift_dot_Kij += vars.shift[i] * shift_dot_Kij[i];

        FOR(j) { shift_dot_covd_shift[i] += vars.shift[j] * covd_shift[i][j]; }
    }

    m_chris_ST = new Tensor<3, data_t, CH_SPACEDIM + 1>;

    (*m_chris_ST)[0][0][0] =
        (rhs.lapse + advec.lapse - shift_dot_shift_dot_Kij) / vars.lapse;

    Tensor<1, data_t> aux_term;
    FOR(i)
    {
        (*m_chris_ST)[0][0][i + 1] =
            (d1.lapse[i] - shift_dot_Kij[i]) / vars.lapse;
        (*m_chris_ST)[0][i + 1][0] = (*m_chris_ST)[0][0][i + 1];

        FOR(j) { (*m_chris_ST)[0][i + 1][j + 1] = -Kij[i][j] / vars.lapse; }

        aux_term[i] = d1.lapse[i] - 2. * shift_dot_Kij[i];
    }
    const Tensor<1, data_t> aux_term_U =
        TensorAlgebra::compute_dot_product(aux_term, metric_UU_spatial);

    FOR(i)
    {
        (*m_chris_ST)[i + 1][0][0] = vars.lapse * aux_term_U[i] -
                                     vars.shift[i] * (*m_chris_ST)[0][0][0] +
                                     rhs.shift[i] + shift_dot_covd_shift[i];

        FOR(j)
        {
            (*m_chris_ST)[i + 1][0][j + 1] =
                -vars.shift[i] * (*m_chris_ST)[0][0][j + 1] -
                vars.lapse * Kij_UL[i][j] + covd_shift[i][j];
            (*m_chris_ST)[i + 1][j + 1][0] = (*m_chris_ST)[i + 1][0][j + 1];

            FOR(k)
            {
                (*m_chris_ST)[i + 1][j + 1][k + 1] =
                    chris_spatial[i][j][k] +
                    vars.shift[i] * Kij[j][k] / vars.lapse;
            }
        }
    }
}
template_GQ void GeometricQuantities_t::compute_Gamma_ST()
{
    if (m_Gamma_ST != nullptr)
        delete m_Gamma_ST;
    m_Gamma_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(
        TensorAlgebra::compute_trace(get_chris_ST(), get_metric_UU_ST()));
}
template_GQ void GeometricQuantities_t::compute_Gamma_L_ST()
{
    if (m_Gamma_L_ST != nullptr)
        delete m_Gamma_L_ST;
    m_Gamma_L_ST = new Tensor<1, data_t, CH_SPACEDIM + 1>(
        TensorAlgebra::compute_dot_product(get_Gamma_ST(), get_metric_ST()));
}
/*
template_GQ void GeometricQuantities_t::compute_d1_Z_L_ST()
{
    if (m_d1_Z_L_ST != nullptr)
        delete m_d1_Z_L_ST;

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

    // const auto &d1_shift_dot_Z =
    //     TensorAlgebra::compute_dot_product(Z, d1.shift, 0, 0);
    // const auto &shift_dot_d1_Z =
    //     TensorAlgebra::compute_dot_product(vars.shift, d1_Z, 0, 0);

    data_t rhs_shift_dot_Z = 0.;
    data_t shift_dot_dt_Z = 0.;
    Tensor<1, data_t> &d1_shift_dot_Z = 0.;
    Tensor<1, data_t> &shift_dot_d1_Z = 0.;

    FOR(i)
    {
        rhs_shift_dot_Z += rhs.shift[i] * Z[i];
        shift_dot_dt_Z += vars.shift[i] * dt_Z[i];

        FOR(j)
        {
            d1_shift_dot_Z[i] += Z[j] * d1.shift[j][i];

            shift_dot_d1_Z[j] += vars.shift[j] * d1_Z[j][i];
        }
    }

    m_d1_Z_L_ST = new Tensor<2, data_t, CH_SPACEDIM + 1>;
    (*m_d1_Z_L_ST)[0][0] = -vars.lapse * rhs.Theta - rhs.lapse * vars.Theta +
                           rhs_shift_dot_Z + shift_dot_dt_Z;

    FOR(i)
    {
        (*m_d1_Z_L_ST)[i + 1][0] = dt_Z[i];
        (*m_d1_Z_L_ST)[0][i + 1] = -vars.lapse * d1.Theta[i] -
                                   d1.lapse[i] * vars.Theta +
                                   d1_shift_dot_Z[i] + shift_dot_d1_Z[i];

        FOR(j) { (*m_d1_Z_L_ST)[i + 1][j + 1] = d1_Z[i][j]; }
    }
}
*/
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template_GQ Tensor<2, data_t>
GeometricQuantities_t::compute_LieD_weyl_electric_part(
    const Tensor<2, Tensor<1, data_t>> &d1Eij,
    const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
    const Tensor<2, data_t> &Bij)
{
    const auto &vars = get_vars();

    const auto &metric_spatial = get_metric_spatial();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const auto &extrinsic_curvature = get_extrinsic_curvature();
    const auto &acceleration = get_acceleration_spatial();
    const auto &chris_spatial = get_chris_spatial();
    const auto &levi_civita_LUU = get_levi_civita_spatial_LUU();

    Tensor<3, data_t> covd_B =
        TensorAlgebra::covariant_derivative(d1Bij, Bij, chris_spatial);

    Tensor<2, data_t> Eij_dot_Kij = TensorAlgebra::compute_dot_product(
        Eij, extrinsic_curvature, metric_UU_spatial);

    data_t Eij_dot_Kij_trace =
        TensorAlgebra::compute_trace(Eij_dot_Kij, metric_UU_spatial);

    Tensor<2, data_t> LieD_E;
    FOR(i, j)
    {
        LieD_E[i][j] = 2. * vars.K * Eij[i][j] +
                       Eij_dot_Kij_trace * metric_spatial[i][j] -
                       5. * Eij_dot_Kij[i][j];
        FOR(k, l)
        {
            LieD_E[i][j] +=
                levi_civita_LUU[i][k][l] * covd_B[k][j][l] -
                2. * levi_civita_LUU[i][k][l] * acceleration[k] * Bij[j][l];
        }
    }

    TensorAlgebra::make_symmetric(LieD_E);

    return LieD_E;
}
template_GQ Tensor<2, data_t>
GeometricQuantities_t::compute_LieD_weyl_magnetic_part(
    const Tensor<2, Tensor<1, data_t>> &d1Eij,
    const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
    const Tensor<2, data_t> &Bij)
{
    const auto &vars = get_vars();

    const auto &metric_spatial = get_metric_spatial();
    const auto &metric_UU_spatial = get_metric_UU_spatial();
    const auto &extrinsic_curvature = get_extrinsic_curvature();
    const auto &acceleration = get_acceleration_spatial();
    const auto &chris_spatial = get_chris_spatial();
    const auto &levi_civita_LUU = get_levi_civita_spatial_LUU();

    Tensor<3, data_t> covd_E =
        TensorAlgebra::covariant_derivative(d1Eij, Eij, chris_spatial);

    Tensor<2, data_t> Bij_dot_Kij = TensorAlgebra::compute_dot_product(
        Bij, extrinsic_curvature, metric_UU_spatial);

    data_t Bij_dot_Kij_trace =
        TensorAlgebra::compute_trace(Bij_dot_Kij, metric_UU_spatial);

    Tensor<2, data_t> LieD_B;
    FOR(i, j)
    {
        LieD_B[i][j] = 2. * vars.K * Bij[i][j] +
                       Bij_dot_Kij_trace * metric_spatial[i][j] -
                       5. * Bij_dot_Kij[i][j];
        FOR(k, l)
        {
            LieD_B[i][j] +=
                -levi_civita_LUU[i][k][l] * covd_E[k][j][l] +
                2. * levi_civita_LUU[i][k][l] * acceleration[k] * Eij[j][l];
        }
    }

    TensorAlgebra::make_symmetric(LieD_B);

    return LieD_B;
}
template_GQ Tensor<2, data_t>
GeometricQuantities_t::compute_dt_weyl_electric_part(
    const Tensor<2, Tensor<1, data_t>> &d1Eij,
    const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
    const Tensor<2, data_t> &Bij, const Tensor<2, data_t> &advec_Eij,
    const Tensor<2, data_t> &advec_Bij)
{
    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();

    const Tensor<2, data_t> LieD_E =
        compute_LieD_weyl_electric_part(d1Eij, d1Bij, Eij, Bij);

    Tensor<2, data_t> dt_E;
    FOR(i, j)
    {
        dt_E[i][j] = advec_Eij[i][j] + LieD_E[i][j] * vars.lapse;
        FOR1(k)
        {
            dt_E[i][j] +=
                Eij[i][k] * d1.shift[k][j] + Eij[j][k] * d1.shift[k][i];
        }
    }
    return dt_E;
}
template_GQ Tensor<2, data_t>
GeometricQuantities_t::compute_dt_weyl_magnetic_part(
    const Tensor<2, Tensor<1, data_t>> &d1Eij,
    const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
    const Tensor<2, data_t> &Bij, const Tensor<2, data_t> &advec_Eij,
    const Tensor<2, data_t> &advec_Bij)
{
    const auto &vars = get_vars();
    const auto &d1 = get_d1_vars();

    const Tensor<2, data_t> LieD_B =
        compute_LieD_weyl_magnetic_part(d1Eij, d1Bij, Eij, Bij);

    Tensor<2, data_t> dt_B;
    FOR(i, j)
    {
        dt_B[i][j] = advec_Bij[i][j] + LieD_B[i][j] * vars.lapse;
        FOR1(k)
        {
            dt_B[i][j] +=
                Bij[i][k] * d1.shift[k][j] + Bij[j][k] * d1.shift[k][i];
        }
    }
    return dt_B;
}
#endif /* GEOMETRICQUANTITIES_IMPL_HPP_ */
