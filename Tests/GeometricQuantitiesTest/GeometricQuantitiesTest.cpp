/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to LICENSE, in Chombo's root directory.
 */
#endif

#include "parstream.H" //Gives us pout()

// #include "GeometricQuantities.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"
#include "ScalarField.hpp"

#include "OldCCZ4.hpp"

template <class data_t>
using MatterVars = MatterCCZ4RHS<ScalarField<>>::Vars<data_t>;
template <class data_t> struct Vars : public MatterVars<data_t>
{
    Tensor<2, data_t> WeylE;
    Tensor<2, data_t> WeylB;

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        MatterVars<data_t>::enum_mapping(mapping_function);

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        // Symmetric 2-tensors
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_E11, c_E33>(), WeylE);
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_B11, c_B33>(), WeylB);
    }
};
template <class data_t>
using Diff2Vars = MatterCCZ4RHS<ScalarField<>>::Diff2Vars<data_t>;

template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank == 0), bool>::type
relative_error_aux(const typename Tensor<rank, data_t, size>::arr_t &t1,
                   const typename Tensor<rank, data_t, size>::arr_t &t2,
                   const std::string &name, double err_max = 1.e-12)
{
    double err = (t1 == 0. || t2 == 0. ? abs(t1 - t2) : abs(t1 / t2 - 1.));
    bool error = err > err_max;
    if (error)
    {
        std::cout << "Error in " << name << ". |";
        if (t1 == 0. || t2 == 0.)
            std::cout << t1 << " - " << t2;
        else
            std::cout << t1 << " / " << t2 << "- 1 ";
        std::cout << "| = " << err << std::endl;
    }
    return error;
}

template <class data_t, int size, int rank>
ALWAYS_INLINE typename std::enable_if<(rank > 0), bool>::type
relative_error_aux(const typename Tensor<rank, data_t, size>::arr_t &t1,
                   const typename Tensor<rank, data_t, size>::arr_t &t2,
                   const std::string &name, double err_max = 1.e-12)
{
    bool out = false;
    for (int i = 0; i < size; ++i)
        out |= relative_error_aux<data_t, size, rank - 1>(
            t1[i], t2[i], name + "[" + std::to_string(i) + "]", err_max);
    return out;
}

template <class data_t, int size, int rank>
ALWAYS_INLINE bool relative_error(const Tensor<rank, data_t, size> &t1,
                                  const Tensor<rank, data_t, size> &t2,
                                  const std::string &name,
                                  double err_max = 1.e-12)
{
    return relative_error_aux<data_t, size, rank>(
        (const typename Tensor<rank, data_t, size>::arr_t &)t1,
        (const typename Tensor<rank, data_t, size>::arr_t &)t2, name, err_max);
}

template <class data_t>
ALWAYS_INLINE bool relative_error(const data_t &t1, const data_t &t2,
                                  const std::string &name,
                                  double err_max = 1.e-12)
{
    return relative_error_aux<data_t, 3, 0>(
        (const typename Tensor<0, data_t>::arr_t &)t1,
        (const typename Tensor<0, data_t>::arr_t &)t2, name, err_max);
}

int runTest(int argc, char *argv[])
{
    int formulation = CCZ4RHS<>::USE_CCZ4;

    CCZ4_params_t<> ccz4_params;

    // Lapse evolution
    ccz4_params.lapse_advec_coeff = 1.;
    ccz4_params.lapse_coeff = 2.;
    ccz4_params.lapse_power = 1.;

    // Shift Evolution
    ccz4_params.shift_advec_coeff = 0.;
    ccz4_params.shift_Gamma_coeff = 0.75;
    ccz4_params.eta = 1.;

    // CCZ4 parameters
    ccz4_params.kappa1 = 0.1;
    ccz4_params.kappa2 = 0.;
    ccz4_params.kappa3 = 1.;
    ccz4_params.covariantZ4 = true;

    int failed = 0;

    //////////////////////////////////////////////////////////////////
    // SET VARIABLES

    Vars<double> vars;
    Vars<Tensor<1, double>> d1;
    Diff2Vars<Tensor<2, double>> d2;
    double cosmological_constant;

#include "randomCCZ4VarsValues.hpp" //Including the auto generated file with variables

    //////////////////////////////////////////////////////////////////
    // SET GEOMETRIC QUANTITIES OBJECT

    // make advection
    Vars<double> advec;
    advec.enum_mapping([&](const int &ivar, double &var) {
        d1.enum_mapping([&](const int &d1_ivar, Tensor<1, double> &d1_var) {
            if (ivar == d1_ivar)
            {
                var = 0.;
                FOR(dir) { var += vars.shift[dir] * d1_var[dir]; }
            }
        });
    });

    // make gauge
    MovingPunctureGauge gauge(ccz4_params);

    GeometricQuantities<double, Vars, Diff2Vars, MovingPunctureGauge> gq(
        vars, d1, d2, "GeometricQuantitiesTest::runTest");
    gq.set_formulation(formulation, ccz4_params);
    gq.set_advection_and_gauge(advec, gauge);
    gq.set_cosmological_constant(cosmological_constant);

    // make EM Tensor
    DefaultPotential potential;
    ScalarField<> sf(potential);
    auto em_tensor_def = sf.compute_emtensor(gq);

    gq.set_em_tensor(em_tensor_def, 1.);

    //////////////////////////////////////////////////////////////////
    // SET MATHEMATICA CALCULATIONS

    emtensor_t<double> em_tensor;

    // spatial conformal
    Tensor<2, double> h_UU;
    chris_t<double> chris;
    Tensor<1, double> Z_U_conformal;
    Tensor<2, double> d1_chris_contracted;
    Tensor<2, double> covd_chi_conformal;
    Tensor<2, double> A_LU;
    Tensor<2, double> A_UU;
    double tr_A2;
    double div_shift;
    Tensor<1, double> Gamma_L;

    // spatial
    Tensor<1, double> shift_L;
    Tensor<2, double> metric_spatial;
    Tensor<2, double> metric_UU_spatial;
    Tensor<2, double> Kij;
    Tensor<3, double> chris_spatial;
    Tensor<1, double> Z_U;
    Tensor<1, double> Z;
    Tensor<2, double> covd_Z;
    Tensor<3, double> d1_Kij;
    Tensor<3, double> covd_Kij;
    Tensor<3, double> levi_civita_spatial;
    Tensor<3, double> levi_civita_spatial_LUU;
    Tensor<2, double> covd_lapse;
    Tensor<1, double> Gamma_spatial;
    Tensor<1, double> Gamma_L_spatial;
    Tensor<1, double> acceleration_spatial;

    // xTensor 3D
    ricci_t<double> ricci;
    ricci_t<double> ricci_1DZ;
    ricci_t<double> ricci_2DZ;
    Tensor<4, double> riemann_conformal_LLLL;
    Tensor<4, double> riemann_spatial_LLLL;
    Tensor<4, double> gauss_codazzi;
    Tensor<3, double> codazzi_mainardi;
    double Ham;
    Tensor<1, double> Mom;
    Tensor<2, double> Eij;
    Tensor<2, double> Bij;
    Vars<double> LIE;
    Vars<double> RHS;
    Tensor<2, double> LieKij;
    Tensor<2, double> eom_double_normal_projection;

    // 4D
    Tensor<2, double, 4> metric_ST;
    Tensor<2, double, 4> projector_LU_ST;
    Tensor<2, double, 4> metric_UU_ST;
    Tensor<1, double, 4> normal_U_ST;
    Tensor<1, double, 4> normal_L_ST;
    Tensor<1, double, 4> shift_ST;
    Tensor<3, double, 4> levi_civita_spatial_ST;
    Tensor<4, double, 4> levi_civita_ST;
    Tensor<1, double, 4> Z_L_ST;
    Tensor<2, double, 4> grad_normal_LL;
    Tensor<1, double, 4> acceleration_ST;
    Tensor<2, double, 4> Tmn;
    double TrTmn;

    // auxiliary to Mathematica for xTensor 4D calculations
    double dtlapse = gq.get_rhs_equations().lapse;
    double dtshift1 = gq.get_rhs_equations().shift[0];
    double dtshift2 = gq.get_rhs_equations().shift[1];
    double dtshift3 = gq.get_rhs_equations().shift[2];

    // xTensor 4D
    Tensor<3, double, 4> chris4D;
    // Tensor<1, double> dtGammaC; // code commented, but working!
    // Tensor<1, double> dtGamma;  // code commented, but working!
    Tensor<1, double> LieZi;
    // Tensor<2, double, 4> d1_Z_L_ST; // code commented, but working!
    Tensor<2, double, 4> covd_Z_L_ST;
    Tensor<4, double, 4> riemann_LLLL_ST;
    Tensor<2, double, 4> ricci_ST;
    double ricci_scalar_ST;
    double ricci_squared;
    double kretschmann;
    Tensor<4, double, 4> weyl_tensor_LLLL;
    double weyl_squared;
    Tensor<1, double, 4> Gamma_ST;
    Tensor<1, double, 4> Gamma_L_ST;
    Tensor<4, double, 4> riemann_LLLU_ST;
    Tensor<4, double, 4> riemann_LULU_ST;

    // extra
    Tensor<2, double> LieD_weyl_electric_part;
    Tensor<2, double> LieD_weyl_magnetic_part;
    Tensor<2, double> dt_weyl_electric_part;
    Tensor<2, double> dt_weyl_magnetic_part;

#include "randomGeometricQuantities.hpp" //Including the auto generated file with calculations

    //////////////////////////////////////////////////////////////////
    // CHECK
    failed |=
        relative_error(gq.get_em_tensor().rho, em_tensor.rho, "emtensor.rho");
    failed |=
        relative_error(gq.get_em_tensor().Si, em_tensor.Si, "emtensor.Si");
    failed |=
        relative_error(gq.get_em_tensor().Sij, em_tensor.Sij, "emtensor.Sij");
    failed |= relative_error(gq.get_em_tensor().S, em_tensor.S, "emtensor.S");

    // spatial conformal
    failed |= relative_error(gq.get_h_UU(), h_UU, "h_UU");
    failed |= relative_error(gq.get_chris().LLL, chris.LLL, "chris.LLL");
    failed |= relative_error(gq.get_chris().ULL, chris.ULL, "chris.ULL");
    failed |= relative_error(gq.get_chris().contracted, chris.contracted,
                             "chris.contracted");
    failed |=
        relative_error(gq.get_Z_U_conformal(), Z_U_conformal, "Z_U_conformal");
    failed |= relative_error(gq.get_d1_chris_contracted(), d1_chris_contracted,
                             "d1_chris_contracted");
    failed |= relative_error(gq.get_covd_chi_conformal(), covd_chi_conformal,
                             "covd_chi_conformal");
    failed |= relative_error(gq.get_A_LU(), A_LU, "A_LU");
    failed |= relative_error(gq.get_A_UU(), A_UU, "A_UU");
    failed |= relative_error(gq.get_tr_A2(), tr_A2, "tr_A2");
    failed |= relative_error(gq.get_div_shift(), div_shift, "div_shift");
    failed |= relative_error(gq.get_Gamma_L(), Gamma_L, "Gamma_L");

    // spatial
    failed |= relative_error(gq.get_shift_L(), shift_L, "shift_L");
    failed |= relative_error(gq.get_metric_spatial(), metric_spatial,
                             "metric_spatial");
    failed |= relative_error(gq.get_metric_UU_spatial(), metric_UU_spatial,
                             "metric_UU_spatial");
    failed |= relative_error(gq.get_extrinsic_curvature(), Kij, "Kij");
    failed |=
        relative_error(gq.get_chris_spatial(), chris_spatial, "chris_spatial");
    failed |= relative_error(gq.get_Z_U(), Z_U, "Z_U");
    failed |= relative_error(gq.get_Z(), Z, "Z");
    failed |= relative_error(gq.get_covd_Z(), covd_Z, "covd_Z");
    failed |= relative_error(gq.get_d1_extrinsic_curvature(), d1_Kij, "d1_Kij");
    failed |=
        relative_error(gq.get_covd_extrinsic_curvature(), covd_Kij, "covd_Kij");
    failed |= relative_error(gq.get_levi_civita_spatial(), levi_civita_spatial,
                             "levi_civita_spatial");
    failed |=
        relative_error(gq.get_levi_civita_spatial_LUU(),
                       levi_civita_spatial_LUU, "levi_civita_spatial_LUU");
    failed |= relative_error(gq.get_covd_lapse(), covd_lapse, "covd_lapse");
    failed |=
        relative_error(gq.get_Gamma_spatial(), Gamma_spatial, "Gamma_spatial");
    failed |= relative_error(gq.get_Gamma_L_spatial(), Gamma_L_spatial,
                             "Gamma_L_spatial");
    failed |= relative_error(gq.get_acceleration_spatial(),
                             acceleration_spatial, "get_acceleration_spatial");

    // xTensor 3D
    failed |= relative_error(gq.get_ricci().LL, ricci.LL, "ricci.LL");
    failed |=
        relative_error(gq.get_ricci().scalar, ricci.scalar, "ricci.scalar");
    failed |=
        relative_error(gq.get_ricci_1DZ().LL, ricci_1DZ.LL, "ricci_1DZ.LL");
    failed |= relative_error(gq.get_ricci_1DZ().scalar, ricci_1DZ.scalar,
                             "ricci_1DZ.scalar");
    failed |=
        relative_error(gq.get_ricci_2DZ().LL, ricci_2DZ.LL, "ricci_2DZ.LL");
    failed |= relative_error(gq.get_ricci_2DZ().scalar, ricci_2DZ.scalar,
                             "ricci_2DZ.scalar");
    failed |= relative_error(gq.get_riemann_conformal_LLLL(),
                             riemann_conformal_LLLL, "riemann_conformal_LLLL");
    failed |= relative_error(gq.get_riemann_spatial_LLLL(),
                             riemann_spatial_LLLL, "riemann_spatial_LLLL");
    failed |=
        relative_error(gq.get_gauss_codazzi(), gauss_codazzi, "gauss_codazzi");
    failed |= relative_error(gq.get_codazzi_mainardi(), codazzi_mainardi,
                             "codazzi_mainardi");
    failed |= relative_error(gq.get_hamiltonian_constraint(), Ham, "Ham");
    failed |= relative_error(gq.get_momentum_constraints(), Mom, "Mom");
    failed |= relative_error(gq.get_weyl_magnetic_part(), Bij, "Bij");
    failed |= relative_error(gq.get_weyl_electric_part(), Eij, "Eij");
    failed |= relative_error(gq.get_lie_derivatives().chi, LIE.chi, "LIE.chi");
    failed |= relative_error(gq.get_lie_derivatives().h, LIE.h, "LIE.h");
    failed |= relative_error(gq.get_lie_derivatives().K, LIE.K, "LIE.K");
    failed |= relative_error(gq.get_lie_derivatives().A, LIE.A, "LIE.A");
    failed |=
        relative_error(gq.get_lie_derivatives().Theta, LIE.Theta, "LIE.Theta");
    failed |=
        relative_error(gq.get_lie_derivatives().Gamma, LIE.Gamma, "LIE.Gamma");
    failed |= relative_error(gq.get_rhs_equations().chi, RHS.chi, "RHS.chi");
    failed |= relative_error(gq.get_rhs_equations().h, RHS.h, "RHS.h");
    failed |= relative_error(gq.get_rhs_equations().K, RHS.K, "RHS.K");
    failed |= relative_error(gq.get_rhs_equations().A, RHS.A, "RHS.A");
    failed |=
        relative_error(gq.get_rhs_equations().Theta, RHS.Theta, "RHS.Theta");
    failed |=
        relative_error(gq.get_rhs_equations().Gamma, RHS.Gamma, "RHS.Gamma");
    failed |= relative_error(gq.get_lie_extrinsic_curvature(), LieKij, "LieKij",
                             1.e-11);
    failed |= relative_error(gq.get_eom_double_normal_projection(),
                             eom_double_normal_projection,
                             "eom_double_normal_projection");

    // 4D
    failed |= relative_error(gq.get_metric_ST(), metric_ST, "metric_ST");
    failed |= relative_error(gq.get_projector_LU_ST(), projector_LU_ST,
                             "projector_LU_ST");
    failed |=
        relative_error(gq.get_metric_UU_ST(), metric_UU_ST, "metric_UU_ST");
    failed |= relative_error(gq.get_normal_U_ST(), normal_U_ST, "normal_U_ST");
    failed |= relative_error(gq.get_normal_L_ST(), normal_L_ST, "normal_L_ST");
    failed |= relative_error(gq.get_shift_ST(), shift_ST, "shift_ST");
    failed |= relative_error(gq.get_levi_civita_ST(), levi_civita_ST,
                             "levi_civita_ST");
    failed |= relative_error(gq.get_levi_civita_spatial_ST(),
                             levi_civita_spatial_ST, "levi_civita_spatial_ST");
    failed |= relative_error(gq.get_Z_L_ST(), Z_L_ST, "Z_L_ST");
    failed |= relative_error(gq.get_grad_normal_LL(), grad_normal_LL,
                             "grad_normal_LL");
    failed |= relative_error(gq.get_acceleration_ST(), acceleration_ST,
                             "get_acceleration_ST");
    failed |= relative_error(gq.get_em_tensor_ST(), Tmn, "Tmn");
    failed |= relative_error(gq.get_em_tensor_trace_ST(), TrTmn, "TrTmn");

    // xTensor 4D
    failed |= relative_error(gq.get_chris_ST(), chris4D, "chris4D");
    // failed |=
    // relative_error(gq.get_dt_chris_contracted(), dtGammaC, "dtGammaC");
    // failed |= relative_error(gq.get_dt_chris_spatial_contracted(),
    // dtGamma, "dtGamma");
    failed |= relative_error(gq.get_lie_Z(), LieZi, "LieZi");
    // failed |= relative_error(gq.get_d1_Z_L_ST(), d1_Z_L_ST, "d1_Z_L_ST");
    failed |= relative_error(gq.get_covd_Z_L_ST(), covd_Z_L_ST, "covd_Z_L_ST");
    failed |= relative_error(gq.get_riemann_LLLL_ST(), riemann_LLLL_ST,
                             "riemann_LLLL_ST", 2.e-11);
    failed |= relative_error(gq.get_riemann_LLLL_ST_v2(), riemann_LLLL_ST,
                             "riemann_LLLL_ST_v2", 2.e-11);
    failed |= relative_error(gq.get_ricci_ST(), ricci_ST, "ricci_ST");
    failed |= relative_error(gq.get_ricci_scalar_ST(), ricci_scalar_ST,
                             "ricci_scalar_ST");
    failed |=
        relative_error(gq.get_ricci_squared(), ricci_squared, "ricci_squared");
    failed |= relative_error(gq.get_kretschmann(), kretschmann, "kretschmann");
    failed |= relative_error(gq.get_kretschmann(), gq.get_riemann_squared(),
                             "riemann_squared");
    failed |= relative_error(gq.get_weyl_tensor_LLLL(), weyl_tensor_LLLL,
                             "weyl_tensor_LLLL", 2.e-11);
    failed |=
        relative_error(gq.get_weyl_squared(), weyl_squared, "weyl_squared");
    failed |= relative_error(gq.get_Gamma_ST(), Gamma_ST, "Gamma_ST");
    failed |= relative_error(gq.get_Gamma_L_ST(), Gamma_L_ST, "Gamma_L_ST");
    failed |= relative_error(gq.get_riemann_LLLU_ST(), riemann_LLLU_ST,
                             "riemann_LLLU_ST", 2.e-11);
    failed |= relative_error(gq.get_riemann_LULU_ST(), riemann_LULU_ST,
                             "riemann_LULU_ST", 2.e-11);

    // extra
    failed |=
        relative_error(gq.compute_LieD_weyl_electric_part(
                           d1.WeylE, d1.WeylB, vars.WeylE, vars.WeylB),
                       LieD_weyl_electric_part, "LieD_weyl_electric_part");
    failed |=
        relative_error(gq.compute_LieD_weyl_magnetic_part(
                           d1.WeylE, d1.WeylB, vars.WeylE, vars.WeylB),
                       LieD_weyl_magnetic_part, "LieD_weyl_magnetic_part");
    failed |= relative_error(
        gq.compute_dt_weyl_electric_part(d1.WeylE, d1.WeylB, vars.WeylE,
                                         vars.WeylB, advec.WeylE, advec.WeylB),
        dt_weyl_electric_part, "dt_weyl_electric_part");
    failed |= relative_error(
        gq.compute_dt_weyl_magnetic_part(d1.WeylE, d1.WeylB, vars.WeylE,
                                         vars.WeylB, advec.WeylE, advec.WeylB),
        dt_weyl_magnetic_part, "dt_weyl_magnetic_part");

    // COMPARE WITH OLD CCZ4 FOR BSSN AS AN EXTRA TEST
    // Also good to test the rest of variables in GeometricQuantities
    gq.set_em_tensor(em_tensor_def, 0.); // set matter to 0 here
    formulation = CCZ4RHS<>::USE_BSSN;
    ccz4_params.kappa1 = 0.;
    ccz4_params.kappa2 = 0.;
    ccz4_params.kappa3 = 0.;
    gq.set_formulation(formulation, ccz4_params);
    vars.Theta = 0.;
    advec.Theta = 0.;
    FOR(i) { d1.Theta[i] = 0.; }

    OldCCZ4 ccz4(ccz4_params, 0., 0., formulation, cosmological_constant);
    Vars<double> old_rhs;
    ccz4.rhs_equation(old_rhs, vars, d1, d2, advec);

    failed |=
        relative_error(gq.get_rhs_equations().chi, old_rhs.chi, "old_rhs.chi");
    failed |= relative_error(gq.get_rhs_equations().h, old_rhs.h, "old_rhs.h");
    failed |= relative_error(gq.get_rhs_equations().K, old_rhs.K, "old_rhs.K");
    failed |= relative_error(gq.get_rhs_equations().A, old_rhs.A, "old_rhs.A");
    failed |= relative_error(gq.get_rhs_equations().Gamma, old_rhs.Gamma,
                             "old_rhs.Gamma");
    failed |= relative_error(gq.get_rhs_equations().Theta, old_rhs.Theta,
                             "old_rhs.Theta");
    failed |= relative_error(gq.get_rhs_equations().lapse, old_rhs.lapse,
                             "old_rhs.lapse");
    failed |= relative_error(gq.get_rhs_equations().shift, old_rhs.shift,
                             "old_rhs.shift");
    failed |= relative_error(gq.get_rhs_equations().B, old_rhs.B, "old_rhs.B");

    return failed;
}

int main(int argc, char *argv[])
{
    // mainSetup(argc, argv);

    int failed = runTest(argc, argv);

    if (failed == 0)
        pout() << "GeometricQuantities test passed." << std::endl;
    else
        pout() << "GeometricQuantities test failed with return code " << failed
               << std::endl;

    // mainFinalize();
    return failed;
}
