/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4_HPP_)
#error "This file should only be included through CCZ4.hpp"
#endif

#ifndef CCZ4_IMPL_HPP_
#define CCZ4_IMPL_HPP_

#define COVARIANTZ4
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

#include "GeometricQuantities.hpp"

#include "Coordinates.hpp"
#include "simd.hpp"
#include <iomanip> //std::setprecision
template <class data_t> bool check(data_t num1, data_t num2)
{
    return (abs(num1) < 1.e14 && abs(num2) < 1.e14 &&
            abs(num1 - num2) < 1e-16) ||
           abs(num1 / num2 - 1.) < 0.000001;
}

inline CCZ4::CCZ4(params_t params, double dx, double sigma, int formulation,
                  double cosmological_constant)
    : m_params(params), m_sigma(sigma), m_formulation(formulation),
      m_cosmological_constant(cosmological_constant), m_deriv(dx)
{
    // A user who wants to use BSSN should also have damping paramters = 0
    if (m_formulation == USE_BSSN)
    {
        if ((m_params.kappa1 != 0.) || (params.kappa2 != 0.) ||
            (params.kappa3 != 0.))
        {
            MayDay::Error("BSSN formulation is selected - CCZ4 kappa values "
                          "should be set to zero in params");
        }
    }
    if (m_formulation > USE_BSSN)
        MayDay::Error("The requested formulation is not supported");
}

#include <time.h>
#include <unistd.h>

template <class data_t> void test()
{
    Tensor<2, data_t, 2> t1;
    int idx = 0.;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            t1[i][j] = idx++;

    Tensor<4, data_t, 2> t2;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                for (int l = 0; l < 2; ++l)
                    t2[i][j][k][l] = idx++;

    Tensor<1, data_t, 2> t3;
    for (int i = 0; i < 2; ++i)
        t3[i] = idx++;

    Tensor<2, data_t, 2> t4 = {2., 1., 1., 2.};
    Tensor<3, data_t, 2> t5 = {2., 1., 1., 2., 2., 1., 1., 2.};
    Tensor<1, Tensor<1, data_t, 2>, 2> t6;
    t6[0][0] = 2.;
    t6[0][0] = 1.;
    t6[0][0] = 1.;
    t6[0][0] = 2.;

    std::cout << "HERE t1 " << std::endl;
    TensorAlgebra::print(t1);
    std::cout << "HERE t2 " << std::endl;
    TensorAlgebra::print(t2);
    std::cout << "HERE t3" << std::endl;
    TensorAlgebra::print(t3);
    std::cout << "HERE t4 " << std::endl;
    TensorAlgebra::print(t4);
    std::cout << "HERE MOVE INDEX " << std::endl;
    TensorAlgebra::print(TensorAlgebra::move_index(t2, 1, 3));
    std::cout << "HERE DOT PRODUCT " << std::endl;
    TensorAlgebra::print(TensorAlgebra::compute_dot_product(t1, t2, 0, 0));
    std::cout << "HERE DOT PRODUCT " << std::endl;
    TensorAlgebra::print(TensorAlgebra::compute_dot_product(t3, t3));
    std::cout << "HERE TRANSPOSE " << std::endl;
    TensorAlgebra::print(TensorAlgebra::transpose(t2, 1, 2));
    std::cout << "HERE TRACE " << std::endl;
    // TensorAlgebra::print(TensorAlgebra::compute_trace(t2));
    TensorAlgebra::print(TensorAlgebra::compute_trace(t2, t4));
    std::cout << "HERE SPATIAL TO ST " << std::endl;
    TensorAlgebra::print(TensorAlgebra::make_spatial_tensor_ST(t4, t3));
    std::cout << std::endl;

    std::cout << "HERE COV. DER. " << std::endl;
    // const std::array<const Tensor<4, data_t, 2> *, 2> tensors = {&t2, &t2};
    // TensorAlgebra::print(TensorAlgebra::add(tensors));
    // std::cout << "HERE8 " << std::endl;
    // TensorAlgebra::print(
    // TensorAlgebra::compute_covariant_derivative(t6, t3, t5, {1}));

    /*    int N = 1000000;
        data_t x = 0.;
        clock_t t = clock();
        for (int i = 0; i < N; ++i)
            x += TensorAlgebra::compute_covariant_derivative(t2, t5, t5,
                                                             {1})[0][0][0][0];
        t = clock() - t;
        printf("It took me %ld clicks (%f seconds).\n", t,
               ((float)t) / CLOCKS_PER_SEC);
        TensorAlgebra::print(
            TensorAlgebra::compute_covariant_derivative(t2, t5, t5, {1}));
        t = clock();
        for (int i = 0; i < N; ++i)
            x += TensorAlgebra::compute_covariant_derivative2(t2, t5, t5,
                                                              {1})[0][0][0][0];
        t = clock() - t;
        printf("It took me %ld clicks (%f seconds).\n", t,
               ((float)t) / CLOCKS_PER_SEC);
        TensorAlgebra::print(
            TensorAlgebra::compute_covariant_derivative2(t2, t5, t5, {1}));

        std::cout << x << std::endl;*/

    MayDay::Error("BUM");
}

template <class data_t> void CCZ4::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    // GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2,
    // m_formulation);
    // const auto emtensor3ST = gq.get_em_tensor_ST();
    // const auto &ricci_2DZ = gq.get_ricci_2DZ();
    // const auto &H = gq.get_hamiltonian_constraint();
    // const auto &M = gq.get_momentum_constraints();
    // const auto &RHS = gq.get_rhs_equations();
    // const auto &T = gq.get_em_tensor_trace_ST();
    // const auto &TmnST = gq.get_em_tensor_effective_ST();
    // const auto &Tmn = gq.get_em_tensor_effective();

    // test<data_t>();
    // MayDay::Error("BUM");

    std::array<data_t, CH_SPACEDIM> center = {0., 0., 0.};
    Coordinates<data_t> coords(current_cell, m_deriv.m_dx, center);
    if (abs(coords.x - 5 - m_deriv.m_dx / 2.) < m_deriv.m_dx / 2. &&
        abs(coords.y - 5 - m_deriv.m_dx / 2.) < m_deriv.m_dx / 2. &&
        abs(coords.z - 5 - m_deriv.m_dx / 2.) < m_deriv.m_dx / 2.)
    {
        std::cout << setprecision(9);
        std::cout << "coords: (" << coords.x << "," << coords.y << ","
                  << coords.z << ")" << std::endl;
        std::cout << "LAPSE" << std::endl;
        std::cout << vars.lapse << std::endl;

        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
        gq.set_formulation(m_formulation, m_params);
        gq.set_advection(advec);

        const auto &chrisST = gq.get_chris_ST();

        const auto &rhs = gq.get_rhs_equations();
        std::cout << "RHS" << std::endl;
        std::cout << rhs.lapse << std::endl;
        TensorAlgebra::print(rhs.shift);
        std::cout << rhs.chi << std::endl;
        TensorAlgebra::print(rhs.h);
        std::cout << rhs.K << std::endl;
        TensorAlgebra::print(rhs.A);
        std::cout << rhs.Theta << std::endl;
        TensorAlgebra::print(rhs.Gamma);
        std::cout << "DONE" << std::endl;
        std::cout << std::endl;

        const auto &d1_Z_L_ST = gq.get_d1_Z_L_ST();
        std::cout << "d1_Z_L_ST" << std::endl;
        TensorAlgebra::print(d1_Z_L_ST);
        std::cout << std::endl;

        const auto &covd_Z_L_ST = gq.get_covd_Z_L_ST();
        const auto &Z_L_ST = gq.get_Z_L_ST();
        const auto covd_Z_L_ST_v2 =
            TensorAlgebra::covariant_derivative(d1_Z_L_ST, Z_L_ST, chrisST);

        std::cout << "covd_Z_L_ST" << std::endl;
        TensorAlgebra::print(covd_Z_L_ST);
        std::cout << "covd_Z_L_ST_v2" << std::endl;
        TensorAlgebra::print(covd_Z_L_ST_v2);

        const auto &n_U_ST = gq.get_normal_L_ST();
        const auto &projector = gq.get_projector_LU_ST();

        Tensor<1, data_t, CH_SPACEDIM + 1> Qi1 =
            TensorAlgebra::compute_dot_product(
                TensorAlgebra::compute_dot_product(n_U_ST, covd_Z_L_ST),
                projector);
        Tensor<1, data_t, CH_SPACEDIM + 1> Qi2 =
            TensorAlgebra::compute_dot_product(
                TensorAlgebra::compute_dot_product(n_U_ST, covd_Z_L_ST_v2),
                projector);

        Tensor<2, data_t, CH_SPACEDIM + 1> Qij1 =
            TensorAlgebra::compute_dot_product(
                TensorAlgebra::compute_dot_product(projector, covd_Z_L_ST, 1,
                                                   0),
                projector, 1, 1);
        Tensor<2, data_t, CH_SPACEDIM + 1> Qij2 =
            TensorAlgebra::compute_dot_product(
                TensorAlgebra::compute_dot_product(projector, covd_Z_L_ST_v2, 1,
                                                   0),
                projector, 1, 1);

        const auto &Z_U = gq.get_Z_U();

        std::cout << "Qis" << std::endl;
        TensorAlgebra::print(Qi1);
        TensorAlgebra::print(Qi2);

        const auto &lie_Z = gq.get_lie_Z();
        const auto &Kij = gq.get_extrinsic_curvature();
        const Tensor<1, data_t> Kij_dot_Z =
            TensorAlgebra::compute_dot_product(Kij, Z_U);

        Tensor<1, data_t> Qi;
        FOR(i)
        {
            Qi[i] =
                Kij_dot_Z[i] + vars.Theta * d1.lapse[i] / vars.lapse + lie_Z[i];
        }

        std::cout << "Qi Real" << std::endl;
        TensorAlgebra::print(Qi);

        std::cout << "Qijs" << std::endl;
        TensorAlgebra::print(Qij1);
        TensorAlgebra::print(Qij2);

        const auto &covd_Z = gq.get_covd_Z();
        Tensor<2, data_t> Qij;
        FOR(i, j) { Qij[i][j] = -Kij[i][j] * vars.Theta + covd_Z[i][j]; }

        std::cout << "Qij Real" << std::endl;
        TensorAlgebra::print(Qij);

        MayDay::Error("BUM");
    }

    Vars<data_t> rhs;
    rhs_equation(rhs, vars, d1, d2, advec);

    m_deriv.add_dissipation(rhs, current_cell, m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void CCZ4::rhs_equation(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                        const vars_t<Tensor<1, data_t>> &d1,
                        const diff2_vars_t<Tensor<2, data_t>> &d2,
                        const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    auto h_UU = compute_inverse_sym(vars.h);
    auto chris = compute_christoffel(d1.h, h_UU);

    Tensor<1, data_t> Z_over_chi;
    Tensor<1, data_t> Z;

    if (m_formulation == USE_BSSN)
    {
        FOR1(i) Z_over_chi[i] = 0.0;
    }
    else
    {
        FOR1(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris.contracted[i]);
    }
    FOR1(i) Z[i] = vars.chi * Z_over_chi[i];

    auto ricci =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    // GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
    // gq.set_formulation(m_formulation, m_params);
    // const auto &ricci2 = gq.get_ricci_2DZ();

    // bool right = true;
    // FOR(i, j) { right &= check(ricci2.LL[i][j], ricci.LL[i][j]); }
    // right &= check(ricci2.scalar, ricci.scalar);

    // // std::cout << right << std::endl;
    // if (!right)
    // {

    //     std::cout << std::setprecision(10);
    //     std::cout << "Rij" << std::endl;
    //     TensorAlgebra::print(ricci.LL);
    //     TensorAlgebra::print(ricci2.LL);
    //     std::cout << "SCALAR" << std::endl;
    //     std::cout << ricci.scalar << std::endl;
    //     std::cout << ricci2.scalar << std::endl;
    //     MayDay::Error("BUM!");
    // }

    data_t divshift = compute_trace(d1.shift);
    data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR2(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR1(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
        covd2lapse[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(vars.A, A_UU);
    rhs.chi = advec.chi +
              (2.0 / GR_SPACEDIM) * vars.chi * (vars.lapse * vars.K - divshift);
    FOR2(i, j)
    {
        rhs.h[i][j] = advec.h[i][j] - 2.0 * vars.lapse * vars.A[i][j] -
                      (2.0 / GR_SPACEDIM) * vars.h[i][j] * divshift;
        FOR1(k)
        {
            rhs.h[i][j] +=
                vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
    }

    Tensor<2, data_t> Adot_TF;
    FOR2(i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, vars.h, h_UU);

    FOR2(i, j)
    {
        rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
                      vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR1(k)
        {
            rhs.A[i][j] +=
                vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i];
            FOR1(l)
            {
                rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
            }
        }
    }

#ifdef COVARIANTZ4
    data_t kappa1_lapse = m_params.kappa1;
#else
    data_t kappa1_lapse = m_params.kappa1 * vars.lapse;
#endif

    if (m_formulation == USE_BSSN)
    {
        rhs.Theta = 0; // ensure the Theta of CCZ4 remains at zero
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs.K = advec.K + vars.lapse * (tr_A2 + vars.K * vars.K / GR_SPACEDIM) -
                tr_covd2lapse;
        rhs.K += -2 * vars.lapse * m_cosmological_constant / (GR_SPACEDIM - 1.);
    }
    else
    {
        rhs.Theta =
            advec.Theta +
            0.5 * vars.lapse *
                (ricci.scalar - tr_A2 +
                 ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * vars.K * vars.K -
                 2 * vars.Theta * vars.K) -
            0.5 * vars.Theta * kappa1_lapse *
                ((GR_SPACEDIM + 1) + m_params.kappa2 * (GR_SPACEDIM - 1)) -
            Z_dot_d1lapse;

        rhs.Theta += -vars.lapse * m_cosmological_constant;
        rhs.K =
            advec.K +
            vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
            kappa1_lapse * GR_SPACEDIM * (1 + m_params.kappa2) * vars.Theta -
            tr_covd2lapse;
        rhs.K += -2 * vars.lapse * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                 m_cosmological_constant;
    }

    Tensor<1, data_t> Gammadot;
    FOR1(i)
    {
        Gammadot[i] = (2.0 / GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2 * m_params.kappa3 * Z_over_chi[i]) -
                           2 * vars.lapse * vars.K * Z_over_chi[i]) -
                      2 * kappa1_lapse * Z_over_chi[i];
        FOR1(j)
        {
            Gammadot[i] +=
                2 * h_UU[i][j] *
                    (vars.lapse * d1.Theta[j] - vars.Theta * d1.lapse[j]) -
                2 * A_UU[i][j] * d1.lapse[j] -
                vars.lapse * ((2 * (GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                                  h_UU[i][j] * d1.K[j] +
                              GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / vars.chi) -
                (chris.contracted[j] + 2 * m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j];

            FOR1(k)
            {
                Gammadot[i] +=
                    2 * vars.lapse * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2.shift[i][j][k] +
                    ((GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
    }

    FOR1(i) { rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i]; }

    const data_t etaDecay = 1.;

    rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                m_params.lapse_coeff * pow(vars.lapse, m_params.lapse_power) *
                    (vars.K - 2 * vars.Theta);
    FOR1(i)
    {
        rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                       m_params.shift_Gamma_coeff * vars.B[i];
        rhs.B[i] = m_params.shift_advec_coeff * advec.B[i] +
                   (1 - m_params.shift_advec_coeff) * advec.Gamma[i] +
                   Gammadot[i] - m_params.eta * etaDecay * vars.B[i];
    }
    /*
        ////////////////////////////////////
        GeometricQuantities<data_t, Vars, Diff2Vars> gq(vars, d1, d2);
        gq.set_formulation(m_formulation, m_params);
        gq.set_advection(advec);
        const auto &rhs2 = gq.get_rhs_equations();
    */
    /*
    bool right = true;
    right &= check(rhs2.chi, rhs.chi);
    FOR(i, j) { right &= check(rhs2.h[i][j], rhs.h[i][j]); }
    right &= check(rhs2.K, rhs.K);
    FOR(i, j) { right &= check(rhs2.A[i][j], rhs.A[i][j]); }
    if (m_formulation == USE_CCZ4)
        right &= check(rhs2.Theta, rhs.Theta);
    FOR(i) { right &= check(rhs2.Gamma[i], rhs.Gamma[i]); }

    // std::cout << right << std::endl;
    if (!right)
    {
        std::cout << std::setprecision(9) << std::endl;
        Tensor<1, data_t> tmp = 0.;
        FOR(i)
        {
            tmp[i] += advec.Gamma[i] +
                      (2.0 / GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2 * m_params.kappa3 * Z_over_chi[i]));
            FOR1(j)
            {
                tmp[i] += -(chris.contracted[j] +
                            2 * m_params.kappa3 * Z_over_chi[j]) *
                          d1.shift[i][j];
            }
        }
        std::cout << "CCZ4 TMP" << std::endl;
        FOR(i) { std::cout << tmp[i] << std::endl; }

        Tensor<1, data_t> tmp2 = TensorAlgebra::lie_derivative(
            advec.Gamma, vars.Gamma, d1.shift, vars.shift, divshift,
            2. / GR_SPACEDIM, {true});
        std::cout << "TMP2" << std::endl;
        FOR(i) { std::cout << tmp2[i] << std::endl; }

        Tensor<1, data_t> tmp3 = TensorAlgebra::lie_derivative(
            advec.Gamma, chris.contracted, d1.shift, vars.shift, divshift,
            2. / GR_SPACEDIM, {true});
        std::cout << "TMP3" << std::endl;
        FOR(i) { std::cout << tmp3[i] << std::endl; }

        std::cout << std::setprecision(9) << std::endl;
        std::cout << "CHI" << std::endl;
        std::cout << rhs2.chi << std::endl;
        std::cout << rhs.chi << std::endl;
        std::cout << std::endl;

        std::cout << "METRIC" << std::endl;
        FOR(i, j) { std::cout << rhs2.h[i][j] << std::endl; }
        std::cout << std::endl;
        FOR(i, j) { std::cout << rhs.h[i][j] << std::endl; }
        std::cout << std::endl;

        std::cout << "K" << std::endl;
        std::cout << rhs2.K << std::endl;
        std::cout << rhs.K << std::endl;
        std::cout << std::endl;

        std::cout << "Aij" << std::endl;
        FOR(i, j) { std::cout << rhs2.A[i][j] << std::endl; }
        std::cout << std::endl;
        FOR(i, j) { std::cout << rhs.A[i][j] << std::endl; }
        std::cout << std::endl;

        std::cout << "Theta" << std::endl;
        std::cout << rhs2.Theta << std::endl;
        std::cout << rhs.Theta << std::endl;
        std::cout << std::endl;

        std::cout << "Gamma" << std::endl;
        FOR(i) { std::cout << rhs2.Gamma[i] << std::endl; }
        std::cout << std::endl;
        FOR(i) { std::cout << rhs.Gamma[i] << std::endl; }
        std::cout << std::endl;

        MayDay::Error("BUM");
    }
    */
}

#endif /* CCZ4_IMPL_HPP_ */
