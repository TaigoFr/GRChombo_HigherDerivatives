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

#include "C2EFT.hpp"
#include "CSystem.hpp"
#include "MatterCCZ4.hpp"

#include "SetupFunctions.hpp" // just to avoid some undefined references of static Derivative::LOCAL and other vars

template <class data_t> using Vars = MatterCCZ4<C2EFT<CSystem>>::Vars<data_t>;
template <class data_t>
using Diff2Vars = MatterCCZ4<C2EFT<CSystem>>::Diff2Vars<data_t>;

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
    int formulation = CCZ4::USE_CCZ4;

    CCZ4::params_t ccz4_params;

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

    C2EFT<CSystem>::params_t hd_params;
    hd_params.epsilon = 0.1;

    CSystem::params_t c_params;
    c_params.tau = 0.1;
    c_params.sigma = 1.;

    int failed = 0;

    //////////////////////////////////////////////////////////////////
    // SET VARIABLES

    Vars<double> vars;
    Vars<Tensor<1, double>> d1;
    Diff2Vars<Tensor<2, double>> d2;

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
                FOR1(dir) { var += vars.shift[dir] * d1_var[dir]; }
            }
        });
    });

    GeometricQuantities<double, Vars, Diff2Vars> gq(vars, d1, d2);
    gq.set_formulation(formulation, ccz4_params);
    gq.set_advection(advec);

    // make EM Tensor
    CSystem Csystem(c_params);
    C2EFT<CSystem> c2eft(Csystem, hd_params, false);
    const auto em_tensor_def = c2eft.compute_emtensor(gq);

    //////////////////////////////////////////////////////////////////
    // SET MATHEMATICA CALCULATIONS

    emtensor_t<double> em_tensor;

#include "randomData.hpp" //Including the auto generated file with calculations

    //////////////////////////////////////////////////////////////////
    // CHECK
    failed |= relative_error(em_tensor_def.rho, em_tensor.rho, "emtensor.rho",
                             2.e-12);
    failed |=
        relative_error(em_tensor_def.Si, em_tensor.Si, "emtensor.Si", 2.e-12);
    failed |= relative_error(em_tensor_def.Sij, em_tensor.Sij, "emtensor.Sij",
                             2.e-12);
    failed |=
        relative_error(em_tensor_def.S, em_tensor.S, "emtensor.S", 2.e-12);

    return failed;
}

int main(int argc, char *argv[])
{
    // mainSetup(argc, argv);

    int failed = runTest(argc, argv);

    if (failed == 0)
        pout() << "CSystem test passed." << std::endl;
    else
        pout() << "CSystem test failed with return code " << failed
               << std::endl;

    // mainFinalize();
    return failed;
}
