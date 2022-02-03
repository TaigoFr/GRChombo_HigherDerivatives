/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include "LorentzBoosts.hpp"
#include "parstream.H"
#include <iomanip>
#include <limits>

bool almost_equal(double value1, double value2,
                  int a_ulp /* units in the last place */)
{
    double diff = abs(value1 - value2);
    const double epsilon = std::numeric_limits<double>::epsilon();
    return (diff < std::max(epsilon * abs(value1 + value2), epsilon) * a_ulp);
}

bool check_value(double calculated_value, double analytic_value,
                 std::string value_name, int a_ulp)
{
    bool similar = almost_equal(calculated_value, analytic_value, a_ulp);

    if (!similar)
    {
        pout() << std::setprecision(16);
        pout() << "Calculated and analytic values for " << value_name
               << " do not agree\n"
               << "calculated value = " << calculated_value << "\n"
               << "analytic value = " << analytic_value << "\n"
               << "difference = " << abs(calculated_value - analytic_value)
               << "\n\n";
    }
    return similar;
}

int runLorentzBoostTest()
{
    // analytic expressions assume boost is in z.
    std::array<double, GR_SPACEDIM> v = {0.0, 0.0, 0.5};
    const double v2 = ArrayTools::norm_squared(v);
    const double lorentz_fac = LorentzBoosts::lorentz_factor(v2);
    Tensor<2, double, CH_SPACETIMEDIM> lorentz_boost_matrix =
        LorentzBoosts::matrix<double>(v);

    /*
    pout() << "lorentz_boost_matrix = \n";
    FOR_ST(i)
    {
        FOR_ST(j)
        {
            pout() << std::setw(10) << lorentz_boost_matrix[i][j] << "\t";
        }
        pout() << "\n";
    }
    */

    const double dx = 2.0;
    IntVect iv;
    iv[0] = 0;
    iv[1] = 0;
    iv[2] = 0;

    Coordinates<double> boosted_coords(iv, dx);
    // pout() << "boosted: " << boosted_coords << "\n";

    Coordinates<double> coords =
        LorentzBoosts::unboost_coords(boosted_coords, v);
    const double r = coords.get_radius();
    const double r3 = r * r * r;
    // pout() << "rest: " << coords << "\n";

    // Isotropic Schwarzschild metric
    double M = 1.0;
    double Omega = 1.0 - M / (2 * r);
    double Omega2 = Omega * Omega;
    double Psi = 1.0 + M / (2 * r);
    double Psi2 = Psi * Psi;
    double Psi3 = Psi2 * Psi;
    double Psi4 = Psi2 * Psi2;
    double Psi6 = Psi4 * Psi2;

    Tensor<2, double, CH_SPACETIMEDIM> g = 0.0;
    FOR(i) { g[i + 1][i + 1] = Psi4; };
    g[0][0] = -Omega2 / Psi2;

    // Metric partial derivatives
    Tensor<1, double> d1_Psi;
    d1_Psi[0] = -0.5 * M * coords.x / r3;
    d1_Psi[1] = -0.5 * M * coords.y / r3;
    d1_Psi[2] = -0.5 * M * coords.z / r3;

    Tensor<3, double, CH_SPACETIMEDIM> dg = 0.0;
    FOR(i, j, k)
    {
        dg[j + 1][k + 1][i + 1] =
            4.0 * Psi3 * static_cast<double>(TensorAlgebra::delta(j, k)) *
            d1_Psi[i];
    }
    FOR(i) { dg[0][0][i + 1] = 4.0 * Omega * d1_Psi[i] / Psi3; }

    auto g_boosted = LorentzBoosts::boost_LL(g, v);
    auto dg_boosted = LorentzBoosts::boost_LLL(dg, v);
    auto adm_vars_boosted =
        TensorAlgebra::adm_vars_from_metric_ST(g_boosted, dg_boosted);

    double analytic_g_zz =
        lorentz_fac * lorentz_fac * (Psi4 - v2 * Omega2 / Psi2);
    double analytic_lapse =
        Omega * Psi2 / (lorentz_fac * sqrt(Psi6 - v2 * Omega2));
    double analytic_shift_z = -v[2] * (Psi6 - Omega2) / (Psi6 - v2 * Omega2);
    double analytic_K_zz =
        -lorentz_fac * lorentz_fac * v[2] *
        (2.0 * Psi6 * (Psi + 2.0 * Omega) - 2.0 * v2 * Omega2) * d1_Psi[2] /
        (Psi * Psi4 * sqrt(Psi6 - v2 * Omega2));

    int status = 0;
    constexpr int ulp = 10;
    status |=
        !check_value(adm_vars_boosted.lapse, analytic_lapse, "lapse", ulp);
    status |= !check_value(adm_vars_boosted.shift[2], analytic_shift_z,
                           "shift_z", ulp);
    status |= !check_value(g_boosted[3][3], analytic_g_zz, "g_zz", ulp);
    status |=
        !check_value(adm_vars_boosted.K_LL[2][2], analytic_K_zz, "K_zz", ulp);

    return status;
}

int main()
{
    int status = runLorentzBoostTest();

    if (status == 0)
        pout() << "LorentzBoost test passed." << std::endl;
    else
        pout() << "LorentzBoost test failed with return code " << status
               << std::endl;

    return status;
}
