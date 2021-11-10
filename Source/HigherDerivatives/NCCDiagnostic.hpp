/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NCCDIAGNOSTIC
#define NCCDIAGNOSTIC

#include "C2EFT.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GeometricQuantities.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"

//! Computes Null Convergence Condition from Stress-Energy Tensor
template <class System> class NCCDiagnostic
{
    // Use the variable definitions in MatterCCZ4RHS
    template <class data_t>
    using Vars = typename MatterCCZ4RHS<C2EFT<System>>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename MatterCCZ4RHS<C2EFT<System>>::template Diff2Vars<data_t>;

  public:
    NCCDiagnostic(const C2EFT<System> &a_matter, double m_dx, int a_formulation,
                  const CCZ4_params_t<> &a_ccz4_params,
                  const std::array<double, CH_SPACEDIM> &a_center,
                  double G_Newton, int a_NCC_plus, int a_NCC_minus,
                  int a_NCC_Z4_plus = -1, int a_NCC_Z4_minus = -1)
        : m_matter(a_matter), m_formulation(a_formulation),
          m_ccz4_params(a_ccz4_params), m_deriv(m_dx), m_center(a_center),
          m_G_Newton(G_Newton), m_NCC_plus(a_NCC_plus),
          m_NCC_minus(a_NCC_minus), m_NCC_Z4_plus(a_NCC_Z4_plus),
          m_NCC_Z4_minus(a_NCC_Z4_minus)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec =
            m_deriv.template advection<Vars>(current_cell, vars.shift);

        Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);

        // make gauge
        MovingPunctureGauge gauge(m_ccz4_params);

        GeometricQuantities<data_t, Vars, Diff2Vars, MovingPunctureGauge> gq(
            vars, d1, d2, "NCCDiagnostic::compute");

        gq.set_formulation(m_formulation, m_ccz4_params);
        gq.set_advection_and_gauge(advec,
                                   gauge); // needed for 'compute_emtensor'
        gq.set_coordinates(coords);

        const auto emtensor = m_matter.compute_emtensor(gq);
        gq.set_em_tensor(emtensor, m_G_Newton);

        Tensor<1, data_t, GR_SPACEDIM + 1> null_vec_plus, null_vec_minus;
        get_null_radial_vector(null_vec_plus, null_vec_minus, vars, coords);

        // check it is null - already checked and it was good!
        // auto &g = gq.get_metric_ST();
        // data_t norm4 =
        // TensorAlgebra::compute_dot_product(null_vec_plus, null_vec_plus, g);
        // CH_assert(simd_compare_lt_any(abs(norm4), 1.e-10));

        auto &em_tensor_ST = gq.get_em_tensor_ST();

        if (m_NCC_plus >= 0)
        {
            data_t NCC_plus = TensorAlgebra::compute_dot_product(
                null_vec_plus, null_vec_plus, em_tensor_ST);
            current_cell.store_vars(NCC_plus, m_NCC_plus);
        }
        if (m_NCC_minus >= 0)
        {
            data_t NCC_minus = TensorAlgebra::compute_dot_product(
                null_vec_minus, null_vec_minus, em_tensor_ST);
            current_cell.store_vars(NCC_minus, m_NCC_minus);
        }

        if (m_NCC_Z4_plus >= 0 || m_NCC_Z4_minus >= 0)
        {
            auto &ricci_ST = gq.get_ricci_ST();

            if (m_NCC_Z4_plus >= 0)
            {

                data_t NCC_Z4_plus = TensorAlgebra::compute_dot_product(
                    null_vec_plus, null_vec_plus, ricci_ST);
                current_cell.store_vars(NCC_Z4_plus, m_NCC_Z4_plus);
            }
            if (m_NCC_Z4_minus >= 0)
            {

                data_t NCC_Z4_minus = TensorAlgebra::compute_dot_product(
                    null_vec_minus, null_vec_minus, ricci_ST);

                current_cell.store_vars(NCC_Z4_minus, m_NCC_Z4_minus);
            }
        }
    }

    // compute null radial vector by summing the timelike normal (usual n^a)
    // to a radial spacelike vector
    template <class data_t>
    void
    get_null_radial_vector(Tensor<1, data_t, GR_SPACEDIM + 1> &null_vec_plus,
                           Tensor<1, data_t, GR_SPACEDIM + 1> &null_vec_minus,
                           const Vars<data_t> &vars,
                           const Coordinates<data_t> &coords) const
    {
        // add normal vector first
        data_t lapse = vars.lapse;
        null_vec_plus[0] = 1. / lapse;
        null_vec_minus[0] = null_vec_plus[0];
        FOR(i)
        {
            null_vec_plus[i + 1] = vars.shift[i] / lapse;
            null_vec_minus[i + 1] = null_vec_plus[i + 1];
        }

        // add spatial vector
        Tensor<1, data_t> spatial_vec = {coords.x, coords.y, coords.z};
        // normalize it first
        data_t norm = 0.;
        FOR(i, j)
        {
            norm += vars.h[i][j] / vars.chi * spatial_vec[i] * spatial_vec[j];
        }
        norm = sqrt(norm);

        FOR(i)
        {
            null_vec_plus[i + 1] += spatial_vec[i] / norm;
            null_vec_minus[i + 1] -= spatial_vec[i] / norm;
        }
    }

  protected:
    const C2EFT<System> &m_matter;
    int m_formulation;
    const CCZ4_params_t<> &m_ccz4_params;
    FourthOrderDerivatives m_deriv;
    std::array<double, CH_SPACEDIM> m_center;
    double m_G_Newton;
    int m_NCC_plus, m_NCC_minus, m_NCC_Z4_plus, m_NCC_Z4_minus;
};

#endif /* NCCDIAGNOSTIC */
