/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCCZ4RHSWITHDIFFUSION_HPP_)
#error                                                                         \
    "This file should only be included through MatterCCZ4RHSWithDiffusion.hpp"
#endif

#ifndef MATTERCCZ4RHSWITHDIFFUSION_IMPL_HPP_
#define MATTERCCZ4RHSWITHDIFFUSION_IMPL_HPP_

template <class matter_t, class gauge_t, class deriv_t>
inline MatterCCZ4RHSWithDiffusion<matter_t, gauge_t, deriv_t>::
    MatterCCZ4RHSWithDiffusion(matter_t a_matter, params_t a_params,
                               diffusion_params_t a_diffusion_params,
                               double a_dx, double a_dt, double a_sigma,
                               const std::array<double, CH_SPACEDIM> a_center,
                               int a_formulation, double a_G_Newton)
    : MatterCCZ4RHS<matter_t, gauge_t, deriv_t>(
          a_matter, a_params, a_dx, a_sigma, a_formulation, a_G_Newton),
      m_diffusion_params(a_diffusion_params), m_dt(a_dt), m_center(a_center)
{
}

template <class matter_t, class gauge_t, class deriv_t>
template <class data_t>
void MatterCCZ4RHSWithDiffusion<matter_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{
    CH_TIME("MatterCCZ4::compute");

    // copy data from chombo gridpoint into local variables
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    // d2 can't be const to do an 'enum_mapping' in 'add_diffusion_terms'
    auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        this->m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

    const Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx,
                                     m_center);

    GeometricQuantities<data_t, Vars, Diff2Vars, gauge_t> gq(
        matter_vars, d1, d2, "MatterCCZ4RHSWithDiffusion::compute");

    gq.set_advection_and_gauge(advec, this->m_gauge);
    gq.set_formulation(this->m_formulation, this->m_params);
    gq.set_cosmological_constant(this->m_cosmological_constant);
    gq.set_coordinates(coords);

    // add evolution of matter fields themselves
    // done before to avoid getting contributions from the EM-tensor
    Vars<data_t> matter_rhs;
    this->my_matter.add_matter_rhs(matter_rhs, gq);

    const auto emtensor = this->my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, this->m_G_Newton);

    gq.compute_rhs_equations(matter_rhs);

    data_t diffCoeffSafe = add_diffusion_terms(matter_rhs, gq);
    this->my_matter.add_diffusion_terms(matter_rhs, gq, diffCoeffSafe);

    // Add dissipation to all terms
    this->m_deriv.add_dissipation(matter_rhs, current_cell, this->m_sigma);

    // Write the matter_rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

template <class matter_t, class gauge_t, class deriv_t>
template <class data_t, template <typename> class rhs_vars_t,
          template <typename> class vars_t,
          template <typename> class diff2_vars_t>
data_t
MatterCCZ4RHSWithDiffusion<matter_t, gauge_t, deriv_t>::add_diffusion_terms(
    rhs_vars_t<data_t> &rhs,
    GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const
{
    // based on arxiv:1512.04532v2 (in Supplemental Material)

    const auto &vars = gq.get_vars();
    const auto &d1 = gq.get_d1_vars();
    const auto &d2 = gq.get_d2_vars();
    const auto &h_UU = gq.get_h_UU();

    data_t diffCoeff = 0;
    FOR(i, j, k) { diffCoeff += pow(d1.h[i][j][k], 2); }

    diffCoeff = pow(sqrt(2. * diffCoeff / ((double)GR_SPACEDIM + 1.) *
                         ((double)GR_SPACEDIM)),
                    m_diffusion_params.lapidusPower) *
                pow(this->m_deriv.m_dx, m_diffusion_params.lapidusPower - 1);
    diffCoeff = 0.5 * m_diffusion_params.lapidusCoeff *
                pow(this->m_deriv.m_dx, 2) * diffCoeff;
    data_t diffCoeffSafe = simd_min(
        diffCoeff, m_diffusion_params.diffCFLFact * pow(this->m_deriv.m_dx, 2) /
                       m_dt); //!< Do not allow the diffusion coefficient to
                              //!< violate the Courant condition

    // Introduce a smooth cutoff:
    auto chi_above_cutoff = sigmoid(vars.chi, m_diffusion_params.chi_width,
                                    m_diffusion_params.chi_threshold);
    diffCoeffSafe *= chi_above_cutoff;

    data_t space_laplace_lapse = 0.;
    data_t space_laplace_chi = 0.;
    Tensor<1, data_t> space_laplace_shift = {0.};
    Tensor<2, data_t> space_laplace_h = {0.};

    FOR(k)
    {
        space_laplace_lapse += d2.lapse[k][k];
        space_laplace_chi += d2.chi[k][k];

        FOR(i)
        {
            space_laplace_shift[i] += d2.shift[i][k][k];
            FOR(j) { space_laplace_h[i][j] += d2.h[i][j][k][k]; }
        }
    }

    data_t tr_space_laplace_h = 0.;
    FOR(i, j) { tr_space_laplace_h += h_UU[i][j] * space_laplace_h[i][j]; }

    rhs.lapse += diffCoeffSafe * space_laplace_lapse;
    rhs.chi += diffCoeffSafe * space_laplace_chi;

    FOR(i)
    {
        rhs.shift[i] += diffCoeffSafe * space_laplace_shift[i];
        FOR(j)
        {
            rhs.h[i][j] += diffCoeffSafe *
                           (space_laplace_h[i][j] -
                            tr_space_laplace_h * vars.h[i][j] / GR_SPACEDIM);
        }
    }

    return diffCoeffSafe;
}

#endif /* MATTERCCZ4RHSWITHDIFFUSION_IMPL_HPP_ */
