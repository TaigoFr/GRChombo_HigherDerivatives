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
                               int a_formulation, double a_G_Newton)
    : MatterCCZ4RHS<matter_t, gauge_t, deriv_t>(
          a_matter, a_params, a_dx, a_sigma, a_formulation, a_G_Newton),
      m_diffusion_params(a_diffusion_params), m_dt(a_dt)
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

    GeometricQuantities<data_t, Vars, Diff2Vars, gauge_t> gq(matter_vars, d1,
                                                             d2);
    gq.set_advection_and_gauge(advec, this->m_gauge);
    gq.set_formulation(this->m_formulation, this->m_params);
    gq.set_cosmological_constant(this->m_cosmological_constant);

    // add evolution of matter fields themselves
    // done before to avoid getting contributions from the EM-tensor
    Vars<data_t> matter_rhs;
    this->my_matter.add_matter_rhs(matter_rhs, gq);

    const auto emtensor = this->my_matter.compute_emtensor(gq);
    gq.set_em_tensor(emtensor, this->m_G_Newton);

    gq.compute_rhs_equations(matter_rhs);

    add_diffusion_terms(matter_rhs, matter_vars, d1, d2);

    // Add dissipation to all terms
    this->m_deriv.add_dissipation(matter_rhs, current_cell, this->m_sigma);

    // Write the matter_rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

/// Adds diffusion terms to the rhs of the GHC equations
template <class matter_t, class gauge_t, class deriv_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void MatterCCZ4RHSWithDiffusion<matter_t, gauge_t, deriv_t>::
    add_diffusion_terms(
        vars_t<data_t> &rhs, //!< Reference to the variables into which the
                             //! output right hand side is written
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        diff2_vars_t<Tensor<2, data_t>> &d2) const
{
    using namespace TensorAlgebra;

    // based on arxiv:1512.04532v2 (in Supplemental Material)

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
    auto chi_above_cutoff =
        sigmoid(vars.chi, -m_diffusion_params.chiCutoff_width,
                m_diffusion_params.chiCutoff);
    diffCoeffSafe *= chi_above_cutoff;

    double diffusion_maximum = m_diffusion_params.diffusion_maximum;
    rhs.enum_mapping([&d2, diffusion_maximum, diffCoeffSafe](const int &ivar,
                                                             data_t &rhs_var) {
        d2.enum_mapping([&rhs_var, ivar, diffusion_maximum, diffCoeffSafe](
                            const int &jvar, const Tensor<2, data_t> &d2_var) {
            if (ivar == jvar)
            {
                data_t space_laplace = 0.;
                FOR(k) { space_laplace += d2_var[k][k]; }

                FOR(i, j)
                {
                    rhs_var +=
                        simd_max(-diffusion_maximum,
                                 simd_min(diffusion_maximum,
                                          diffCoeffSafe * space_laplace));
                }
            }
        });
    });
}

#endif /* MATTERCCZ4RHSWITHDIFFUSION_IMPL_HPP_ */
