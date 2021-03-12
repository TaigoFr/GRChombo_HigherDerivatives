/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIFFUSIONDIAGNOSTIC
#define DIFFUSIONDIAGNOSTIC

#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHSWithDiffusion.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

template <class matter_t>
class DiffusionDiagnostic
    : MatterCCZ4RHSWithDiffusion<matter_t, MovingPunctureGauge,
                                 FourthOrderDerivatives>
{
  public:
    // Use this alias for the same template instantiation as this class
    using MatterCCZ4 = MatterCCZ4RHSWithDiffusion<matter_t, MovingPunctureGauge,
                                                  FourthOrderDerivatives>;

    using params_t = typename MatterCCZ4::params_t;

    template <class data_t>
    using Vars = typename MatterCCZ4::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars = typename MatterCCZ4::template Diff2Vars<data_t>;

  public:
    DiffusionDiagnostic(matter_t &a_matter, params_t a_params,
                        diffusion_params_t a_diffusion_params, double a_dx,
                        double a_dt, double a_sigma, int a_formulation,
                        double a_G_Newton, int a_diffusion_var,
                        int a_diffusion_rhs_var,
                        int a_det_conformal_metric = -1)
        : MatterCCZ4RHSWithDiffusion<matter_t, MovingPunctureGauge,
                                     FourthOrderDerivatives>(
              a_matter, a_params, a_diffusion_params, a_dx, a_dt, a_sigma,
              a_formulation, a_G_Newton),
          m_diffusion_var(a_diffusion_var),
          m_diffusion_rhs_var(a_diffusion_rhs_var),
          m_det_conformal_metric(a_det_conformal_metric)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        GeometricQuantities<data_t, Vars, Diff2Vars, MovingPunctureGauge> gq(
            vars, d1, d2);
        gq.set_advection_and_gauge(advec, this->m_gauge);
        gq.set_formulation(this->m_formulation, this->m_params);

        Vars<data_t> rhs;

        // not needed as we are only looking at the RHS of the metric for the
        // moment my_matter.add_matter_rhs(rhs, gq); const auto emtensor =
        // my_matter.compute_emtensor(gq); gq.set_em_tensor(emtensor,
        // this->m_G_Newton);

        gq.compute_rhs_equations(rhs);

        Vars<data_t> diffusion;
        this->template add_diffusion_terms(diffusion, gq);

        // for the moment only saving for vars.chi, just to see how it looks
        current_cell.store_vars(diffusion.chi, m_diffusion_var);
        current_cell.store_vars(rhs.chi, m_diffusion_rhs_var);

        if (m_det_conformal_metric >= 0)
        {
            data_t det_h = TensorAlgebra::compute_determinant(vars.h);
            current_cell.store_vars(det_h, m_det_conformal_metric);
        }
    }

  protected:
    int m_diffusion_var, m_diffusion_rhs_var, m_det_conformal_metric;
};

#endif /* DIFFUSIONDIAGNOSTIC */
