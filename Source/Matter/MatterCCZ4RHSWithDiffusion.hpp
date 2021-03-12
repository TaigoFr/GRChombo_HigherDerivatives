/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERCCZ4RHSWITHDIFFUSION_HPP_
#define MATTERCCZ4RHSWITHDIFFUSION_HPP_

#include "MatterCCZ4RHS.hpp"

struct diffusion_params_t
{
    double diffCFLFact;  //!< diffusion terms CFL factor cutoff
    double lapidusCoeff; //!< Lapidus coefficient
    double lapidusPower; //!< Lapidus power
    double chiCutoff;    //!< Cut off for diffusion terms
    double chiCutoff_width;
};

template <class matter_t, class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives>
class MatterCCZ4RHSWithDiffusion
    : public MatterCCZ4RHS<matter_t, gauge_t, deriv_t>
{
  public:
    // Use this alias for the same template instantiation as this class
    using MatterCCZ4 = MatterCCZ4RHS<matter_t, gauge_t, deriv_t>;

    using params_t = typename MatterCCZ4::params_t;

    template <class data_t>
    using Vars = typename MatterCCZ4::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars = typename MatterCCZ4::template Diff2Vars<data_t>;

    /// Constructor
    MatterCCZ4RHSWithDiffusion(matter_t a_matter, params_t a_params,
                               diffusion_params_t a_diffusion_params,
                               double a_dx, double a_dt, double a_sigma,
                               int a_formulation = CCZ4::USE_CCZ4,
                               double a_G_Newton = 1.0);

    /// Compute function
    /** This function orchestrates the calculation of the rhs for one specific
     * grid cell. This function is called by the BoxLoops::loop for each grid
     * cell; there should rarely be a need to call it directly.
     */
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    diffusion_params_t m_diffusion_params;
    double m_dt;

    template <class data_t, template <typename> class rhs_vars_t,
              template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    data_t add_diffusion_terms(
        rhs_vars_t<data_t> &rhs, //!< Reference to the variables into which the
                                 //! output right hand side is written
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    // output is <10^{-k} for x>t(1+k+w) and >1-10^{-k} for x<t(1-k*w)
    template <class data_t>
    inline static data_t sigmoid(data_t x, double width, double threshold)
    {
        return 1. / (1. + pow((data_t)10., (x / threshold - 1.) / width));
    }
};

#include "MatterCCZ4RHSWithDiffusion.impl.hpp"

#endif /* MATTERCCZ4RHSWITHDIFFUSION_HPP_ */
