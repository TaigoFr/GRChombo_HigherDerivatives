/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef C2EFT_HPP_
#define C2EFT_HPP_

#include "GeometricQuantities.hpp"
#include "MatterCCZ4RHSWithDiffusion.hpp"
#include "Tensor.hpp"

template <class System> class C2EFT
{
  public:
    struct params_t : public ah_chi_params_t
    {
        double epsilon;
        double epsilon_final;        
        double time_epsilon;        

        // weak field parameters
        double weak_field_threshold;
        double weak_field_width;
    };

    template <class data_t> using Vars = typename System::template Vars<data_t>;
    template <class data_t>
    using Diff2Vars = typename System::template Diff2Vars<data_t>;

    //!  Constructor of class C2EFT, inputs are the matter parameters.
    C2EFT(System &a_system, params_t a_params, bool apply_weak_field)
        : m_system(a_system), m_params(a_params),
          m_apply_weak_field(apply_weak_field)
    {
    }

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    emtensor_t<data_t> compute_emtensor(
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void compute_emtensor_4D(
        Tensor<2, data_t, CH_SPACETIMEDIM> &Tmn,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< value of the RHS for all vars
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    data_t weak_field_var(
        const emtensor_t<data_t> &emtensor,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    static data_t weak_field_condition(
        const data_t &weak_field_var,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
        const C2EFT<System>::params_t &pm);

    template <class data_t, template <typename> class rhs_vars_t,
              template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void add_diffusion_terms(
        rhs_vars_t<data_t> &rhs, //!< Reference to the variables into which the
                                 //! output right hand side is written
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
        data_t diffCoeffSafe) const;

    params_t m_params;

  private:
    System &m_system;
    bool m_apply_weak_field;

    // output is <10^{-k} for x>t(1+k+w) and >1-10^{-k} for x<t(1-k*w)
    template <class data_t>
    inline static data_t sigmoid(data_t x, double width, double threshold)
    {
        return 1. / (1. + pow((data_t)10., (x / threshold - 1.) / width));
    }
};

#include "C2EFT.impl.hpp"

#endif /* C2EFT_HPP_ */
