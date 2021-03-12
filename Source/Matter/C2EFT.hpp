/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef C2EFT_HPP_
#define C2EFT_HPP_

#include "GeometricQuantities.hpp"
#include "Tensor.hpp"

template <class System> class C2EFT
{
  public:
    struct params_t
    {
        double epsilon;

        // weak field parameters
        double chi_threshold;
        double chi_width;
        double weak_field_threshold;
        double weak_field_width;
    };

    template <class data_t> using Vars = typename System::template Vars<data_t>;
    template <class data_t>
    using Diff2Vars = typename System::template Diff2Vars<data_t>;

    //!  Constructor of class C2EFT, inputs are the matter parameters.
    C2EFT(System a_system, params_t a_params)
        : m_system(a_system), m_params(a_params)
    {
    }

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    emtensor_t<data_t>
    compute_emtensor(GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq,
                     bool apply_weak_field = true) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void compute_emtensor_4D(
        Tensor<2, data_t, CH_SPACEDIM + 1> &Tmn,
        GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const;

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< value of the RHS for all vars
        GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    data_t weak_field_condition(
        const data_t &weak_field_var,
        GeometricQuantities<data_t, vars_t, diff2_vars_t> &gq) const;

  private:
    System m_system;
    params_t m_params;

    template <class data_t>
    inline static data_t sigmoid(data_t x, double width, double threshold)
    {
        return 1. / (1. + exp(2. / width * (x / threshold - 1.)));
    }
};

#include "C2EFT.impl.hpp"

#endif /* C2EFT_HPP_ */
