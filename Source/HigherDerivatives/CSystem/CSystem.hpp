/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CSYSTEM_HPP_
#define CSYSTEM_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

#include "GeometricQuantities.hpp"

class CSystem
{
  public:
    struct params_t
    {
        double tau;   // for the evolution equation of C
        double sigma; // for the wave operator in the evolution of C
        bool use_only_time_derivatives; // make the physical C a static solution
        bool rescale_tau_sigma_by_lapse; // for when using only time derivatives
        bool add_advection;              // for when using only time derivatives
    };

    //!  Constructor of class CSystem, inputs are the matter parameters.
    CSystem(params_t a_params) : m_params(a_params) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t C;
        data_t dCdt;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_C, C);
            VarsTools::define_enum_mapping(mapping_function, c_dCdt, dCdt);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t C;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_C, C);
        }
    };

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void compute_C(
        data_t &C, Tensor<1, data_t, CH_SPACEDIM + 1> &d1_C,
        Tensor<2, data_t, CH_SPACEDIM + 1> &d2_C,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void compute_Riemann(
        Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LLLU,
        Tensor<4, data_t, CH_SPACEDIM + 1> &riemann_LULU,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< value of the RHS for all vars
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    template <class data_t, template <typename> class rhs_vars_t,
              template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void add_diffusion_terms(
        rhs_vars_t<data_t> &rhs, //!< Reference to the variables into which the
                                 //! output right hand side is written
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
        data_t diffCoeffSafe) const;

  private:
    params_t m_params;
};

#include "CSystem.impl.hpp"

#endif /* CSYSTEM_HPP_ */
