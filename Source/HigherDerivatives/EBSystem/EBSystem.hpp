/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EBSYSTEM_HPP_
#define EBSYSTEM_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

#include "GeometricQuantities.hpp"

class EBSystem
{
  public:
    struct params_t
    {
        double tau; // for the evolution equation of Eij and Bij
    };

    //!  Constructor of class EBSystem, inputs are the matter parameters.
    EBSystem(params_t a_params) : m_params(a_params) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        // adding 'ij' because B already exists from Gamma-driver
        Tensor<2, data_t> Eij;
        Tensor<2, data_t> Bij;
        // need to store physical ones to get their spatial derivatives
        Tensor<2, data_t> Ephys;
        Tensor<2, data_t> Bphys;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            // Symmetric 2-tensors
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_E11, c_E33>(), Eij);
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_B11, c_B33>(), Bij);

            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_Ephys11, c_Ephys33>(), Ephys);
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_Bphys11, c_Bphys33>(), Bphys);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        // adding 'ij' because B already exists from Gamma-driver
        Tensor<2, data_t> Eij;
        Tensor<2, data_t> Bij;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            // Symmetric 2-tensors
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_E11, c_E33>(), Eij);
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_B11, c_B33>(), Bij);
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

  private:
    params_t m_params;
};

#include "EBSystem.impl.hpp"

#endif /* EBSYSTEM_HPP_ */
