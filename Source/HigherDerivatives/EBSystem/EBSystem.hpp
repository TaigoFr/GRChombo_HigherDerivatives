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

#include "C2EFT.hpp" //added for Box_transition

class EBSystem
{
  public:
    struct params_t
    {
        double tau; // for the evolution equation of Eij and Bij
        bool rescale_tau_by_lapse;
        int rescale_sigma_by_lapse; // (0, 1 or 2 for lapse or lapse^2
                                    // rescaling)

        int version; // v1 is the 1st order eq, v2 is the 2nd order eq.,
                     // v3 is laplacian based equation

        double sigma;           // only for v2 or v3
        int advection_type;     // only for v2
                                // (0, 1 or 2 for simple advection or Luis'
                                // advection proposal)
                                // or in system v3 (1 or 0 for on / off)
        double advection_coeff; // for advection_type 1 and 2 in v2

        bool Box_transition; // Could be use to active a transition for the
                             // advection velocities

        bool use_last_index_raised;

        bool use_tau_radial_decay;
        double tau_asymptotic;   // if 'use_tau_radial_decay'
        double tau_decay_length; // if 'use_tau_radial_decay'
    };

    //!  Constructor of class EBSystem, inputs are the matter parameters.
    EBSystem(params_t a_params) : m_params(a_params) {}

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        // adding 'ij' because B already exists from Gamma-driver
        Tensor<2, data_t> Eij;
        Tensor<2, data_t> Bij;
        // v1: need to store physical ones to get their spatial derivatives
        // v2: need to store the first order variables for the 2nd order eq.
        Tensor<2, data_t> Eaux;
        Tensor<2, data_t> Baux;

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
                mapping_function, GRInterval<c_Eaux11, c_Eaux33>(), Eaux);
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_Baux11, c_Baux33>(), Baux);
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
    void
    compute_C(data_t &C, Tensor<1, data_t, CH_SPACETIMEDIM> &d1_C,
              Tensor<2, data_t, CH_SPACETIMEDIM> &d2_C,
              GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
              const C2EFT<EBSystem>::params_t &pm) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t>
    void compute_Riemann(
        Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LLLU,
        Tensor<4, data_t, CH_SPACETIMEDIM> &riemann_LULU,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq) const;

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< value of the RHS for all vars
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
        const C2EFT<EBSystem>::params_t &pm) const;

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

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t, class gauge_t,
              template <typename> class rhs_vars_t>
    void compute_d2_Eij_and_Bij(
        Tensor<2, Tensor<2, data_t, CH_SPACETIMEDIM>> &d2_Eij,
        Tensor<2, Tensor<2, data_t, CH_SPACETIMEDIM>> &d2_Bij,
        GeometricQuantities<data_t, vars_t, diff2_vars_t, gauge_t> &gq,
        rhs_vars_t<data_t> &rhs,
        Tensor<2, Tensor<1, data_t, CH_SPACETIMEDIM>> &d1_h) const;

    template <class data_t> data_t tau_from_radius(const data_t &radius) const
    {
        return simd_min((data_t)m_params.tau,
                        (m_params.tau - m_params.tau_asymptotic) *
                                m_params.tau_decay_length / radius +
                            m_params.tau_asymptotic);
    }
};

#include "EBSystem.impl.hpp"

#endif /* EBSYSTEM_HPP_ */
