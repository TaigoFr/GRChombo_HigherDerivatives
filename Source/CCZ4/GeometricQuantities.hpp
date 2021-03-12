/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GEOMETRICQUANTITIES_HPP_
#define GEOMETRICQUANTITIES_HPP_

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
template <class data_t, int size = CH_SPACEDIM> struct emtensor_t
{
    Tensor<2, data_t, size> Sij; //!< S_ij = T_ij
    Tensor<1, data_t, size> Si;  //!< S_i = T_ia_n^a
    data_t S;                    //!< S = S^i_i
    data_t rho;                  //!< rho = T_ab n^a n^b
};

template <class data_t, int size = CH_SPACEDIM> struct ricci_t
{
    Tensor<2, data_t, size> LL; // Ricci with two indices down
    data_t scalar;              // Ricci scalar
};

#include "CCZ4.hpp" // need 'formulations'

//!  Calculates the several spatial and spacetime geometric quantities
/*!
   This class calculates geometric quantities such as the 3-Ricci, the
   4-Weyl, the Kretschmann scalar, etc.. It's made such that calculations
   are minimized, re-using previously calculated quantities whenever
   possible. To be used inside compute functions.
*/
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
class GeometricQuantities
{
  public:
    using Vars = vars_t<data_t>;
    using Diff1Vars = vars_t<Tensor<1, data_t>>;
    using Diff2Vars = diff2_vars_t<Tensor<2, data_t>>;

    // enum Output
    // { // types of output
    //     Kretschmann,
    //     Riemann_max,
    //     Weyl_squared,
    //     Ricci_scalar,
    //     Ricci_max
    // };

    //! Constructor of class GeometricQuantities
    GeometricQuantities();
    GeometricQuantities(const Vars &a_vars, const Diff1Vars &a_d1_vars,
                        const Diff2Vars &a_d2_vars);
    ~GeometricQuantities();

    // disable copy and move semantics, to prevent copying pointers
    GeometricQuantities(const GeometricQuantities &) = delete;
    GeometricQuantities &operator=(const GeometricQuantities &) = delete;

    //! The compute member which calculates added geometric quantities at each
    //! point on the grid
    // template <class data_t> void compute(Cell<data_t> current_cell) const;

    // allow the user not to set d2 if not needed.
    void set_vars(const Vars &a_vars);
    void set_d1_vars(const Diff1Vars &a_d1_vars);
    void set_d2_vars(const Diff2Vars &a_d2_vars);
    void set_advection(const Vars &a_advection);
    void set_formulation(int formulation, const CCZ4::params_t &a_ccz4_params);
    void set_em_tensor(const emtensor_t<data_t> &a_em_tensor, double G_Newton);
    void set_cosmological_constant(double cosmological_constant);

    void set_all_vars(const Vars &a_vars, const Diff1Vars &a_d1_vars,
                      const Diff2Vars &a_d2_vars);

    //////// SET BY USER ////////
    // vars, d1, d2 and emtensor are assumed to stay in scope
    const double get_formulation() const;
    const double get_cosmological_constant() const;
    const Vars &get_vars() const;
    const Diff1Vars &get_d1_vars() const;
    const Diff2Vars &get_d2_vars() const;
    const Vars &get_advection() const;
    const emtensor_t<data_t> &get_em_tensor() const;

    //////// spatial ////////
    // spatial conformal
    const Tensor<2, data_t> &get_h_UU();
    const chris_t<data_t> &get_chris();
    const Tensor<1, data_t> &get_Z_U_conformal();
    const Tensor<2, data_t> &get_d1_chris_contracted();
    const Tensor<2, data_t> &get_covd_chi_conformal();
    const Tensor<4, data_t> &get_riemann_conformal_LLLL();
    const Tensor<2, data_t> &get_A_LU();
    const Tensor<2, data_t> &get_A_UU();
    const data_t &get_tr_A2();

    // spatial non-conformal
    const Tensor<1, data_t> &get_shift_L();
    const Tensor<2, data_t> &get_metric_spatial();
    const Tensor<2, data_t> &get_metric_UU_spatial();
    const Tensor<2, data_t> &get_extrinsic_curvature();
    const Tensor<3, data_t> &get_chris_spatial();
    const Tensor<1, data_t> &get_Z_U();
    const Tensor<1, data_t> &get_Z();
    const Tensor<2, data_t> &get_covd_Z();
    const Tensor<3, data_t> &get_d1_extrinsic_curvature();
    const Tensor<3, data_t> &get_covd_extrinsic_curvature();
    const Tensor<3, data_t> &get_levi_civita_spatial();     // LLL
    const Tensor<3, data_t> &get_levi_civita_spatial_LUU(); // LUU
    const Tensor<2, data_t> &get_covd_lapse();
    const ricci_t<data_t> &get_ricci_qDZ(int q);
    const ricci_t<data_t> &get_ricci(); // original formula
    const ricci_t<data_t> &get_ricci_1DZ();
    // uses evolved variables, not calculated Gammas
    const ricci_t<data_t> &get_ricci_2DZ();
    const Tensor<2, data_t> &get_weyl_magnetic_part();
    const Tensor<4, data_t> &get_riemann_spatial_LLLL();
    const Tensor<4, data_t> &get_gauss_codazzi();
    const Tensor<3, data_t> &get_codazzi_mainardi();

    // EM-tensor dependent
    const Tensor<1, data_t> &get_momentum_constraints();
    const Tensor<2, data_t> &get_weyl_electric_part();
    const Tensor<1, data_t> &get_lie_Z();

    // EOM dependent
    const emtensor_t<data_t> &
    get_em_tensor_effective(); // with cosmological constant
    const data_t &get_hamiltonian_constraint();
    const Vars &get_lie_derivatives(); // Lie derivatives of non-gauge variables
    const Tensor<2, data_t> &get_lie_extrinsic_curvature();
    const Tensor<2, data_t> &get_eom_double_normal_projection();

    // Advection dependent
    const Vars &get_rhs_equations();
    // const Tensor<1, data_t> &get_dt_chris_contracted(); // commented, but
    // working!
    // const Tensor<1, data_t> &get_dt_chris_spatial_contracted(); //
    // commented, but working!

    //////// ST ////////
    // ST
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_metric_ST();
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_projector_LU_ST();
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_metric_UU_ST();
    const Tensor<1, data_t, CH_SPACEDIM + 1> &get_normal_U_ST();
    const Tensor<1, data_t, CH_SPACEDIM + 1> &get_normal_L_ST();
    const Tensor<1, data_t, CH_SPACEDIM + 1> &get_shift_ST();
    const Tensor<3, data_t, CH_SPACEDIM + 1> &
    get_levi_civita_spatial_ST();                                   // LLL
    const Tensor<4, data_t, CH_SPACEDIM + 1> &get_levi_civita_ST(); // LLL
    const Tensor<1, data_t, CH_SPACEDIM + 1> &get_Z_L_ST();
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_grad_normal_LL();
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_covd_Z_L_ST();

    // EM-tensor dependent
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_em_tensor_ST();
    const data_t &get_em_tensor_trace_ST();
    const Tensor<4, data_t, CH_SPACEDIM + 1> &get_weyl_tensor_LLLL();
    const data_t &get_weyl_squared();

    // EOM dependent
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_em_tensor_effective_ST();
    const data_t &get_em_tensor_effective_trace_ST();
    const Tensor<4, data_t, CH_SPACEDIM + 1> &get_riemann_LLLL_ST();
    const Tensor<4, data_t, CH_SPACEDIM + 1> &get_riemann_LLLL_ST_v2();
    const Tensor<2, data_t, CH_SPACEDIM + 1> &get_ricci_ST();
    const data_t &get_ricci_scalar_ST();
    const data_t &get_ricci_squared();
    const data_t &get_kretschmann();
    const data_t &get_riemann_squared(); // alternative to Kretschmann

    // Advection dependent
    // const Tensor<3, data_t, CH_SPACEDIM + 1> &get_chris_ST(); // commented,
    // but working!
    // const Tensor<2, data_t, CH_SPACEDIM + 1> &get_d1_Z_L_ST();
    // // commented, but working!

    //////// EXTRA ////////
    ricci_t<data_t> compute_ricci_qDZ(int q); // computes Rij + q * DiZj
    // do we always need the scalar? Should it be separated?

    // extra ideas: extrinsic curvature ST? Acceleration? ST christoffels?

    void clean(); // clean all computations made so far
  protected:
    void set_all_to_null();
    void clean_em_tensor_dependent(); // clean EM-tensor dependent
    void clean_eom_dependent();       // clean EOM dependent variables
    void clean_advection_dependent(); // clean EOM dependent variables

    int m_formulation;
    const CCZ4::params_t *m_ccz4_params;
    double m_16_pi_G_Newton;
    double m_cosmological_constant;

    //////// SET BY USER ////////
    const Vars *m_vars;
    const Diff1Vars *m_d1_vars;
    const Diff2Vars *m_d2_vars;
    const Vars *m_advection;
    const emtensor_t<data_t> *m_em_tensor;

    //////// spatial ////////
    // spatial conformal
    const Tensor<2, data_t> *m_h_UU;
    const chris_t<data_t> *m_chris;
    const Tensor<1, data_t> *m_Z_U_conformal;
    const Tensor<2, data_t> *m_d1_chris_contracted;
    const Tensor<2, data_t> *m_covd_chi_conformal;
    const Tensor<4, data_t> *m_riemann_conformal_LLLL;
    const Tensor<2, data_t> *m_A_LU;
    const Tensor<2, data_t> *m_A_UU;
    const data_t *m_tr_A2;

    // spatial non-conformal
    const Tensor<2, data_t> *m_metric_spatial;
    const Tensor<2, data_t> *m_metric_UU_spatial;
    const Tensor<1, data_t> *m_shift_L;
    const Tensor<2, data_t> *m_extrinsic_curvature;
    const Tensor<3, data_t> *m_chris_spatial;
    const Tensor<1, data_t> *m_Z_U;
    const Tensor<1, data_t> *m_Z;
    const Tensor<2, data_t> *m_covd_Z;
    const ricci_t<data_t> *m_ricci;
    const ricci_t<data_t> *m_ricci_1DZ;
    const ricci_t<data_t> *m_ricci_2DZ;
    const Tensor<3, data_t> *m_d1_extrinsic_curvature;
    const Tensor<3, data_t> *m_covd_extrinsic_curvature;
    const Tensor<3, data_t> *m_levi_civita_spatial;
    const Tensor<3, data_t> *m_levi_civita_spatial_LUU;
    const Tensor<2, data_t> *m_weyl_magnetic_part;
    const Tensor<4, data_t> *m_riemann_spatial_LLLL;
    const Tensor<4, data_t> *m_gauss_codazzi;
    const Tensor<3, data_t> *m_codazzi_mainardi;
    const Tensor<2, data_t> *m_covd_lapse;

    // EM-tensor dependent
    const Tensor<1, data_t> *m_momentum_constraints;
    const Tensor<2, data_t> *m_weyl_electric_part;
    const Tensor<1, data_t> *m_lie_Z;

    // EOM dependent
    const emtensor_t<data_t>
        *m_em_tensor_effective; // with cosmological constant
    const data_t *m_hamiltonian_constraint;
    const Vars *m_lie_derivatives;
    const Tensor<2, data_t> *m_lie_extrinsic_curvature;
    const Tensor<2, data_t> *m_eom_double_normal_projection;

    // Advection dependent
    const Vars *m_rhs_equations;
    // const Tensor<1, data_t> *m_dt_chris_contracted;
    // const Tensor<1, data_t> *m_dt_chris_spatial_contracted;

    //////// ST ////////
    // ST
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_metric_ST;
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_projector_LU_ST;
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_metric_UU_ST;
    const Tensor<1, data_t, CH_SPACEDIM + 1> *m_normal_U_ST;
    const Tensor<1, data_t, CH_SPACEDIM + 1> *m_normal_L_ST;
    const Tensor<1, data_t, CH_SPACEDIM + 1> *m_shift_ST;
    const Tensor<3, data_t, CH_SPACEDIM + 1> *m_levi_civita_spatial_ST;
    const Tensor<4, data_t, CH_SPACEDIM + 1> *m_levi_civita_ST;
    const Tensor<1, data_t, CH_SPACEDIM + 1> *m_Z_L_ST;
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_covd_Z_L_ST;
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_grad_normal_LL;

    // EM-tensor
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_em_tensor_ST;
    const data_t *m_em_tensor_trace_ST;
    const Tensor<4, data_t, CH_SPACEDIM + 1> *m_weyl_tensor_LLLL;
    const data_t *m_weyl_squared;

    // EOM dependent
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_em_tensor_effective_ST;
    const data_t *m_em_tensor_effective_trace_ST;
    const Tensor<4, data_t, CH_SPACEDIM + 1> *m_riemann_LLLL_ST;
    const Tensor<4, data_t, CH_SPACEDIM + 1> *m_riemann_LLLL_ST_v2;
    const Tensor<2, data_t, CH_SPACEDIM + 1> *m_ricci_ST;
    const data_t *m_ricci_scalar_ST;
    const data_t *m_ricci_squared;
    const data_t *m_kretschmann;
    const data_t *m_riemann_squared;

    // Advection dependent
    // const Tensor<3, data_t, CH_SPACEDIM + 1> *m_chris_ST;
    // const Tensor<2, data_t, CH_SPACEDIM + 1> *m_d1_Z_L_ST;

    //////// spatial ////////
    // spatial conformal
    void compute_h_UU();
    void compute_chris();
    void compute_Z_U_conformal();
    void compute_d1_chris_contracted();
    void compute_covd_chi_conformal();
    void compute_riemann_conformal_LLLL();
    void compute_A_LU();
    void compute_A_UU();
    void compute_tr_A2();

    // spatial non-conformal
    void compute_metric_spatial();
    void compute_metric_UU_spatial();
    void compute_shift_L();
    void compute_extrinsic_curvature();
    void compute_em_tensor();
    void compute_chris_spatial();
    void compute_Z_U();
    void compute_Z();
    void compute_covd_Z();
    void compute_ricci();
    void compute_ricci_1DZ();
    void compute_ricci_2DZ();
    void compute_d1_extrinsic_curvature();
    void compute_covd_extrinsic_curvature();
    void compute_levi_civita_spatial();
    void compute_levi_civita_spatial_LUU();
    void compute_weyl_magnetic_part();
    void compute_riemann_spatial_LLLL();
    void compute_gauss_codazzi();
    void compute_codazzi_mainardi();
    void compute_covd_lapse();

    // EM-tensor dependent
    void compute_momentum_constraints();
    void compute_weyl_electric_part();
    void compute_lie_Z();

    // EOM dependent
    void compute_em_tensor_effective();
    void compute_hamiltonian_constraint();
    void compute_lie_derivatives();
    void compute_lie_extrinsic_curvature();
    void compute_eom_double_normal_projection();

    // Advection dependent
    void compute_rhs_equations();
    // void compute_dt_chris_contracted();
    // void compute_dt_chris_spatial_contracted();

    //////// ST ////////
    // ST
    void compute_metric_ST();
    void compute_projector_LU_ST();
    void compute_metric_UU_ST();
    void compute_normal_U_ST();
    void compute_normal_L_ST();
    void compute_shift_ST();
    void compute_levi_civita_spatial_ST();
    void compute_levi_civita_ST();
    void compute_Z_L_ST();
    void compute_covd_Z_L_ST();
    void compute_grad_normal_LL();

    // EM-tensor
    void compute_em_tensor_ST();
    void compute_em_tensor_trace_ST();
    void compute_weyl_tensor_LLLL();
    void compute_weyl_squared();

    // EOM dependent
    void compute_em_tensor_effective_ST();
    void compute_em_tensor_effective_trace_ST();
    void compute_riemann_LLLL_ST();
    void compute_riemann_LLLL_ST_v2();
    void compute_ricci_ST();
    void compute_ricci_scalar_ST();
    void compute_ricci_squared();
    void compute_kretschmann();
    void compute_riemann_squared();

    // Advection dependent
    // void compute_chris_ST();
    // void compute_d1_Z_L_ST();
};

#include "GeometricQuantities.impl.hpp"

#endif /* GEOMETRICQUANTITIES_HPP_ */
