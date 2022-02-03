/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CCZ4RHS.hpp" // need 'formulations'; should be included before GeometricQuantities (before the #ifndef)

#ifndef GEOMETRICQUANTITIES_HPP_
#define GEOMETRICQUANTITIES_HPP_

#include "Coordinates.hpp"

// similar to CH_assert in Chombo CH_assert.H
#ifndef NDEBUG
#define assert_with_label(cond, label)                                         \
    if (!(cond))                                                               \
    {                                                                          \
        CH_XD::MayDay::Abort(                                                  \
            (__FILE__ ":" CH_assert_xstr(__LINE__) ": Assertion `" #cond       \
                                                   "' failed. Called from " +  \
             label)                                                            \
                .c_str());                                                     \
    }
#else
#define assert_with_label(cond, label) (void)0
#endif

// auxiliary class as default template for when GeometricQuantities doesn't need
// a gauge
class EmptyGauge
{
  public:
    struct params_t
    {
    };
};

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

//!  Calculates the several spatial and spacetime geometric quantities
/*!
   This class calculates geometric quantities such as the 3-Ricci, the
   4-Weyl, the Kretschmann scalar, etc.. It's made such that calculations
   are minimized, re-using previously calculated quantities whenever
   possible. To be used inside compute functions.
*/
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t, class gauge_t = EmptyGauge>
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

    // a_label is simply used to track where the GeometricQuantities object
    // was created in case of errors
    // Not adding a default label "" to avoid forgetting to add one
    //! Constructor of class GeometricQuantities
    GeometricQuantities(const std::string &a_label);
    GeometricQuantities(const Vars &a_vars, const Diff1Vars &a_d1_vars,
                        const Diff2Vars &a_d2_vars, const std::string &a_label);
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
    void set_advection_and_gauge(const Vars &a_advection,
                                 const gauge_t &a_gauge);
    void set_formulation(int formulation, const CCZ4_params_t<> &a_ccz4_params);
    void set_em_tensor(const emtensor_t<data_t> &a_em_tensor, double G_Newton);
    void set_cosmological_constant(double cosmological_constant);
    void set_coordinates(const Coordinates<data_t> &a_coords);
    /*
    Cheat sheet (L is \Lambda, k = 16 \pi G):
     k T_mn -> -2 L g_mn
     k T    -> -2 (d+1) L
     k rho  -> 2 L
     k S_i  -> 0
     k S_ij -> -2 L h_ij
     k S    -> -2 d L
    */

    void set_all_vars(const Vars &a_vars, const Diff1Vars &a_d1_vars,
                      const Diff2Vars &a_d2_vars);

    //////// SET BY USER ////////
    // vars, d1, d2 and emtensor are assumed to stay in scope
    int get_formulation() const;
    const CCZ4_params_t<> &get_formulation_params() const;
    double get_cosmological_constant() const;
    const Vars &get_vars() const;
    const Diff1Vars &get_d1_vars() const;
    const Diff2Vars &get_d2_vars() const;
    const Vars &get_advection() const;
    const gauge_t &get_gauge() const;
    const emtensor_t<data_t> &get_em_tensor() const;
    const Coordinates<data_t> &get_coordinates() const;

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
    const data_t &get_div_shift();
    const Tensor<1, data_t> &get_Gamma_L();

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
    const Tensor<2, data_t> &get_weyl_magnetic_part();
    const Tensor<4, data_t> &get_riemann_spatial_LLLL();
    const Tensor<4, data_t> &get_gauss_codazzi();
    const Tensor<3, data_t> &get_codazzi_mainardi();
    const Tensor<1, data_t> &get_Gamma_spatial();
    const Tensor<1, data_t> &get_Gamma_L_spatial();
    const Tensor<1, data_t> &get_acceleration_spatial();

    // EM-tensor dependent
    const Tensor<1, data_t> &get_momentum_constraints();
    const Tensor<1, data_t> &get_lie_Z();

    // EOM / Formulation dependent
    const ricci_t<data_t> &get_ricci_qDZ(int q);
    const ricci_t<data_t> &get_ricci(); // original formula
    const ricci_t<data_t> &get_ricci_1DZ();
    // uses evolved variables, not calculated Gammas
    const ricci_t<data_t> &get_ricci_2DZ();
    const Tensor<2, data_t> &get_weyl_electric_part();
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
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_metric_ST();
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_projector_LU_ST();
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_metric_UU_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_normal_U_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_normal_L_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_shift_ST();
    const Tensor<3, data_t, CH_SPACETIMEDIM> &
    get_levi_civita_spatial_ST();                                   // LLL
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_levi_civita_ST(); // LLL
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_Z_L_ST();
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_grad_normal_LL();
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_covd_Z_L_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_acceleration_ST();

    // EM-tensor dependent
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_em_tensor_ST();
    const data_t &get_em_tensor_trace_ST();
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_weyl_tensor_LLLL();
    const data_t &get_weyl_squared();

    // EOM / Formulation dependent
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_riemann_LLLL_ST();
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_riemann_LLLL_ST_v2();
    const Tensor<2, data_t, CH_SPACETIMEDIM> &get_ricci_ST();
    const data_t &get_ricci_scalar_ST();
    const data_t &get_ricci_squared();
    const data_t &get_kretschmann();
    const data_t &get_riemann_squared(); // alternative to Kretschmann
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_riemann_LLLU_ST();
    const Tensor<4, data_t, CH_SPACETIMEDIM> &get_riemann_LULU_ST();

    // Advection dependent
    const Tensor<3, data_t, CH_SPACETIMEDIM> &get_chris_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_Gamma_ST();
    const Tensor<1, data_t, CH_SPACETIMEDIM> &get_Gamma_L_ST();
    // const Tensor<2, data_t, CH_SPACETIMEDIM> &get_d1_Z_L_ST();
    // // commented, but working!

    //////// EXTRA ////////
    ricci_t<data_t> compute_ricci_qDZ(int q); // computes Rij + q * DiZj
    void compute_rhs_equations(Vars &);
    void compute_rhs_equations_no_gauge(Vars &);
    Tensor<4, data_t, CH_SPACETIMEDIM>
    compute_weyl_tensor_LLLL(const Tensor<2, data_t> &Eij,
                             const Tensor<2, data_t> &Bij);

    //////// NON-STANDARD ////////
    // ignores matter, assumes E is variable
    Tensor<2, data_t>
    compute_LieD_weyl_electric_part(const Tensor<2, Tensor<1, data_t>> &d1Eij,
                                    const Tensor<2, Tensor<1, data_t>> &d1Bij,
                                    const Tensor<2, data_t> &Eij,
                                    const Tensor<2, data_t> &Bij);
    Tensor<2, data_t>
    compute_LieD_weyl_magnetic_part(const Tensor<2, Tensor<1, data_t>> &d1Eij,
                                    const Tensor<2, Tensor<1, data_t>> &d1Bij,
                                    const Tensor<2, data_t> &Eij,
                                    const Tensor<2, data_t> &Bij);
    // require advection:
    Tensor<2, data_t> compute_dt_weyl_electric_part(
        const Tensor<2, Tensor<1, data_t>> &d1Eij,
        const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
        const Tensor<2, data_t> &Bij, const Tensor<2, data_t> &advec_Eij,
        const Tensor<2, data_t> &advec_Bij);
    Tensor<2, data_t> compute_dt_weyl_magnetic_part(
        const Tensor<2, Tensor<1, data_t>> &d1Eij,
        const Tensor<2, Tensor<1, data_t>> &d1Bij, const Tensor<2, data_t> &Eij,
        const Tensor<2, data_t> &Bij, const Tensor<2, data_t> &advec_Eij,
        const Tensor<2, data_t> &advec_Bij);

    // extra ideas: extrinsic curvature ST? Acceleration?

    void clean(); // clean all computations made so far
  protected:
    void set_all_to_null();
    void clean_em_tensor_dependent();           // clean EM-tensor dependent
    void clean_eom_dependent();                 // clean EOM dependent variables
    void clean_advection_and_gauge_dependent(); // clean EOM dependent variables

    std::string m_label;
    int m_formulation;
    const CCZ4_params_t<> *m_ccz4_params;
    double m_16_pi_G_Newton;
    double m_cosmological_constant;

    //////// SET BY USER ////////
    const Vars *m_vars;
    const Diff1Vars *m_d1_vars;
    const Diff2Vars *m_d2_vars;
    const Vars *m_advection;
    const gauge_t *m_gauge;
    const emtensor_t<data_t> *m_em_tensor;
    const Coordinates<data_t> *m_coords;

    //////// spatial ////////
    // spatial conformal
    Tensor<2, data_t> *m_h_UU;
    chris_t<data_t> *m_chris;
    Tensor<1, data_t> *m_Z_U_conformal;
    Tensor<2, data_t> *m_d1_chris_contracted;
    Tensor<2, data_t> *m_covd_chi_conformal;
    Tensor<4, data_t> *m_riemann_conformal_LLLL;
    Tensor<2, data_t> *m_A_LU;
    Tensor<2, data_t> *m_A_UU;
    data_t *m_tr_A2;
    data_t *m_div_shift;
    Tensor<1, data_t> *m_Gamma_L;

    // spatial non-conformal
    Tensor<2, data_t> *m_metric_spatial;
    Tensor<2, data_t> *m_metric_UU_spatial;
    Tensor<1, data_t> *m_shift_L;
    Tensor<2, data_t> *m_extrinsic_curvature;
    Tensor<3, data_t> *m_chris_spatial;
    Tensor<1, data_t> *m_Z_U;
    Tensor<1, data_t> *m_Z;
    Tensor<2, data_t> *m_covd_Z;
    Tensor<3, data_t> *m_d1_extrinsic_curvature;
    Tensor<3, data_t> *m_covd_extrinsic_curvature;
    Tensor<3, data_t> *m_levi_civita_spatial;
    Tensor<3, data_t> *m_levi_civita_spatial_LUU;
    Tensor<2, data_t> *m_weyl_magnetic_part;
    Tensor<4, data_t> *m_riemann_spatial_LLLL;
    Tensor<4, data_t> *m_gauss_codazzi;
    Tensor<3, data_t> *m_codazzi_mainardi;
    Tensor<2, data_t> *m_covd_lapse;
    Tensor<1, data_t> *m_Gamma_spatial;
    Tensor<1, data_t> *m_Gamma_L_spatial;
    Tensor<1, data_t> *m_acceleration_spatial;

    // EM-tensor dependent
    Tensor<1, data_t> *m_momentum_constraints;
    Tensor<2, data_t> *m_weyl_electric_part;
    Tensor<1, data_t> *m_lie_Z;

    // EOM dependent
    ricci_t<data_t> *m_ricci;
    ricci_t<data_t> *m_ricci_1DZ;
    ricci_t<data_t> *m_ricci_2DZ;
    data_t *m_hamiltonian_constraint;
    Vars *m_lie_derivatives;
    Tensor<2, data_t> *m_lie_extrinsic_curvature;
    Tensor<2, data_t> *m_eom_double_normal_projection;

    // Advection dependent
    Vars *m_rhs_equations;
    // Tensor<1, data_t> *m_dt_chris_contracted;
    // Tensor<1, data_t> *m_dt_chris_spatial_contracted;

    //////// ST ////////
    // ST
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_metric_ST;
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_projector_LU_ST;
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_metric_UU_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_normal_U_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_normal_L_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_shift_ST;
    Tensor<3, data_t, CH_SPACETIMEDIM> *m_levi_civita_spatial_ST;
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_levi_civita_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_Z_L_ST;
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_grad_normal_LL;
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_covd_Z_L_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_acceleration_ST;

    // EM-tensor
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_em_tensor_ST;
    data_t *m_em_tensor_trace_ST;
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_weyl_tensor_LLLL;
    data_t *m_weyl_squared;

    // EOM dependent
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_riemann_LLLL_ST;
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_riemann_LLLL_ST_v2;
    Tensor<2, data_t, CH_SPACETIMEDIM> *m_ricci_ST;
    data_t *m_ricci_scalar_ST;
    data_t *m_ricci_squared;
    data_t *m_kretschmann;
    data_t *m_riemann_squared;
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_riemann_LLLU_ST;
    Tensor<4, data_t, CH_SPACETIMEDIM> *m_riemann_LULU_ST;

    // Advection dependent
    Tensor<3, data_t, CH_SPACETIMEDIM> *m_chris_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_Gamma_ST;
    Tensor<1, data_t, CH_SPACETIMEDIM> *m_Gamma_L_ST;
    // Tensor<2, data_t, CH_SPACETIMEDIM> *m_d1_Z_L_ST;

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
    void compute_div_shift();
    void compute_Gamma_L();

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
    void compute_d1_extrinsic_curvature();
    void compute_covd_extrinsic_curvature();
    void compute_levi_civita_spatial();
    void compute_levi_civita_spatial_LUU();
    void compute_weyl_magnetic_part();
    void compute_riemann_spatial_LLLL();
    void compute_gauss_codazzi();
    void compute_codazzi_mainardi();
    void compute_covd_lapse();
    void compute_Gamma_spatial();
    void compute_Gamma_L_spatial();
    void compute_acceleration_spatial();

    // EM-tensor dependent
    void compute_momentum_constraints();
    void compute_weyl_electric_part();
    void compute_lie_Z();

    // EOM dependent
    void compute_ricci();
    void compute_ricci_1DZ();
    void compute_ricci_2DZ();
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
    void compute_grad_normal_LL();
    void compute_covd_Z_L_ST();
    void compute_acceleration_ST();

    // EM-tensor
    void compute_em_tensor_ST();
    void compute_em_tensor_trace_ST();
    void compute_weyl_tensor_LLLL();
    void compute_weyl_squared();

    // EOM dependent
    void compute_riemann_LLLL_ST();
    void compute_riemann_LLLL_ST_v2();
    void compute_ricci_ST();
    void compute_ricci_scalar_ST();
    void compute_ricci_squared();
    void compute_kretschmann();
    void compute_riemann_squared();
    void compute_riemann_LLLU_ST();
    void compute_riemann_LULU_ST();

    // Advection dependent
    void compute_chris_ST();
    void compute_Gamma_ST();
    void compute_Gamma_L_ST();
    // void compute_d1_Z_L_ST();
};

#include "GeometricQuantities.impl.hpp"

#endif /* GEOMETRICQUANTITIES_HPP_ */
