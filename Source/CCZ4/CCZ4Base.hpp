/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4BASE_HPP_
#define CCZ4BASE_HPP_

/// Base parameter struct for CCZ4
/** This struct collects the gauge independent CCZ4 parameters i.e. the damping
 * ones
 */
struct CCZ4_base_params_t
{
    double kappa1;    //!< Damping parameter kappa1 as in arXiv:1106.2254
    double kappa2;    //!< Damping parameter kappa2 as in arXiv:1106.2254
    double kappa3;    //!< Damping parameter kappa3 as in arXiv:1106.2254
    bool covariantZ4; //!< if true, replace kappa1->kappa1/lapse as in
                      //!<  arXiv:1307.7391 eq. 27
};

enum CCZ4Formulation
{
    USE_CCZ4,
    USE_BSSN
};

#endif /* CCZ4BASE_HPP_ */
