/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,
    c_Mom,

#ifdef USE_EBSYSTEM
    c_E_diff,
    c_B_diff,
#elif USE_CSYSTEM
    c_C_diff,
    c_Cphys,
#endif

    c_WeakFieldVar,
    c_WeakFieldCondition,
    c_WeakFieldVar_after_WFC,

    c_NCC_plus,
    c_NCC_minus,
    c_NCC_Z4_plus,
    c_NCC_Z4_minus,

    c_rhs_chi,
    c_diffusion_chi,
    c_det_h,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",
    "Mom",

#ifdef USE_EBSYSTEM
    "E_diff",
    "B_diff",
#elif USE_CSYSTEM
    "C_diff",
    "Cphys",
#endif

    "WeakFieldVar",
    "WeakFieldCondition",
    "WeakFieldVar_after_WFC",

    "NCC_plus",
    "NCC_minus",
    "NCC_Z4_plus",
    "NCC_Z4_minus",

    "rhs_chi",
    "diffusion_chi",
    "det_h"};
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
