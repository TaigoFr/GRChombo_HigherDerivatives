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

    c_Ephys11,
    c_Ephys12,
    c_Ephys13,
    c_Ephys22,
    c_Ephys23,
    c_Ephys33,

    c_Bphys11,
    c_Bphys12,
    c_Bphys13,
    c_Bphys22,
    c_Bphys23,
    c_Bphys33,
#elif USE_CSYSTEM
    c_C_diff,
    c_Cphys,
#endif

    c_WeakFieldVar,
    c_WeakFieldCondition,
    c_WeakFieldVar_after_WFC,
    c_WeakField_over_Kretschmann,

    c_NCC_plus,
    c_NCC_minus,
    c_NCC_Z4_plus,
    c_NCC_Z4_minus,

    c_rhs_chi,
    c_diffusion_chi,
    c_det_h,

    c_Madm,
    c_Jadm,

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

    "Ephys11",
    "Ephys12",
    "Ephys13",
    "Ephys22",
    "Ephys23",
    "Ephys33",

    "Bphys11",
    "Bphys12",
    "Bphys13",
    "Bphys22",
    "Bphys23",
    "Bphys33",
#elif USE_CSYSTEM
    "C_diff",
    "Cphys",
#endif

    "WeakFieldVar",
    "WeakFieldCondition",
    "WeakFieldVar_after_WFC",
    "WeakField_over_Kretschmann",

    "NCC_plus",
    "NCC_minus",
    "NCC_Z4_plus",
    "NCC_Z4_minus",

    "rhs_chi",
    "diffusion_chi",
    "det_h",

    "M_adm",
    "J_adm"};
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
