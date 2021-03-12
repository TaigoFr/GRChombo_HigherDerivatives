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

    c_C_diff,

    c_WeakFieldVar,
    c_WeakFieldCondition,
    c_WeakFieldVar_after_WFC,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",
    "Mom",
    "C_diff",
    "WeakFieldVar",
    "WeakFieldCondition",
    "WeakFieldVar_after_WFC"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
