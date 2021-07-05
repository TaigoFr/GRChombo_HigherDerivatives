/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

// assign an enum to each variable
enum
{
// Note that it is important that the first enum value is set to 1 more than
// the last CCZ4 var enum
#ifdef USE_EBSYSTEM
    c_E11 = NUM_CCZ4_VARS,
    c_E12,
    c_E13,
    c_E22,
    c_E23,
    c_E33,

    c_B11,
    c_B12,
    c_B13,
    c_B22,
    c_B23,
    c_B33,

    c_Eaux11,
    c_Eaux12,
    c_Eaux13,
    c_Eaux22,
    c_Eaux23,
    c_Eaux33,

    c_Baux11,
    c_Baux12,
    c_Baux13,
    c_Baux22,
    c_Baux23,
    c_Baux33,

#elif USE_CSYSTEM
    c_C = NUM_CCZ4_VARS,
    c_dCdt,
#else
#error "Please define either USE_CSYSTEM or USE_EBSYSTEM"
#endif

    NUM_VARS
};

namespace UserVariables
{
#ifdef USE_EBSYSTEM
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {
        "E11",    "E12",    "E13",    "E22",    "E23",    "E33",

        "B11",    "B12",    "B13",    "B22",    "B23",    "B33",

        "Eaux11", "Eaux12", "Eaux13", "Eaux22", "Eaux23", "Eaux33",

        "Baux11", "Baux12", "Baux13", "Baux22", "Baux23", "Baux33"};

#elif USE_CSYSTEM
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"C", "dCdt"};
#endif

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
