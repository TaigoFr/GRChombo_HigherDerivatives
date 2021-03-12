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

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"E11",
                           "E12",
                           "E13",
                           "E22",
                           "E23",
                           "E33",

                           "B11",
                           "B12",
                           "B13",
                           "B22",
                           "B23",
                           "B33"

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
                           "Bphys33"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
