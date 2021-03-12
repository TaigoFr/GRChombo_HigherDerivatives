/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's rootSchectory.
 */

#include "GeometryTestsLevel.hpp"
#include "BoxLoops.hpp"
#include "CCZ4.hpp"
#include "ChiTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "Constraints.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemoval.hpp"

// Initial data
#include "GammaCalculator.hpp"
#include "SchBH_KerrSchild.hpp"

void GeometryTestsLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void GeometryTestsLevel::initialData()
{
    CH_TIME("GeometryTestsLevel::initialData");
    if (m_verbosity)
        pout() << "GeometryTestsLevel::initialData " << m_level << endl;

    // First set everything to zero then calculate initial data  Get the Sch
    // solution in the variables, then calculate the \tilde\Gamma^i numerically
    // as these are non zero and not calculated in the Sch ICs
    BoxLoops::loop(make_compute_pack(SetValue(0.), SchBH(m_p.sch_params, m_dx)),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

void GeometryTestsLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS, disable_simd());
}

void GeometryTestsLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side and set constraints to zero to avoid
    // undefined values
    BoxLoops::loop(CCZ4(m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS, disable_simd());
}

void GeometryTestsLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void GeometryTestsLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}