/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "HigherDerivativesLevel.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "ChiRelaxation.hpp"
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "NewMatterConstraintsWithGauge.hpp"

// For tag cells
#include "ChiTaggingCriterion.hpp"

// Problem specific includes
#include "C2EFT.hpp"
#include "CSystem.hpp"
#include "CdiffDiagnostic.hpp"
#include "ComputeC.hpp"
#include "NCCDiagnostic.hpp"
#include "WeakFieldConditionDiagnostic.hpp"

// BH ID defined here:
#include "SimulationParameters.hpp"

// Things to do at each advance step, after the RK4 is calculated
void HigherDerivativesLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void HigherDerivativesLevel::initialData()
{
    CH_TIME("HigherDerivativesLevel::initialData");
    if (m_verbosity)
        pout() << "HigherDerivativesLevel::initialData " << m_level << endl;

    // Set up the compute class for the SchwKS initial data
    InitialData id(m_p.id_params, m_dx);

    // First set everything to zero then desired initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.0), id), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(
        make_compute_pack(GammaCalculator(m_dx),
                          ComputeC(m_dx, m_p.formulation, m_p.ccz4_params)),
        m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void HigherDerivativesLevel::prePlotLevel()
{
    fillAllGhosts();
    CSystem Csystem(m_p.c_params);
    C2EFT<CSystem> c2eft(Csystem, m_p.hd_params);

    MatterConstraints<C2EFT<CSystem>> constraints(
        c2eft, m_dx, m_p.G_Newton, m_p.formulation, m_p.ccz4_params, c_Ham,
        Interval(c_Mom, c_Mom));
    CdiffDiagnostic Cdiff(m_dx, m_p.formulation, m_p.ccz4_params);
    WeakFieldConditionDiagnostic<CSystem> weakField(
        c2eft, m_dx, m_p.formulation, m_p.ccz4_params);
    NCCDiagnostic<CSystem> ncc(c2eft, m_dx, m_p.formulation, m_p.ccz4_params,
                               m_p.center, m_p.G_Newton, c_NCC_plus,
                               c_NCC_minus, c_NCC_Z4_plus, c_NCC_Z4_minus);

    // these are a lot of diagnostics and it would be more efficient to compute
    // them all in a single class, but it's nicer for now to have them separate.
    // Niceness over performance for now, as these are just Plotfiles
    BoxLoops::loop(make_compute_pack(constraints, Cdiff, weakField, ncc),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void HigherDerivativesLevel::specificEvalRHS(GRLevelData &a_soln,
                                             GRLevelData &a_rhs,
                                             const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    CSystem Csystem(m_p.c_params);
    C2EFT<CSystem> c2eft(Csystem, m_p.hd_params);
    MatterCCZ4<C2EFT<CSystem>> my_ccz4_matter(
        c2eft, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void HigherDerivativesLevel::specificUpdateODE(GRLevelData &a_soln,
                                               const GRLevelData &a_rhs,
                                               Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void HigherDerivativesLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}

void HigherDerivativesLevel::specificPostTimeStep()
{
    CH_TIME("HigherDerivativesLevel::specificPostTimeStep");
#ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        prePlotLevel();
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
}