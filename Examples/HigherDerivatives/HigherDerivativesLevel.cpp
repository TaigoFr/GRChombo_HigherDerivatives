/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "HigherDerivativesLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "NewMatterConstraintsWithGauge.hpp"

// For tag cells
#include "ChiTaggingCriterion.hpp"

// Problem specific includes
#include "C2EFT.hpp"
#include "ChiRelaxation.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"

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
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
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

    // fillAllGhosts();
    // BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
    // EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void HigherDerivativesLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(MatterConstraints<C2EFT>(C2EFT(m_p.hd_params), m_dx,
                                            m_p.G_Newton, m_p.formulation,
                                            m_p.ccz4_params, c_Ham,
                                            Interval(c_Mom, c_Mom)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void HigherDerivativesLevel::specificEvalRHS(GRLevelData &a_soln,
                                             GRLevelData &a_rhs,
                                             const double a_time)
{

    // Relaxation function for chi - this will eventually be done separately
    // with hdf5 as input
    if (m_time < m_p.relaxtime)
    {
        // Calculate chi relaxation right hand side
        // Note this assumes conformal chi and Mom constraint trivially
        // satisfied.
        // No evolution in other variables, which are assumed to
        // satisfy constraints per initial conditions
        ChiRelaxation<C2EFT> relaxation(C2EFT(m_p.hd_params), m_dx,
                                        m_p.relaxspeed, m_p.G_Newton);
        SetValue set_other_values_zero(0.0, Interval(c_h11, NUM_VARS - 1));
        auto compute_pack1 =
            make_compute_pack(relaxation, set_other_values_zero);
        BoxLoops::loop(compute_pack1, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else
    {

        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(),
                              PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
            a_soln, a_soln, INCLUDE_GHOST_CELLS);

        // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
        MatterCCZ4<C2EFT> my_ccz4_matter(C2EFT(m_p.hd_params), m_p.ccz4_params,
                                         m_dx, m_p.sigma, m_p.formulation,
                                         m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
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
