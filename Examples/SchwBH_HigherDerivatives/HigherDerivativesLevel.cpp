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
#include "MatterCCZ4RHSWithDiffusion.hpp"

// For constraints calculation
#include "NewMatterConstraintsWithGauge.hpp"

// For tag cells
#include "ChiExtractionTaggingCriterion.hpp"

// BH ID and System defined here:
#include "SimulationParameters.hpp"

// Problem specific includes
#include "C2EFT.hpp"
#include "DiffusionDiagnostic.hpp"
#include "NCCDiagnostic.hpp"
#include "WeakFieldConditionDiagnostic.hpp"

#ifdef USE_EBSYSTEM
#include "ComputeEB.hpp"
#include "EBdiffDiagnostic.hpp"
#elif USE_CSYSTEM
#include "CDiagnostics.hpp"
#endif

#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"

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
    BoxLoops::loop(make_compute_pack(SetValue(0.), id), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();

#ifdef USE_EBSYSTEM
    ComputeEB compute(m_dx, m_p.formulation, m_p.ccz4_params,
                      Interval(c_E11, c_E33), Interval(c_B11, c_B33),
                      m_p.system_params.use_last_index_raised);
    BoxLoops::loop(make_compute_pack(GammaCalculator(m_dx), compute),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

    if (m_p.system_params.version == 2)
    {
        bool compute_time_derivatives = true;
        ComputeEB compute2(
            m_dx, m_p.formulation, m_p.ccz4_params,
            Interval(c_Eaux11, c_Eaux33), Interval(c_Baux11, c_Baux33),
            m_p.system_params.use_last_index_raised, compute_time_derivatives);

        BoxLoops::loop(make_compute_pack(compute2), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS);
    }
#elif USE_CSYSTEM
    CDiagnostics compute(m_dx, m_p.formulation, m_p.ccz4_params, c_C, -1);
    BoxLoops::loop(make_compute_pack(GammaCalculator(m_dx), compute),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
#endif

#ifdef USE_AHFINDER
    // Diagnostics needed for AHFinder (calculate anyway as AHFinder not yet
    // setup)
    computeDiagnostics();
#endif
}

// Things to do before outputting a plot file
void HigherDerivativesLevel::prePlotLevel()
{
#ifdef USE_AHFINDER
    // already calculated in 'specificPostTimeStep' or in 'initialData'
    if ((m_time == 0. && m_p.AH_activate) ||
        m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        return;
#endif

    computeDiagnostics();
}

// Things to do before outputting a plot file
void HigherDerivativesLevel::computeDiagnostics()
{
    fillAllGhosts();
    bool apply_weak_field = false;
    System EBsystem(m_p.system_params);
    C2EFT<System> c2eft(EBsystem, m_p.hd_params, apply_weak_field);

    MatterConstraints<C2EFT<System>> constraints(
        c2eft, m_dx, m_p.G_Newton, m_p.formulation, m_p.ccz4_params, c_Ham,
        Interval(c_Mom, c_Mom));

#ifdef USE_EBSYSTEM
    EBdiffDiagnostic diff(m_dx, m_p.formulation, m_p.ccz4_params,
                          m_p.system_params.use_last_index_raised);
#elif USE_CSYSTEM
    CDiagnostics diff(m_dx, m_p.formulation, m_p.ccz4_params, c_Cphys,
                      c_C_diff);
#endif

    WeakFieldConditionDiagnostic<System> weakField(c2eft, m_dx, m_p.formulation,
                                                   m_p.ccz4_params);
    NCCDiagnostic<System> ncc(c2eft, m_dx, m_p.formulation, m_p.ccz4_params,
                              m_p.center, m_p.G_Newton, c_NCC_plus, c_NCC_minus,
                              c_NCC_Z4_plus, c_NCC_Z4_minus);
    DiffusionDiagnostic<C2EFT<System>> diffusion(
        c2eft, m_p.ccz4_params, m_p.diffusion_params, m_dx, m_dt, m_p.sigma,
        m_p.formulation, m_p.G_Newton, c_diffusion_chi, c_rhs_chi, c_det_h);

#ifdef USE_EBSYSTEM
    if (m_p.system_params.version == 1)
    {
        // need to first store Eaux and Baux on the grid with the physical E & B
        BoxLoops::loop(ComputeEB(m_dx, m_p.formulation, m_p.ccz4_params,
                                 Interval(c_Eaux11, c_Eaux33),
                                 Interval(c_Baux11, c_Baux33),
                                 m_p.system_params.use_last_index_raised),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    }
#endif

    BoxLoops::loop(
        make_compute_pack(constraints, diff, weakField, ncc, diffusion),
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

#ifdef USE_EBSYSTEM
    if (m_p.system_params.version == 1)
    {
        // need to first store Eaux and Baux on the grid with the physical E & B
        // don't include the above BoxLoops in this one because this one
        // excludes ghost cells and the fillAllGhosts will not fill the outer
        // boundaries for sommerfeld BC
        BoxLoops::loop(ComputeEB(m_dx, m_p.formulation, m_p.ccz4_params,
                                 Interval(c_Eaux11, c_Eaux33),
                                 Interval(c_Baux11, c_Baux33),
                                 m_p.system_params.use_last_index_raised),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
        fillAllGhosts();
    }
#endif

    // m_p.hd_params.update_min_chi(a_time, m_p.id_params.spin);
    double spin = 0.;
    m_p.hd_params.update_min_chi(a_time, spin);
    m_p.diffusion_params.update_min_chi(a_time, spin);

    bool apply_weak_field = true;
    System EBsystem(m_p.system_params);
    C2EFT<System> c2eft(EBsystem, m_p.hd_params, apply_weak_field);
    MatterCCZ4RHSWithDiffusion<C2EFT<System>> my_ccz4_matter(
        c2eft, m_p.ccz4_params, m_p.diffusion_params, m_dx, m_dt, m_p.sigma,
        m_p.formulation, m_p.G_Newton);
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
    BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                 m_p.extraction_params,
                                                 m_p.activate_extraction),
                   current_state, tagging_criterion);
}

void HigherDerivativesLevel::specificPostTimeStep()
{
    CH_TIME("HigherDerivativesLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    // Do the extraction on the min extraction level
    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_adm = at_level_timestep_multiple(min_level);
        if (calculate_adm)
        {
            // Populate the ADM Mass and Spin values on the grid
            fillAllGhosts();
            BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx,
                                         c_Madm, c_Jadm),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            if (m_level == min_level)
            {
                CH_TIME("ADMExtraction");
                // Now refresh the interpolator and do the interpolation
                m_gr_amr.m_interpolator->refresh();
                ADMQuantitiesExtraction my_extraction(
                    m_p.extraction_params, m_dt, m_time, first_step,
                    m_restart_time, c_Madm, c_Jadm);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

#ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        computeDiagnostics();
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif

    if (m_level == 0)
    {
        pout() << "Used estimation of AH at chi = "
               << m_p.hd_params.chi_ignore_threshold << std::endl;
        pout() << "Applied excision at chi ~< " << m_p.hd_params.chi_threshold
               << std::endl;
    }
}
