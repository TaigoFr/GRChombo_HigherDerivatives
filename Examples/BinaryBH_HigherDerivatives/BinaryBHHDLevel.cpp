/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "BinaryBHHDLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// For RHS update
#include "MatterCCZ4RHSWithDiffusion.hpp"

// For constraints calculation
#include "NewMatterConstraintsWithGauge.hpp"

// For tag cells
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"

// BH ID and System defined here:
#include "BinaryBH.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes
#include "C2EFT.hpp"
// #include "DiffusionDiagnostic.hpp"
#include "NCCDiagnostic.hpp"
#include "WeakFieldConditionDiagnostic.hpp"

#ifdef USE_EBSYSTEM
#include "ComputeEB.hpp"
#include "EBdiffDiagnostic.hpp"
#elif USE_CSYSTEM
#include "CDiagnostics.hpp"
#endif

// Things to do at each advance step, after the RK4 is calculated
void BinaryBHHDLevel::specificAdvance()
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
void BinaryBHHDLevel::initialData()
{
    CH_TIME("BinaryBHHDLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBHHDLevel::initialData " << m_level << endl;

    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // First set everything to zero then desired initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();

    // 'GammaCalculator(m_dx)' not needed for binaries (conformally flag initial
#ifdef USE_EBSYSTEM
    ComputeEB compute(m_dx, m_p.formulation, m_p.ccz4_params,
                      Interval(c_E11, c_E33), Interval(c_B11, c_B33));
    if (m_p.system_params.version == 2)
    {
        bool compute_time_derivatives = true;
        ComputeEB compute2(m_dx, m_p.formulation, m_p.ccz4_params,
                           Interval(c_Eaux11, c_Eaux33),
                           Interval(c_Baux11, c_Baux33),
                           compute_time_derivatives);
        BoxLoops::loop(make_compute_pack(compute, compute2), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS);
    }
    else
    {
        BoxLoops::loop(compute, m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    }
#elif USE_CSYSTEM
    CDiagnostics compute(m_dx, m_p.formulation, m_p.ccz4_params, c_C, -1);
    BoxLoops::loop(compute, m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
#endif

#ifdef USE_AHFINDER
    // Diagnostics needed for AHFinder (calculate anyway as AHFinder not yet
    // setup)
    computeDiagnostics();
#endif
}

// Things to do before outputting a plot file
void BinaryBHHDLevel::prePlotLevel()
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
void BinaryBHHDLevel::computeDiagnostics()
{
    fillAllGhosts();
    bool apply_weak_field = false;
    System EBsystem(m_p.system_params);
    C2EFT<System> c2eft(EBsystem, m_p.hd_params, apply_weak_field);

    MatterConstraints<C2EFT<System>> constraints(
        c2eft, m_dx, m_p.G_Newton, m_p.formulation, m_p.ccz4_params, c_Ham,
        Interval(c_Mom, c_Mom));

#ifdef USE_EBSYSTEM
    EBdiffDiagnostic diff(m_dx, m_p.formulation, m_p.ccz4_params);
#elif USE_CSYSTEM
    CDiagnostics diff(m_dx, m_p.formulation, m_p.ccz4_params, c_Cphys,
                      c_C_diff);
#endif

    WeakFieldConditionDiagnostic<System> weakField(c2eft, m_dx, m_p.formulation,
                                                   m_p.ccz4_params);
    NCCDiagnostic<System> ncc(c2eft, m_dx, m_p.formulation, m_p.ccz4_params,
                              m_p.center, m_p.G_Newton, c_NCC_plus, c_NCC_minus,
                              c_NCC_Z4_plus, c_NCC_Z4_minus);
    // DiffusionDiagnostic<C2EFT<System>> diffusion(
    //     c2eft, m_p.ccz4_params, m_p.diffusion_params, m_dx, m_dt, m_p.sigma,
    //     m_p.formulation, m_p.G_Newton, c_diffusion_h11, c_rhs_h11);

#ifdef USE_EBSYSTEM
    if (m_p.system_params.version == 1)
    {
        BoxLoops::loop(ComputeEB(m_dx, m_p.formulation, m_p.ccz4_params,
                                 Interval(c_Eaux11, c_Eaux33),
                                 Interval(c_Baux11, c_Baux33)),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    }
#endif

    BoxLoops::loop(
        make_compute_pack(constraints, diff, weakField, ncc /*, diffusion*/),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void BinaryBHHDLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
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
        // don't include the above BoxLoops in this one because this one
        // excludes ghost cells and the fillAllGhosts will not fill the outer
        // boundaries for sommerfeld BC
        BoxLoops::loop(ComputeEB(m_dx, m_p.formulation, m_p.ccz4_params,
                                 Interval(c_Eaux11, c_Eaux33),
                                 Interval(c_Baux11, c_Baux33)),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
        fillAllGhosts();
    }
#endif

    double spin = 0.;
    m_p.hd_params.update_min_chi(a_time, spin);
    m_p.diffusion_params.update_min_chi(a_time, spin);

    bool apply_weak_field = true;
    System EBsystem(m_p.system_params);
    C2EFT<System> c2eft(EBsystem, m_p.hd_params, apply_weak_field);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHSWithDiffusion<C2EFT<System>, MovingPunctureGauge,
                                   FourthOrderDerivatives>
            my_ccz4_matter(c2eft, m_p.ccz4_params, m_p.diffusion_params, m_dx,
                           m_dt, m_p.sigma, m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHSWithDiffusion<C2EFT<System>, MovingPunctureGauge,
                                   SixthOrderDerivatives>
            my_ccz4_matter(c2eft, m_p.ccz4_params, m_p.diffusion_params, m_dx,
                           m_dt, m_p.sigma, m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void BinaryBHHDLevel::specificUpdateODE(GRLevelData &a_soln,
                                        const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BinaryBHHDLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

void BinaryBHHDLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state)
{
    if (m_p.track_punctures)
    {
        const vector<double> puncture_masses = {m_p.bh1_params.mass,
                                                m_p.bh2_params.mass};
        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void BinaryBHHDLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHHDLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(Weyl4(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();

        bool apply_weak_field = false;
        System EBsystem(m_p.system_params);
        C2EFT<System> c2eft(EBsystem, m_p.hd_params, apply_weak_field);
        MatterConstraints<C2EFT<System>> constraints(
            c2eft, m_dx, m_p.G_Newton, m_p.formulation, m_p.ccz4_params, c_Ham,
            Interval(c_Mom, c_Mom));

        BoxLoops::loop(constraints, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom, c_Mom));
            SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                         m_restart_time, SmallDataIO::APPEND,
                                         first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }

#ifdef USE_AHFINDER
    // if print is on and there are Diagnostics to write, calculate them!
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        computeDiagnostics();
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        if (m_p.AH_set_origins_to_punctures && m_p.track_punctures)
        {
            m_bh_amr.m_ah_finder.set_origins(
                m_bh_amr.m_puncture_tracker.get_puncture_coords());
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif

    if (m_level == 0 && m_time != 0.)
    {
        pout() << "Used estimation of AH at chi = "
               << m_p.hd_params.chi_ignore_threshold << std::endl;
        pout() << "Applied excision at chi ~< " << m_p.hd_params.chi_threshold
               << std::endl;
    }
}
