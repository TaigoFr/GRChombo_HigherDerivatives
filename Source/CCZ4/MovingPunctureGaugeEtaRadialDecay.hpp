/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGEETARADIALDECAY_HPP_
#define MOVINGPUNCTUREGAUGEETARADIALDECAY_HPP_

#include "DimensionDefinitions.hpp"
#include "GeometricQuantities.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/
class MovingPunctureGaugeEtaRadialDecay
{
  public:
    struct params_t : MovingPunctureGauge::params_t
    {
        double eta_sigmoid_decay = 17.;
        double eta_sigmoid_chi_threshold = 0.92;
        double eta_asymptotic = 0.1;
    };

    const params_t m_params;

    MovingPunctureGaugeEtaRadialDecay(const params_t &a_params)
        : m_params(a_params)
    {
    }

    template <class data_t, template <typename> class vars_rhs_t,
              template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    inline void
    rhs_gauge(vars_rhs_t<data_t> &rhs,
              GeometricQuantities<data_t, vars_t, diff2_vars_t,
                                  MovingPunctureGaugeEtaRadialDecay> &gq) const
    {
        const auto &vars = gq.get_vars();
        const auto &advec = gq.get_advection();

        data_t eta = (m_params.eta_asymptotic - m_params.eta) *
                         sigmoid(vars.chi, m_params.eta_sigmoid_decay,
                                 m_params.eta_sigmoid_chi_threshold) +
                     m_params.eta;

        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                    m_params.lapse_coeff *
                        pow(vars.lapse, m_params.lapse_power) *
                        (vars.K - 2 * vars.Theta);
        FOR(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B[i];
            rhs.B[i] = m_params.shift_advec_coeff * advec.B[i] -
                       m_params.shift_advec_coeff * advec.Gamma[i] +
                       rhs.Gamma[i] - eta * vars.B[i];
        }
    }

    // output is <10^{-k} for x>t(1+k+w) and >1-10^{-k} for x<t(1-k*w)
    template <class data_t>
    inline static data_t sigmoid(data_t x, double decay, double threshold)
    {
        return 1. / (1. + pow((data_t)10., -decay * (x / threshold - 1.)));
    }
};

#endif /* MOVINGPUNCTUREGAUGEETARADIALDECAY_HPP_ */
