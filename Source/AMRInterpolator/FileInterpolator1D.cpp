/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "MayDay.H"

// Other includes
#include "FileInterpolator1D.hpp"
#include "GRParmParse.hpp"

// Chombo namespace
#include "UsingNamespace.H"

const int FileInterpolator1D::dim = 1;
const int FileInterpolator1D::order = 4;

FileInterpolator1D::FileInterpolator1D(const std::string &Nlabel,
                                       const std::string &xlabel,
                                       const std::string &ylabel,
                                       double a_y_default_if_missing_file,
                                       bool a_error_if_extrapolating)
    : m_y_default_if_missing_file(a_y_default_if_missing_file),
      m_error_if_extrapolating(a_error_if_extrapolating), interpolator(source)
{
    GRParmParse pp;

    std::vector<double> xs, ys;
    pp.load(Nlabel.c_str(), m_N, 0);
    pp.load(xlabel.c_str(), xs, m_N, {0.});
    pp.load(ylabel.c_str(), ys, m_N, {0.});

    if (m_N >= order) // 4 at least for the Lagrange interpolation of 4th order
    {
        m_dx = xs[1] - xs[0];
        m_x_min = xs[0];
        m_x_max = xs[m_N - 1];

        // quick tests the 'dx' is constant
        double error = 1.e-10;
        CH_assert(fabs(xs[m_N - 1] - xs[m_N - 2] - m_dx) < error);
        CH_assert(fabs(m_x_min + m_dx * (m_N - 1) - m_x_max) < error);

        source.set({m_N}, {m_dx});
        box.set({m_N}, ys);
    }
    else
        MayDay::Warning(
            "No file found or too little points for interpolation.");
}

double FileInterpolator1D::interpolate(double x, int derivative)
{
    if (m_N < order)
        return m_y_default_if_missing_file;

    if (x < m_x_min || x > m_x_max)
    {
        std::string msg =
            "FileInterpolator1D::interpolate - interpolating at " +
            std::to_string(x) + " outside of valid range [" +
            std::to_string(m_x_min) + "," + std::to_string(m_x_max) + "]";
        if (m_error_if_extrapolating)
            MayDay::Abort(msg.c_str());
        else
            MayDay::Warning(msg.c_str());
    }

    interpolator.setup({derivative}, {(x - m_x_min) / m_dx});
    double answer = interpolator.interpData(box);

    // printf("\nx=%.10lf\ty=%.10lf", x, answer);

    return answer;
}
