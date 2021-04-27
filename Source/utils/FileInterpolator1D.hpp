/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FILEINTERPOLATOR1D_HPP_
#define FILEINTERPOLATOR1D_HPP_

#include <string>
#include <vector>

#include "Lagrange.hpp"
#include "SimpleArrayBox.hpp"
#include "SimpleInterpSource.hpp"

//! Calculate 'x's, 'y's elsewhere, import them to here
//! and this class interpolates with 4th order Lagrange
class FileInterpolator1D
{
  public:
    FileInterpolator1D(const std::string &Nlabel, const std::string &xlabel,
                       const std::string &ylabel,
                       double a_y_default_if_missing_file = 0.,
                       bool a_error_if_extrapolating = false);

    double interpolate(double x, int derivative = 0);

  private:
    static const int dim;
    static const int order;

    int m_N;
    double m_dx, m_x_min, m_x_max; // assume dx = constant
    double m_y_default_if_missing_file;
    bool m_error_if_extrapolating;

    SimpleInterpSource<1> source;
    SimpleArrayBox<1> box;
    Lagrange<4, 1> interpolator;
};

#endif /*  FILEINTERPOLATOR1D_HPP_ */
