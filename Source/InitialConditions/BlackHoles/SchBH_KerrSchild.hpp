/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCHBH_KERRSCHILD_HPP_
#define SCHBH_KERRSCHILD_HPP_

#include "Cell.hpp"
#include <array>

class SchBH
{
  public:
    struct params_t
    {
        double mass;
        std::array<double, CH_SPACEDIM> center;
    };

    SchBH(params_t a_params, double a_dx) : m_dx(a_dx), m_params(a_params) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    params_t m_params;
};

#include "SchBH_KerrSchild.impl.hpp"

#endif /* SCHBH_KERRSCHILD_HPP_ */
