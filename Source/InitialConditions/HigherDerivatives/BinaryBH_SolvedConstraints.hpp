/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBH_SOLVEDCONSTRAINTS_HPP_
#define BINARYBH_SOLVEDCONSTRAINTS_HPP_

#include "BoostedSchwarzschild_SolvedConstraints.hpp"
#include "Cell.hpp"

class BinaryBH_SolvedConstraints
{
  protected:
    double m_dx;
    BoostedSchwarzschild_SolvedConstraints m_bh1;
    BoostedSchwarzschild_SolvedConstraints m_bh2;

  public:
    using params_t = typename BoostedSchwarzschild_SolvedConstraints::params_t;

    BinaryBH_SolvedConstraints(params_t a_bh1_params, params_t a_bh2_params,
                               double a_dx);

    template <class data_t> void compute(Cell<data_t> current_cell) const;
};

#include "BinaryBH_SolvedConstraints.impl.hpp"

#endif /* BINARYBH_SOLVEDCONSTRAINTS_HPP_ */
