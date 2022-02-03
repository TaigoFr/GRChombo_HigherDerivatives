/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBH_SOLVEDCONSTRAINTS_HPP_
#define BINARYBH_SOLVEDCONSTRAINTS_HPP_

#include "Cell.hpp"
#include "Schwarzschild_SolvedConstraints.hpp"

class BinaryBH_SolvedConstraints
{
  protected:
    double m_dx;
    Schwarzschild_SolvedConstraints m_bh1;
    Schwarzschild_SolvedConstraints m_bh2;

  public:
    using params_t = typename Schwarzschild_SolvedConstraints::params_t;

    BinaryBH_SolvedConstraints(params_t a_bh1_params, params_t a_bh2_params,
                               double a_dx);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t, template <typename> class vars_t>
    void compute_vars_superposition(vars_t<data_t> &vars,
                                    const vars_t<data_t> &vars1,
                                    const vars_t<data_t> &vars2) const;
};

#include "BinaryBH_SolvedConstraints.impl.hpp"

#endif /* BINARYBH_SOLVEDCONSTRAINTS_HPP_ */
