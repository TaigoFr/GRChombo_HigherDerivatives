#ifndef BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS
#define BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS

#include "Schwarzschild_SolvedConstraints.hpp"
#include "TensorAlgebra.hpp"

class BoostedSchwarzschild_SolvedConstraints
    : public Schwarzschild_SolvedConstraints
{
  public:
    struct params_t : Schwarzschild_SolvedConstraints::params_t
    {
        std::array<double, CH_SPACEDIM> boost_velocity;
    };

    BoostedSchwarzschild_SolvedConstraints(params_t a_params, const double a_dx,
                                           const std::string &append = "")
        : Schwarzschild_SolvedConstraints(a_params, a_dx, append),
          m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t>
    adm_metric_t<data_t>
    compute_adm_boosted_vars(const Coordinates<data_t> &boosted_coords) const;

    template <class data_t>
    void
    compute_rest_spacetime_metric(Tensor<2, data_t, CH_SPACETIMEDIM> &g,
                                  Tensor<3, data_t, CH_SPACETIMEDIM> &dg,
                                  const Coordinates<data_t> &rest_coords) const;

    // there will be a m_params copy for the non-boosted params
    // we could access it via this->Schwarzschild_SolvedConstraints::m_params
    // but no need
    const params_t m_params;
};

#include "BoostedSchwarzschild_SolvedConstraints.impl.hpp"

#endif /* BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS*/
