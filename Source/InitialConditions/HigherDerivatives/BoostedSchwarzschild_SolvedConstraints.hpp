#ifndef BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS
#define BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS

#include "Schwarzschild_SolvedConstraints.hpp"
#include "TensorAlgebra.hpp"

class BoostedSchwarzschild_SolvedConstraints
    : public Schwarzschild_SolvedConstraints
{
  public:
    BoostedSchwarzschild_SolvedConstraints(
        params_t a_params, const double a_dx, const std::string &append = "",
        const std::array<double, CH_SPACEDIM> &a_boost_velocity = {0.})
        : Schwarzschild_SolvedConstraints(a_params, a_dx, append),
          m_boost_velocity(a_boost_velocity)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t>
    TensorAlgebra::adm_metric_t<data_t>
    compute_adm_boosted_vars(const Coordinates<data_t> &boosted_coords) const;

    template <class data_t>
    void
    compute_rest_spacetime_metric(Tensor<2, data_t, CH_SPACETIMEDIM> &g,
                                  Tensor<3, data_t, CH_SPACETIMEDIM> &dg,
                                  const Coordinates<data_t> &rest_coords) const;

    template <class data_t, template <typename> class vars_t>
    void compute_conformal_variables(
        vars_t<data_t> &vars,
        TensorAlgebra::adm_metric_t<data_t> adm_vars) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_boost_velocity;
};

#include "BoostedSchwarzschild_SolvedConstraints.impl.hpp"

#endif /* BOOSTEDSCHWARZSCHILD_SOLVEDCONSTRAINTS*/
