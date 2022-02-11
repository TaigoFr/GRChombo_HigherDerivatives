#ifndef SCHWARZSCHILD_SOLVEDCONSTRAINTS
#define SCHWARZSCHILD_SOLVEDCONSTRAINTS

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FileInterpolator1D.hpp"
#include "Tensor.hpp"

class Schwarzschild_SolvedConstraints
{
  public:
    //! Stuct for the params of the Schw BH
    struct params_t
    {
        double mass;                            //!<< The mass of the Schw BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the Schw BH
    };

    static double get_bh_radius(params_t a_params)
    {
        return a_params.mass / 2.;
    }

  protected:
    double m_dx;
    const params_t m_params;

    mutable FileInterpolator1D file_psi, file_Krr;

  public:
    Schwarzschild_SolvedConstraints(params_t a_params, const double a_dx,
                                    const std::string &append = "")
        : m_dx(a_dx), m_params(a_params),
          file_psi("Npoints" + append, "rs" + append, "psi" + append,
                   -1.), // negative such that we identify when there is no
                         // file and set the usual chi=(1+M/(2r))^4
                         // (see 'fill_from_files')
          file_Krr("Npoints" + append, "rs" + append, "Krr" + append, 0.)
    {
        file_psi.allow_extrapolation(true);
        file_Krr.allow_extrapolation(true);
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    template <class data_t, template <typename> class vars_t>
    void compute_vars(vars_t<data_t> &vars,
                      const Coordinates<data_t> &coords) const;

    template <class data_t>
    void fill_from_files(data_t &chi, Tensor<2, data_t> &A,
                         const Coordinates<data_t> &coords) const;
};

#include "Schwarzschild_SolvedConstraints.impl.hpp"

#endif /* SCHWARZSCHILD_SOLVEDCONSTRAINTS*/
