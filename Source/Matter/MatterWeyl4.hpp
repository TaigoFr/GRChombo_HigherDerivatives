/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERWEYL4_HPP_
#define MATTERWEYL4_HPP_

#include "MatterCCZ4.hpp"
#include "Weyl4.hpp"

//!  Calculates the Weyl4 scalar for spacetimes with matter content
/*!
   This class calculates the Weyl4 scalar real and im parts. It inherits from
   the Weyl4 class and adds in the matter terms as appropriate depending on the
   formulation
*/
template <class matter_t> class MatterWeyl4 : public Weyl4
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    MatterWeyl4(const matter_t &a_matter,
                const std::array<double, CH_SPACEDIM> a_center,
                const double a_dx,
                const int a_formulation = CCZ4Formulation::USE_CCZ4,
                double a_G_Newton = 1.0)
        : Weyl4(a_center, a_dx, a_formulation), m_matter(a_matter),
          m_G_Newton(a_G_Newton)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const matter_t &m_matter; //!< The matter object, e.g. a scalar field
    const double m_G_Newton;  //!< Newton's constant, set to one by default
};

#include "MatterWeyl4.impl.hpp"

#endif /* MATTERWEYL4_HPP_ */
