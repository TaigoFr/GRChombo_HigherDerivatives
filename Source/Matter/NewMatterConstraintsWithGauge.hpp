/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NEWMATTERCONSTRAINTSWITHGAUGE_HPP_
#define NEWMATTERCONSTRAINTSWITHGAUGE_HPP_

#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "NewConstraints.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

// Same as NewMatterConstraints, but needs gauge variables, advection,
// formulation and emtensor set

//!  Calculates the Hamiltonian and Momentum constraints with matter fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the matter terms. It adds in the matter terms for a given
   matter class matter_t, which must provide it with the Energy Momentum Tensor.
   For an example of a matter_t class see ScalarField. \sa Constraints(),
   ScalarField()
*/
template <class matter_t> class MatterConstraints : public Constraints
{
  public:
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    template <class data_t>
    struct MatterMetricVarsWithGauge : public CCZ4::Vars<data_t>,
                                       public MatterVars<data_t>
    {
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4::Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    template <class data_t>
    struct MatterDiff2MetricVarsWithGauge : public CCZ4::Diff2Vars<data_t>,
                                            public MatterDiff2Vars<data_t>
    {
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4::Diff2Vars<data_t>::enum_mapping(mapping_function);
            MatterDiff2Vars<data_t>::enum_mapping(mapping_function);
        }
    };

    // added formulation, needed for higher derivative EM-tensors through
    // GeometricQuantities class
    //! Constructor of class MatterConstraints
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    MatterConstraints(const matter_t a_matter, double dx, double G_Newton,
                      int formulation, CCZ4::params_t a_params, int a_c_Ham,
                      const Interval &a_c_Moms, int a_c_Ham_abs_terms = -1,
                      const Interval &a_c_Moms_abs_terms = Interval());

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    matter_t my_matter; //!< The matter object, e.g. a scalar field
    double m_G_Newton;  //!< Newton's constant, set to one by default.
    int m_formulation;
    CCZ4::params_t m_params;
};

#include "NewMatterConstraintsWithGauge.impl.hpp"

#endif /* NEWMATTERCONSTRAINTSWITHGAUGE_HPP_ */
