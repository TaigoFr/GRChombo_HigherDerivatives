/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LORENTZBOOSTS_HPP_
#define LORENTZBOOSTS_HPP_

#include "AlwaysInline.hpp"
#include "ArrayTools.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

namespace LorentzBoosts
{

// returns the Lorentz factor often called gamma
ALWAYS_INLINE double lorentz_factor(double speed_sq)
{
    return 1.0 / sqrt(1.0 - speed_sq);
}

// returns the Lorentz boost matrix (Jacobian) for provided velocity
// Note the reversed sign convention
template <class data_t>
static Tensor<2, data_t, CH_SPACETIMEDIM>
matrix(const std::array<double, GR_SPACEDIM> &a_velocity)
{
    double speed_sq = ArrayTools::norm_squared(a_velocity);
    double lorentz_fac = lorentz_factor(speed_sq);

    Tensor<2, data_t, CH_SPACETIMEDIM> out;
    // if boost is zero, this term is irrelevant.
    double speed_sq_inverse = 1.0 / std::max(speed_sq, 1.0e-8);

    // spatial components
    FOR(i, j)
    {
        out[i + 1][j + 1] =
            TensorAlgebra::delta(i, j) + (lorentz_fac - 1.0) * a_velocity[i] *
                                             a_velocity[j] * speed_sq_inverse;
    }

    // mixed components
    FOR(i)
    {
        // reverse sign convention here
        out[i + 1][0] = lorentz_fac * a_velocity[i];
        out[0][i + 1] = out[i + 1][0];
    }

    // time component
    out[0][0] = lorentz_fac;

    return out;
}

// returns the inverse Lorentz boost matrix (Jacobian) for provided velocity
// Note the reversed sign convention
template <class data_t>
static Tensor<2, data_t, CH_SPACETIMEDIM>
inverse_matrix(const std::array<double, GR_SPACEDIM> &a_velocity)
{
    std::array<double, GR_SPACEDIM> minus_velocity;
    FOR(i) { minus_velocity[i] = -a_velocity[i]; }
    return matrix<data_t>(minus_velocity);
}

// unboost spatial coordinates to convert boosted spatial coordinates on grid
// to rest frame spatial coordinates
template <class data_t>
static Coordinates<data_t>
unboost_coords(const Coordinates<data_t> &a_boosted_coords,
               const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto inverse_boost_matrix = inverse_matrix<data_t>(a_velocity);

    Tensor<1, data_t> boosted_coords_tensor;
    boosted_coords_tensor[0] = a_boosted_coords.x;
    boosted_coords_tensor[1] = a_boosted_coords.y;
    boosted_coords_tensor[2] = a_boosted_coords.z;

    Tensor<1, data_t> rest_coords_tensor = 0.0;
    FOR(i, j)
    {
        rest_coords_tensor[i] +=
            inverse_boost_matrix[i + 1][j + 1] * boosted_coords_tensor[j];
    }

    // copy boosted coords to keep center
    Coordinates<data_t> rest_coords = a_boosted_coords;
    rest_coords.x = rest_coords_tensor[0];
    rest_coords.y = rest_coords_tensor[1];
    rest_coords.z = rest_coords_tensor[2];

    return rest_coords;
}

// since the new initial slice will be t_boosted = 0, this function returns the
// value of t_rest in the boosted frame
template <class data_t>
data_t get_t_rest(const Coordinates<data_t> &a_boosted_coords,
                  const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto inverse_boost_matrix = inverse_matrix<data_t>(a_velocity);

    Tensor<1, data_t> boosted_coords_tensor;
    boosted_coords_tensor[0] = a_boosted_coords.x;
    boosted_coords_tensor[1] = a_boosted_coords.y;
    boosted_coords_tensor[2] = a_boosted_coords.z;

    data_t out;
    FOR(i) { out += inverse_boost_matrix[0][i + 1] * boosted_coords_tensor[i]; }

    return out;
}

// boost a type (1,0) tensor
// Tensor should already be evaluated at boosted coordinates
template <class data_t>
static Tensor<1, data_t, CH_SPACETIMEDIM>
boost_U(const Tensor<1, data_t, CH_SPACETIMEDIM> &a_rest_U,
        const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto boost_matrix = matrix<data_t>(a_velocity);

    Tensor<1, data_t, CH_SPACETIMEDIM> out = 0.0;
    FOR_ST(i, j) { out[i] += boost_matrix[i][j] * a_rest_U[j]; }

    return out;
}

// boost a type (2,0) tensor
// Tensor should already be evaluated at boosted coordinates
template <class data_t>
static Tensor<2, data_t, CH_SPACETIMEDIM>
boost_UU(const Tensor<2, data_t, CH_SPACETIMEDIM> &a_rest_UU,
         const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto boost_matrix = matrix<data_t>(a_velocity);

    Tensor<2, data_t, CH_SPACETIMEDIM> out = 0.0;
    FOR_ST(i, j, k, l)
    {
        out[i][j] += boost_matrix[i][k] * boost_matrix[j][k] * a_rest_UU[k][l];
    }

    return out;
}

// boost a type (0,1) tensor
// Tensor should already be evaluated at boosted coordinates
template <class data_t>
static Tensor<1, data_t, CH_SPACETIMEDIM>
boost_L(const Tensor<1, data_t, CH_SPACETIMEDIM> &a_rest_L,
        const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto inverse_boost_matrix = inverse_matrix<data_t>(a_velocity);

    Tensor<1, data_t, CH_SPACETIMEDIM> out = 0.0;
    FOR_ST(i, j) { out[i] += a_rest_L[j] * inverse_boost_matrix[j][i]; }
    return out;
}

// boost a type (0,2) tensor
// Tensor should already be evaluated at boosted coordinates
template <class data_t>
static Tensor<2, data_t, CH_SPACETIMEDIM>
boost_LL(const Tensor<2, data_t, CH_SPACETIMEDIM> &a_rest_LL,
         const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto inverse_boost_matrix = inverse_matrix<data_t>(a_velocity);

    Tensor<2, data_t, CH_SPACETIMEDIM> out = 0.0;
    FOR_ST(i, j, k, l)
    {
        out[i][j] += a_rest_LL[k][l] * inverse_boost_matrix[k][i] *
                     inverse_boost_matrix[l][j];
    }

    return out;
}

// boost a type (0,3) tensor
// Tensor should already be evaluated at boosted coordinates
template <class data_t>
static Tensor<3, data_t, CH_SPACETIMEDIM>
boost_LLL(const Tensor<3, data_t, CH_SPACETIMEDIM> &a_rest_LLL,
          const std::array<double, GR_SPACEDIM> &a_velocity)
{
    auto inverse_boost_matrix = inverse_matrix<data_t>(a_velocity);

    Tensor<3, data_t, CH_SPACETIMEDIM> out = 0.0;
    FOR_ST(i, j, k, l)
    {
        FOR_ST(m, n)
        {
            out[i][j][k] += a_rest_LLL[l][m][n] * inverse_boost_matrix[l][i] *
                            inverse_boost_matrix[m][j] *
                            inverse_boost_matrix[n][k];
        }
    }

    return out;
}

} // namespace LorentzBoosts
#endif /* LORENTZBOOSTS_HPP_ */
