/**
 * @file particle_utilities.h
 * @brief Header file for definitions of utility functions acting on particles.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef PARTICLE_UTILITIES_H
#define PARTICLE_UTILITIES_H

namespace utilities
{
    /**
     * @brief Type alias for a three-vector.
     * @details A three-vector is a tuple of three doubles representing the
     * x, y, and z components of a vector in 3D space. In the context of this
     * package, this is used to represent the momentum and position of
     * particles. The three-vector is used to simplify the implementation of
     * functions which operate on particles.
     */
    typedef std::tuple<double, double, double> three_vector;

    /**
     * @brief Convert a caf::Proxy<float[3]> to a three_vector;
     * @details caf::Proxy is an annoying type to deal with sometimes.
     * Here we define a conversion between the two.
     * @param proxyArr the caf::Proxy<float[3]> to convert
     * @return the three_vector made of the caf::Proxy's three entries
     */
    three_vector to_three_vector(const caf::Proxy<float[3]>& proxyArr)
    {
      return std::make_tuple(proxyArr[0], proxyArr[1], proxyArr[2]);
    }

    /**
     * @brief Calculates the dot product of two three-vectors.
     * @details The dot product of two three-vectors is calculated as the sum
     * of the products of the corresponding components of the two vectors.
     * @param a the first three-vector.
     * @param b the second three-vector.
     * @return the dot product of the two three-vectors.
     */
    double dot_product(const three_vector & a, const three_vector & b)
    {
        return std::get<0>(a)*std::get<0>(b) + std::get<1>(a)*std::get<1>(b) + std::get<2>(a)*std::get<2>(b);
    }

    /**
     * @brief Calculates the magnitude of a three-vector.
     * @details The magnitude of a three-vector is calculated as the square root
     * of the sum of the squares of the components of the vector.
     * @param a the three-vector to calculate the magnitude of.
     * @return the magnitude of the three-vector.
     */
    double magnitude(const three_vector & a)
    {
        return std::sqrt(std::pow(std::get<0>(a), 2) + std::pow(std::get<1>(a), 2) + std::pow(std::get<2>(a), 2));
    }

    /**
     * @brief Implements the addition operator for three-vectors.
     * @details The addition operator for three-vectors is implemented as the
     * addition of the corresponding components of the two vectors.
     * @param a the first three-vector.
     * @param b the second three-vector.
     * @return the sum of the two three-vectors.
     */
    three_vector add(const three_vector & a, const three_vector & b)
    {
        return std::make_tuple(std::get<0>(a) + std::get<0>(b), std::get<1>(a) + std::get<1>(b), std::get<2>(a) + std::get<2>(b));
    }

    /**
     * @brief Implements the subtraction operator for three-vectors.
     * @details The subtraction operator for three-vectors is implemented as the
     * subtraction of the corresponding components of the two vectors.
     * @param a the first three-vector.
     * @param b the second three-vector.
     * @return the difference of the two three-vectors.
     */
    three_vector subtract(const three_vector & a, const three_vector & b)
    {
        return std::make_tuple(std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b), std::get<2>(a) - std::get<2>(b));
    }

    /**
     * @brief Implements the normalization operator for three-vectors.
     * @details The normalization operator for three-vectors is implemented as
     * the division of each component of the vector by the magnitude of the vector.
     * @param a the three-vector to normalize.
     * @return the normalized three-vector.
     */
    three_vector normalize(const three_vector & a)
    {
        double mag = magnitude(a);
        return std::make_tuple(std::get<0>(a)/mag, std::get<1>(a)/mag, std::get<2>(a)/mag);
    }

    /**
     * @brief Calculates the transverse component of the particle's momentum.
     * @details The transverse component of the momentum is calculated with
     * respect to either the BNB beam direction or the NuMI beam direction. In
     * the case of BNB, this is the z-axis. For NuMI, we have make an
     * assumption that the neutrino direction is defined by the vector between
     * the target and the interaction vertex. The position of the NuMI target
     * in detector coordinates is given by the vector (315.120380, 33.644912,
     * 733.632532) in meters. We then need to offset this by the position of
     * the interaction vertex to get the best estimate of the neutrino direction.
     * A preprocessor macro is used to define the beam direction.
     * @param p the momentum of the particle as a tuple (three-vector).
     * @param vtx the position of the interaction vertex as a tuple (three-vector).
     * @return the transverse momentum (three-vector) of the particle as a
     * tuple.
     */
    three_vector transverse_momentum(three_vector & p, three_vector & vtx)
    {
        three_vector unit;
        if constexpr(!BEAM_IS_NUMI)
            unit = std::make_tuple(0, 0, 1);
        else
        {
            three_vector beam = std::make_tuple(315.120380 + std::get<0>(vtx), 33.644912 + std::get<1>(vtx), 733.632532 + std::get<2>(vtx));
            unit = normalize(beam);
        }
        double scale = dot_product(p, unit);
        return subtract(p, std::make_tuple(scale*std::get<0>(unit), scale*std::get<1>(unit), scale*std::get<2>(unit)));
    }

    /**
     * @brief Calculates the longitudinal component of the particle's momentum.
     * @details The longitudinal component of the momentum is calculated with
     * respect to either the BNB beam direction or the NuMI beam direction. In
     * the case of BNB, this is the z-axis. For NuMI, we have make an
     * assumption that the neutrino direction is defined by the vector between
     * the target and the interaction vertex. The position of the NuMI target
     * in detector coordinates is given by the vector (315.120380, 33.644912,
     * 733.632532) in meters. We then need to offset this by the position of
     * the interaction vertex to get the best estimate of the neutrino
     * direction. A preprocessor macro is used to define the beam direction.
     * @param p the momentum of the particle as a tuple (three-vector).
     * @param vtx the position of the interaction vertex as a tuple (three-vector).
     * @return the longitudinal momentum (three-vector) of the particle as a
     * tuple.
     */
    three_vector longitudinal_momentum(three_vector & p, three_vector & vtx)
    {
        three_vector unit;
        if constexpr(!BEAM_IS_NUMI)
            unit = std::make_tuple(0, 0, 1);
        else
        {
            three_vector beam = std::make_tuple(315.120380 + std::get<0>(vtx), 33.644912 + std::get<1>(vtx), 733.632532 + std::get<2>(vtx));
            unit = normalize(beam);
        }
        double scale = dot_product(p, unit);
        return std::make_tuple(scale*std::get<0>(unit), scale*std::get<1>(unit), scale*std::get<2>(unit));
    }
}
#endif // PARTICLE_UTILITIES_H
