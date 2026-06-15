/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: MDR.cpp
 * Description:
 *    This file contains the implementation of the class "MDR" described in MDR.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../../include/h_files/reconciliation.hpp"


namespace reconciliation
{
    
    MDR::MDR(int dim, bool print_flag) : print_flag(print_flag)
    {
        if ((dim == 1) || (dim == 2) || (dim == 4) || (dim == 8))
            dimension = dim;
        else
        {
            throw std::runtime_error("Error: Reconciliation dimension must be 1, 2, 4 or 8.");
        }

        if (print_flag)
        {
            std::cout << "MDR dimension is set to " << dimension << std::endl;
        }
    }

    std::vector<double> MDR::normalize(std::vector<double> &x) const
    {
        std::vector<double> normalization_factor(x.size(), 0.0);
#pragma omp parallel for
        for (size_t i = 0; i < normalization_factor.size() / dimension; i++)
        {
            double temp_sum = 0;

            for (size_t k = 0; k < dimension; k++)
            {
                temp_sum += x[i * dimension + k] * x[i * dimension + k]; //
            }

            double sqrt_temp_sum = std::sqrt(temp_sum);

            for (size_t k = 0; k < dimension; k++)
            {
                normalization_factor[i * dimension + k] = sqrt_temp_sum;
                x[i * dimension + k] /= normalization_factor[i * dimension + k];
            }
        }

        return normalization_factor;
    }

    template <typename T1, typename T2>
    std::vector<double> MDR::multiplication_Alice(const T1 &channel_message, const T2 &q_states) const
    {
        std::vector<double> post_MDR_sequence_Alice(q_states.size(), 0.0);
#pragma omp parallel for
        for (size_t i = 0; i < q_states.size() / dimension; i++)
        {
            std::vector<double> temp_channel_message(dimension, 0.0);
            std::vector<double> temp_quantum_states(dimension, 0.0);

            std::copy_n(channel_message.begin() + i * dimension, dimension, temp_channel_message.begin());
            std::copy_n(q_states.begin() + i * dimension, dimension, temp_quantum_states.begin());

            auto temp_message = compute_cayley_dickson_construction(temp_channel_message, compute_higher_dimension_conjugate(temp_quantum_states));

            std::copy_n(temp_message.begin(), dimension, post_MDR_sequence_Alice.begin() + i * dimension);
        }

        if (print_flag)
        {
            std::cout << "MDR is completed." << std::endl;
        }

        return post_MDR_sequence_Alice;
    }

    template <typename T1, typename T2>
    std::vector<double> MDR::multiplication_Bob(const T1 &received_quantum_states, const T2 &QRNG_output_bipolar) const
    {
        std::vector<double> post_MDR_sequence_Bob(QRNG_output_bipolar.size(), 0.0);
#pragma omp parallel for
        for (size_t i = 0; i < QRNG_output_bipolar.size() / dimension; i++)
        {
            std::vector<double> temp_quantum_states(dimension, 0.0);
            std::vector<double> temp_QRNG_output(dimension, 0.0);

            std::copy_n(received_quantum_states.begin() + i * dimension, dimension, temp_quantum_states.begin());
            std::copy_n(QRNG_output_bipolar.begin() + i * dimension, dimension, temp_QRNG_output.begin());

            auto temp_message = compute_cayley_dickson_construction(temp_QRNG_output, temp_quantum_states);
            std::copy_n(temp_message.begin(), dimension, post_MDR_sequence_Bob.begin() + i * dimension);
        }

        if (print_flag)
        {
            std::cout << "MDR is completed." << std::endl;
        }

        return post_MDR_sequence_Bob;
    }

    template <typename T>
    std::vector<double> MDR::compute_higher_dimension_conjugate(const T &x) const
    {
        std::vector<double> xstar = x;
        std::transform(xstar.begin(), xstar.end(), xstar.begin(), std::negate<double>());
        xstar[0] = -xstar[0];
        return xstar;
    }

    std::vector<double> MDR::compute_cayley_dickson_construction(const std::vector<double> &x, const std::vector<double> &y) const
    {
        size_t n = x.size();
        std::vector<double> z(n, 0.0);

        if (n == 1)
        {
            z[0] = x[0] * y[0];
            return z;
        }

        size_t m = n / 2;

        std::vector<double> a(x.begin(), x.begin() + m);
        std::vector<double> b(x.begin() + m, x.end());
        std::vector<double> c(y.begin(), y.begin() + m);
        std::vector<double> d(y.begin() + m, y.end());

        std::vector<double> ac = compute_cayley_dickson_construction(a, c);
        std::vector<double> db = compute_cayley_dickson_construction(compute_higher_dimension_conjugate(d), b);
        std::vector<double> da = compute_cayley_dickson_construction(d, a);
        std::vector<double> bc = compute_cayley_dickson_construction(b, compute_higher_dimension_conjugate(c));

        std::transform(ac.begin(), ac.end(), db.begin(), z.begin(), std::minus<double>());
        std::transform(da.begin(), da.end(), bc.begin(), z.begin() + m, std::plus<double>());

        return z;
    }


}
