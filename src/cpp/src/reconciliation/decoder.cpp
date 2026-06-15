/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: decoder.cpp
 * Description:
 *    This file contains the implementation of the class "decoder" described in decoder.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../../include/h_files/reconciliation.hpp"

constexpr double sign(double x) { return (x >= 0.0) ? +1.0 : -1.0; }

constexpr double BP_INFINITY = 9e5;
constexpr double BP_EPSILON = 1e-16;

constexpr double TANH_MAX = 19.07;
constexpr double ATANH_MAX = 1.0;

constexpr size_t TANH_TABLE_SIZE = 16384;

namespace reconciliation
{
    decoder::decoder(size_t NoI, bool layered_flag, bool fast_flag, const std::string &H_file_name, int lifting_factor, bool print_flag)
        : print_flag(print_flag), layered_flag(layered_flag), lifting_factor(lifting_factor), max_iter(NoI), fast_decoder(fast_flag), inverse_stepsize_tanh(static_cast<double>(TANH_TABLE_SIZE) / TANH_MAX), inverse_stepsize_atanh(static_cast<double>(TANH_TABLE_SIZE) / ATANH_MAX)
    {
        read_H(H_file_name);

        if (print_flag == 1)
            std::cout << "Decoder initialized. Fast algorithm : " << std::boolalpha << fast_decoder << std::endl;

        if (fast_flag)
        {
            tanh_table.resize(TANH_TABLE_SIZE, 0);
            atanh_table.resize(TANH_TABLE_SIZE, 0);

            double stepsize_tanh = TANH_MAX / TANH_TABLE_SIZE;
            double stepsize_atanh = ATANH_MAX / TANH_TABLE_SIZE;

            for (size_t i = 0; i < TANH_TABLE_SIZE; i++)
            {
                tanh_table[i] = std::tanh((i + 0.5) * stepsize_tanh);
                atanh_table[i] = std::atanh((i + 0.5) * stepsize_atanh);
            }
        }
    }

    size_t decoder::get_N() const
    {
        return N;
    }

    size_t decoder::get_M() const
    {
        return M;
    }

    void decoder::set_rate(double rate)
    {

        N = static_cast<size_t>(std::floor(4.0 * this->lifting_factor / rate));
        M = N - 4 * this->lifting_factor;

        if ((M > CN_degrees.size()) || (N > number_of_VNs))
        {
            throw std::runtime_error("Error: Rate is lower than the minimum supported rate 0.01.");
        }
        else if (rate > 0.2)
        {
            std::cerr << "Warning: Rate is higher than the maximum supported rate 0.2." << std::endl;
        }

        this->rate = 1.0 - static_cast<double>(M) / N;
    }

    double decoder::get_rate() const
    {
        return rate;
    }

    template <typename T>
    std::vector<uint8_t> decoder::get_syndrome(const T &word) const
    {
        std::vector<uint8_t> syndrome(M, 0);

        for (size_t CN_idx = 0; CN_idx < M; CN_idx++)
        {
            uint8_t sum = 0;
            for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
            {
                sum += word[Connected_VNs[CN_idx * max_CN_degree + edge_idx]];
            }

            syndrome[CN_idx] = (sum % 2);
        }
        return syndrome;
    }

    template <typename T1, typename T2>
    std::tuple<std::vector<uint8_t>, int, bool> decoder::SPA_decoder(const T1 &llr_in, const T2 &syndrome) const
    {
        // allocate memory for the llr_out
        std::vector<double> llr_VN_sum_prev(N, 0);
        std::vector<double> llr_VN_sum_next(N, 0);
        std::vector<double> llr_edge(Connected_VNs.size(), 0);
        std::vector<double> llr_CN_input(max_CN_degree, 0);
        std::vector<size_t> zero_LLR_idx(max_CN_degree, 0);
        std::vector<uint8_t> hard_decision_out(N, 0);
        std::vector<uint8_t> calculated_syndrome(M, 0);
        bool syndrome_match_flag = false;

        std::copy(llr_in.begin(), llr_in.end(), llr_VN_sum_prev.begin());
        std::copy(llr_in.begin(), llr_in.end(), llr_VN_sum_next.begin());

        double tanh_edge_multiplication = 0;
        double tanh_edge_multiplication_one_edge = 0;
        int final_iteration = 0;

        for (size_t iteration = 0; iteration < max_iter; iteration++)
        {
            // go over all the check nodes
            for (size_t CN_idx = 0; CN_idx < M; CN_idx++)
            {
                int sign_multiplication = 1;
                tanh_edge_multiplication = 1;

                std::vector<double> tanh_edge_vector(CN_degrees[CN_idx], 0);
                int zero_llr_count = 0;

                // go over all the connected edges
                for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
                {

                    // get the llr of the connected edge
                    if (!layered_flag)
                        llr_CN_input[edge_idx] = llr_VN_sum_prev[Connected_VNs[CN_idx * max_CN_degree + edge_idx]] - llr_edge[CN_idx * max_CN_degree + edge_idx];
                    else
                    {
                        llr_CN_input[edge_idx] = llr_VN_sum_next[Connected_VNs[CN_idx * max_CN_degree + edge_idx]] - llr_edge[CN_idx * max_CN_degree + edge_idx];
                        llr_VN_sum_next[Connected_VNs[CN_idx * max_CN_degree + edge_idx]] -= llr_edge[CN_idx * max_CN_degree + edge_idx];
                    }

                    sign_multiplication *= ((llr_CN_input[edge_idx] > 0) ? 1 : -1);
                    if (llr_CN_input[edge_idx] != 0)
                        tanh_edge_vector[edge_idx] = tanh_lookup(0.5 * std::abs(llr_CN_input[edge_idx]));
                    else
                    {
                        zero_LLR_idx[zero_llr_count++] = edge_idx;
                        if (zero_llr_count > 1)
                            break;
                        continue;
                    }

                    tanh_edge_multiplication *= tanh_edge_vector[edge_idx];
                }

                sign_multiplication *= (-2 * syndrome[CN_idx] + 1); // syndrome decoder adjustment

                // go over all the connected edges
                if (zero_llr_count == 0)
                {
                    for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
                    {
                        tanh_edge_multiplication_one_edge = atanh_lookup(tanh_edge_multiplication / tanh_edge_vector[edge_idx]);

                        llr_edge[CN_idx * max_CN_degree + edge_idx] = ((tanh_edge_multiplication_one_edge) * (sign_multiplication * (llr_CN_input[edge_idx] > 0 ? 1 : -1)));

                        // update the llr of the connected variable node
                        llr_VN_sum_next[Connected_VNs[CN_idx * max_CN_degree + edge_idx]] += llr_edge[CN_idx * max_CN_degree + edge_idx];
                    }
                }
                else if (zero_llr_count == 1)
                {
                    for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
                    {

                        if (edge_idx == zero_LLR_idx[0])
                        {
                            llr_edge[CN_idx * max_CN_degree + edge_idx] = atanh_lookup(tanh_edge_multiplication) * (sign_multiplication * (llr_CN_input[edge_idx] > 0 ? 1 : -1));

                            // update the llr of the connected variable node
                            llr_VN_sum_next[Connected_VNs[CN_idx * max_CN_degree + edge_idx]] += llr_edge[CN_idx * max_CN_degree + edge_idx];
                        }
                        else
                        {
                            llr_edge[CN_idx * max_CN_degree + edge_idx] = 0.0;
                        }
                    }
                }
                else
                {
                    std::fill(llr_edge.begin() + CN_idx * max_CN_degree, llr_edge.begin() + (CN_idx + 1) * max_CN_degree, 0.0);
                }
            }

            // update the llr_VN_sum_prev and llr_VN_sum_next
            if (!layered_flag)
            {
                std::copy(llr_VN_sum_next.begin(), llr_VN_sum_next.end(), llr_VN_sum_prev.begin());
                std::copy(llr_in.begin(), llr_in.end(), llr_VN_sum_next.begin());
            }

            // do the hard decision and check if the decoder is converged
            for (size_t VN_idx = 0; VN_idx < N; VN_idx++)
            {
                hard_decision_out[VN_idx] = (((!layered_flag) ? llr_VN_sum_prev[VN_idx] : llr_VN_sum_next[VN_idx]) > 0 ? 0 : 1);
            }

            calculated_syndrome = get_syndrome(hard_decision_out);

            final_iteration = iteration + 1;

            if (calculated_syndrome.size() == syndrome.size() && std::equal(calculated_syndrome.begin(), calculated_syndrome.end(), syndrome.begin()))
            {
                syndrome_match_flag = true;
                break;
            }
        }

        return std::make_tuple(hard_decision_out, final_iteration, syndrome_match_flag);
    }

    template <typename T1, typename T2>
    int decoder::get_decoding_error_count(const T1 &decoded_message, const T2 &true_message) const
    {
        return std::inner_product(decoded_message.begin(), decoded_message.end(), true_message.begin(), 0,
                                  std::plus<int>(), std::not_equal_to<int>());
    }

    double decoder::tanh_lookup(double a) const
    {
        if (fast_decoder)
        {
            int sign = a >= 0 ? 1 : -1;

            if (a * sign > TANH_MAX)
            {
                return sign;
            }

            size_t index = static_cast<int>(a * sign * inverse_stepsize_tanh);

            if (index < TANH_TABLE_SIZE)
            {
                return tanh_table[index] * sign;
            }

            return sign;
        }
        else
        {
            return std::tanh(a);
        }
    }

    double decoder::atanh_lookup(double a) const
    {
        size_t sign = a >= 0 ? 1 : -1;

        if (fast_decoder)
        {
            if (a * sign > 1 - BP_EPSILON)
            {
                return sign * BP_INFINITY;
            }

            int index = static_cast<int>(a * sign * inverse_stepsize_atanh);

            if (index < TANH_TABLE_SIZE)
            {
                return 2 * atanh_table[index] * sign;
            }

            return 2 * std::atanh(a);
        }
        else
        {
            if (a * sign >= 1.0)
                return 2 * sign * TANH_MAX;
            else
                return (std::log1p(a) - std::log1p(-a));
        }
    }

    void decoder::read_H(const std::string &H_file_name)
    {

        auto read_file = [](const std::string &file_name, std::vector<size_t> &output)
        {
            std::ifstream is(file_name);
            if (!is.is_open())
            {
                throw std::runtime_error("Error: File " + file_name + " not found.");
            }

            std::string read_value;
            std::getline(is, read_value);
            size_t first_value = std::stoul(read_value);

            while (std::getline(is, read_value, ' '))
            {
                output.push_back(std::stoul(read_value));
            }

            return first_value;
        };

        std::string vn_file = PCM_DIR "/VN_connections_" + H_file_name + ".txt";
        number_of_VNs = read_file(vn_file, Connected_VNs);
        std::transform(Connected_VNs.begin(), Connected_VNs.end(), Connected_VNs.begin(), [](size_t x)
                       { return x - 1; });

        std::string cn_file = PCM_DIR "/CN_degrees_" + H_file_name + ".txt";
        max_CN_degree = read_file(cn_file, CN_degrees);
    }


}
