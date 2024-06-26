/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: reconciliation.cpp
 * Description: 
 *    This file contains the implementation of the classes/methods described in reconciliation.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 ********************************************************************/

#include "../include/h_files/reconciliation.hpp"

#define sign(x) (((x) >= 0.0) ? +1.0 : -1.0)

#define BP_INFINITY 9e5
#define BP_EPSILON 1e-16

#define TANH_MAX 19.07
#define ATANH_MAX 1.0

#define TANH_TABLE_SIZE 16384

namespace reconciliation
{
    decoder::decoder(size_t NoI, int layered_flag, int fast_flag,  const std::string H_file_name, int lifting_factor, int print_flag)
    {
        this->print_flag = print_flag;
        this->layered_flag = layered_flag;
        read_H(H_file_name);

        this->lifting_factor = lifting_factor;
        max_iter = NoI;
        fast_decoder = fast_flag;
        this->inverse_stepsize_tanh = (double)(TANH_TABLE_SIZE) / (TANH_MAX);
        this->inverse_stepsize_atanh = (double)(TANH_TABLE_SIZE) / (ATANH_MAX);

        if (print_flag == 1)
            std::cout << "Decoder initialized. Fast algorithm : " << fast_decoder << std::endl;

        if (fast_flag == 1)
        {
            tanh_table.resize(TANH_TABLE_SIZE, 0);
            atanh_table.resize(TANH_TABLE_SIZE, 0);

            double stepsize_tanh = (double)(TANH_MAX) / TANH_TABLE_SIZE;
            double stepsize_atanh = (double)(ATANH_MAX) / TANH_TABLE_SIZE;

            for (size_t i = 0; i < TANH_TABLE_SIZE; i++)
            {
                tanh_table[i] = std::tanh((i + 0.5) * stepsize_tanh);
                atanh_table[i] = std::atanh((i + 0.5) * stepsize_atanh);
            }
        }
    }

    size_t decoder::get_N()
    {
        return N;
    }

    size_t decoder::get_M()
    {
        return M;
    }

    void decoder::set_rate(double rate)
    {
        
        N = std::floor(4.0 * this->lifting_factor / rate);
        M = N - 4 * this->lifting_factor;
        

        if ((M > CN_degrees.size()) || (N > number_of_VNs))
        {
            std::cout << "Error: Rate is lower than the minimum supported rate 0.01." << std::endl;
            exit(1);
        }
        else if (rate > 0.2)
        {
            std::cout << "Warning: Rate is higher than the maximum supported rate 0.2." << std::endl;
        }

        double code_rate = 1.0 - (double)(M) / N;

        this->rate = code_rate;

    }

    double decoder::get_rate()
    {
        return rate;
    }

    template <typename T>
    std::vector<uint8_t> decoder::get_syndrome(T word) const
    {
        // this function checks the syndrome of the input word.
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
    std::tuple<std::vector<uint8_t>, int, int> decoder::SPA_decoder(const T1 &llr_in, const T2 &syndrome) const
    {
        // allocate memory for the llr_out
        std::vector<double> llr_VN_sum_prev(N, 0);
        std::vector<double> llr_VN_sum_next(N, 0);
        std::vector<double> llr_edge(Connected_VNs.size(), 0);
        std::vector<double> llr_CN_input(max_CN_degree, 0);
        std::vector<size_t> zero_LLR_idx(max_CN_degree, 0);
        std::vector<uint8_t> hard_decision_out(N, 0);
        std::vector<uint8_t> calculated_syndrome(M, 0);
        int syndrome_match_flag = 0;
        //int layered_flag = 1;

        // initialize the llr_VN_sum_prev
        for (size_t VN_idx = 0; VN_idx < N; VN_idx++)
        {
            llr_VN_sum_prev[VN_idx] = (llr_in[VN_idx]);
        }
        // initialize the llr_edge
        for (size_t edge_idx = 0; edge_idx < Connected_VNs.size(); edge_idx++)
        {
            llr_edge[edge_idx] = 0;
        }
        // initialize the llr_CN_input
        for (size_t CN_idx = 0; CN_idx < max_CN_degree; CN_idx++)
        {
            llr_CN_input[CN_idx] = 0;
        }
        // initialize the llr_VN_sum_next
        for (size_t VN_idx = 0; VN_idx < N; VN_idx++)
        {
            llr_VN_sum_next[VN_idx] = llr_in[VN_idx];
        }

        double tanh_edge_multiplication = 0;
        double tanh_edge_multiplication_one_edge = 0;
        int final_iteration = 0;

        for (size_t iteration = 0; iteration < max_iter; iteration++)
        {
            // go over all the check nodes
            for (size_t CN_idx = 0; CN_idx < M; CN_idx++)
            {
                int sign_multiplication = 1;
                // go over all the connected edges

                tanh_edge_multiplication = 1;
                std::vector<double> tanh_edge_vector(CN_degrees[CN_idx], 0);
                int zero_llr_count = 0;

                for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
                {

                    // get the llr of the connected edge
                    if (layered_flag == 0)
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
                        zero_LLR_idx[zero_llr_count] = edge_idx;
                        zero_llr_count++;
                        if (zero_llr_count > 1)
                            break;
                        continue;
                    }

                    // sum all the llr of the connected edges
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
                    for (size_t edge_idx = 0; edge_idx < CN_degrees[CN_idx]; edge_idx++)
                        llr_edge[CN_idx * max_CN_degree + edge_idx] = 0.0;
                }
            }

            // update the llr_VN_sum_prev and llr_VN_sum_next
            if (layered_flag == 0)
            {
                for (size_t VN_idx = 0; VN_idx < N; VN_idx++)
                {
                    llr_VN_sum_prev[VN_idx] = (llr_VN_sum_next[VN_idx]);    
                    llr_VN_sum_next[VN_idx] = llr_in[VN_idx];
                }
            }

            // do the hard decision and check if the decoder is converged
            for (size_t VN_idx = 0; VN_idx < N; VN_idx++)
            {
                hard_decision_out[VN_idx] = (((layered_flag==0) ? llr_VN_sum_prev[VN_idx] : llr_VN_sum_next[VN_idx]) > 0 ? 0 : 1);
            }

            calculated_syndrome = get_syndrome(hard_decision_out);

            final_iteration = iteration + 1;

            if (calculated_syndrome.size() == syndrome.size() && std::equal(calculated_syndrome.begin(), calculated_syndrome.end(), syndrome.begin()))
            {
                syndrome_match_flag = 1;
                break;
            }
        }

        return std::make_tuple(hard_decision_out, final_iteration, syndrome_match_flag);
    }

    template <typename T1, typename T2>
    int decoder::get_decoding_error_count(const T1 &decoded_message, const T2 &true_message) const
    {
        int error_count = 0;
        for (size_t i = 0; i < decoded_message.size(); i++)
        {
            if (decoded_message[i] != true_message[i])
            {
                error_count++;
            }
        }
        return error_count;
    }

    double decoder::tanh_lookup(double a) const
    {
        if (fast_decoder == 1)
        {
            int sign = a >= 0 ? 1 : -1;

            if (a * sign > TANH_MAX)
            {
                return sign;
            }

            int index = (int)(a * sign * inverse_stepsize_tanh);

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
        if (fast_decoder == 1)
        {
            int sign = a >= 0 ? 1 : -1;

            if (a * sign > 1 - BP_EPSILON)
            {
                return sign * BP_INFINITY;
            }

            int index = (int)(a * sign * inverse_stepsize_atanh);

            if (index < TANH_TABLE_SIZE)
            {
                return 2 * atanh_table[index] * sign;
            }

            return 2 * std::atanh(a);
        }
        else
        {
            if (a >= 1.0)
                return 2 * (TANH_MAX);
            else
                return (std::log1p(a) - std::log1p(-a));
        }
    }

    void decoder::read_H(const std::string H_file_name)
    {

        {
            std::ostringstream buffer;

            buffer << PCM_DIR <<"/VN_connections_" << H_file_name << ".txt";

            std::string full_file_name;
            full_file_name = buffer.str();
            std::ifstream is(full_file_name);
            if (is.is_open())
            {
                std::string read_value;
                std::getline(is, read_value);
                number_of_VNs = std::stoi(read_value);
                while (std::getline(is, read_value, ' '))
                {
                    Connected_VNs.push_back(std::stoi(read_value) - 1);
                }
            }else {
                std::cout << "Error: File VN_connections_" << H_file_name <<  ".txt not found." << std::endl;
                exit(1);
            }
        }

        {
            std::ostringstream buffer;

            buffer <<  PCM_DIR <<"/CN_degrees_" << H_file_name << ".txt";

            std::string full_file_name;
            full_file_name = buffer.str();
            std::ifstream is(full_file_name);
            if (is.is_open())
            {
                std::string read_value;
                std::getline(is, read_value);
                max_CN_degree = std::stoi(read_value);
                while (std::getline(is, read_value, ' '))
                {
                    CN_degrees.push_back(std::stoi(read_value));
                }
            } else {
                std::cout << "Error: File CN_degrees_" << H_file_name <<  ".txt not found." << std::endl;
                exit(1);
            }
        }
    }

    MDR::MDR(int dim, int print_flag)
    {
        this->print_flag = print_flag;
        if ((dim == 1) || (dim == 2) || (dim == 4) || (dim == 8))
            dimension = dim;
        else
        {
            std::cout << "Error: Reconciliation dimension must be 1, 2, 4 or 8" << std::endl;
            exit(1);
        }

        if (print_flag == 1)
        {
            std::cout << "MDR dimension is set to " << dimension << std::endl;
        }
    }

    std::vector<double> MDR::normalize(std::vector<double> &x) const
    {
        std::vector<double> normalization_factor(x.size(), 0);
#pragma omp parallel for
        for (size_t i = 0; i < normalization_factor.size() / dimension; i++)
        {
            double temp_sum = 0;

            for (size_t k = 0; k < dimension; k++)
            {
                temp_sum += x[i * dimension + k] * x[i * dimension + k]; //
            }

            for (size_t k = 0; k < dimension; k++)
            {
                normalization_factor[i * dimension + k] = std::sqrt(temp_sum);
                x[i * dimension + k] = x[i * dimension + k] / normalization_factor[i * dimension + k];
            }
        }

        return normalization_factor;
    }

    template <typename T1, typename T2>
    std::vector<double> MDR::multiplication_Alice(const T1 &channel_message, const T2 &q_states) const
    {
        std::vector<double> post_MDR_sequence_Alice(channel_message.size(), 0);
#pragma omp parallel for
        for (size_t i = 0; i < channel_message.size() / dimension; i++)
        {
            std::vector<double> temp_channel_message(dimension, 0);
            std::vector<double> temp_quantum_states(dimension, 0);

            for (size_t k = 0; k < dimension; k++)
            {
                temp_channel_message[k] = channel_message[i * dimension + k];
                temp_quantum_states[k] = q_states[i * dimension + k];
            }

            auto temp_message = compute_cayley_dickson_construction(temp_channel_message, compute_higher_dimension_conjugate(temp_quantum_states));

            for (size_t k = 0; k < dimension; k++)
            {
                post_MDR_sequence_Alice[i * dimension + k] = temp_message[k];
            }
        }

        if (print_flag == 1)
        {
            std::cout << "MDR is completed." << std::endl;
        }

        return post_MDR_sequence_Alice;
    }

    template <typename T1, typename T2>
    std::vector<double> MDR::multiplication_Bob(const T1 &received_quantum_states, const T2 &QRNG_output_bipolar) const
    {
        std::vector<double> post_MDR_sequence_Bob(received_quantum_states.size(), 0);
#pragma omp parallel for
        for (size_t i = 0; i < received_quantum_states.size() / dimension; i++)
        {
            std::vector<double> temp_quantum_states(dimension, 0);
            std::vector<double> temp_QRNG_output(dimension, 0);

            for (size_t k = 0; k < dimension; k++)
            {
                temp_quantum_states[k] = received_quantum_states[i * dimension + k];
                temp_QRNG_output[k] = QRNG_output_bipolar[i * dimension + k];
            }

            auto temp_message = compute_cayley_dickson_construction(temp_QRNG_output, temp_quantum_states);
            for (size_t k = 0; k < dimension; k++)
            {
                post_MDR_sequence_Bob[i * dimension + k] = temp_message[k];
            }
        }

        if (print_flag == 1)
        {
            std::cout << "MDR is completed." << std::endl;
        }

        return post_MDR_sequence_Bob;
    }

    template <typename T>
    std::vector<double> MDR::compute_higher_dimension_conjugate(const T &x) const
    {
        std::vector<double> xstar = x;
        for (size_t i = 0; i < xstar.size(); i++)
        {
            xstar[i] = -xstar[i];
        }
        xstar[0] = -xstar[0];
        return xstar;
    }

    std::vector<double> MDR::compute_cayley_dickson_construction(const std::vector<double> &x, const std::vector<double> &y) const
    {
        int n = x.size();
        std::vector<double> z(n, 0.0);

        if (n == 1)
        {
            z[0] = x[0] * y[0];
            return z;
        }

        int m = n / 2;

        std::vector<double> a(x.begin(), x.begin() + m);
        std::vector<double> b(x.begin() + m, x.end());
        std::vector<double> c(y.begin(), y.begin() + m);
        std::vector<double> d(y.begin() + m, y.end());

        std::vector<double> ac = compute_cayley_dickson_construction(a, c);
        std::vector<double> db = compute_cayley_dickson_construction(compute_higher_dimension_conjugate(d), b);
        std::vector<double> da = compute_cayley_dickson_construction(d, a);
        std::vector<double> bc = compute_cayley_dickson_construction(b, compute_higher_dimension_conjugate(c));

        for (int i = 0; i < m; i++)
        {
            z[i] = ac[i] - db[i];
        }
        for (int i = m; i < n; i++)
        {
            z[i] = da[i - m] + bc[i - m];
        }

        return z;
    }

    QRNG::QRNG(size_t seq_length, size_t NoF, int print_flag)
        : sequence_length(seq_length), num_of_frames(NoF), uniform_dist(0, 1), print_flag(print_flag)
    {
    }

    std::vector<std::vector<int8_t>> QRNG::generate_random_sequence()
    {

        std::vector<std::vector<int8_t>> QRNG_output(num_of_frames, std::vector<int8_t>(sequence_length));

#pragma omp parallel
        {
            thread_local auto rng = generate_randomly_seeded_engine();

#pragma omp parallel for
            for (size_t i = 0; i < num_of_frames; i++)
            {
                for (size_t j = 0; j < sequence_length; j++)
                    QRNG_output[i][j] = uniform_dist(rng);
            }
        }

        if (print_flag == 1)
        {
            std::cout << "QRNG output is generated." << std::endl;
        }

        return QRNG_output;
    }

    std::mt19937 QRNG::generate_randomly_seeded_engine()
    {
        std::vector<unsigned int> random_data(std::mt19937::state_size);
        std::random_device source;
        std::generate(begin(random_data), end(random_data), [&]()
                      { return source(); });
        std::seed_seq seeds(begin(random_data), end(random_data));
        std::mt19937 seededEngine(seeds);
        return seededEngine;
    }

   CRC::CRC(uint32_t polynomial) : poly(polynomial) {}

    template <typename T>
    std::vector<int> CRC::get_CRC_checksum(const T &data)
    {
        // uint32_t poly = 0xEDB88320; // CRC-32 polynomial
        uint32_t crc = 0xFFFFFFFF; // Initialize CRC with all bits set

        for (size_t i = 0; i < data.size(); ++i)
        {
            crc ^= data[i];
            for (size_t j = 0; j < 8; ++j)
            {
                if (crc & 1)
                {
                    crc = (crc >> 1) ^ poly;
                }
                else
                {
                    crc >>= 1;
                }
            }
        }

        crc ^= 0xFFFFFFFF; // Finalize the CRC

        std::vector<int> result(32);
        for (int i = 0; i < 32; ++i)
        {
            result[31 - i] = (crc >> i) & 1;
        }

        return result; // Return the CRC result
    }

    std::vector<double> binary_to_bipolar(const std::vector<std::vector<int8_t>> &QRNG_output, const size_t N, const size_t num_of_decoding_frames, const int MDR_dim)
    {

        std::vector<double> QRNG_output_bipolar(num_of_decoding_frames * N, 0);

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                QRNG_output_bipolar[i * N + j] = (-2 * QRNG_output[i][j] + 1) / std::sqrt(MDR_dim);
            }
        }

        return QRNG_output_bipolar;
    }

    std::vector<std::vector<double>> get_LLR(std::vector<double> synthetic_channel_output, const std::vector<double> &alice_normalization_vector, const std::vector<double> &bob_normalization_vector, double noise_variance, int MDR_dim, size_t num_of_decoding_frames, int N, int print_flag)
    {

#pragma omp parallel for
        for (size_t i = 0; i < synthetic_channel_output.size(); i++)
        {
            double Lc = 2 * alice_normalization_vector[i] * bob_normalization_vector[i] / noise_variance / std::sqrt(MDR_dim);

            synthetic_channel_output[i] = Lc * synthetic_channel_output[i];
        }

        // Reshape the LLRs to the LDPC frames

        std::vector<std::vector<double>> LLR(num_of_decoding_frames, std::vector<double>(N, 0));

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                LLR[i][j] = synthetic_channel_output[i * N + j];
            }
        }

        if (print_flag == 1)
        {
            std::cout << "LLR values are calculated." << std::endl;
        }

        return LLR;
    }

    template <typename T>
    std::vector<std::vector<uint8_t>> get_CRC_values(const std::vector<std::vector<T>> &decoded_frame, CRC &crc, size_t num_of_decoding_frames)
    {

        std::vector<std::vector<uint8_t>> CRC_decoded(num_of_decoding_frames, std::vector<uint8_t>(32, 0));

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            auto temp = crc.get_CRC_checksum(decoded_frame[i]);
            std::copy(temp.begin(), temp.end(), CRC_decoded[i].begin());
        }

        return CRC_decoded;
    }

    std::tuple<double,int> check_CRC(const std::vector<std::vector<uint8_t>> &CRC_decoded, const std::vector<std::vector<uint8_t>> &CRC_QRNG, const std::vector<int> &syndrome_flags, utilities::statistics &stats, size_t num_of_decoding_frames, int print_flag)
    {
        std::vector<double> CRC_mismatch(num_of_decoding_frames, 0);
        int undetected_error = 0;
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {

            if (syndrome_flags[i] == 1) // It is a Codeword
            {
                for (size_t j = 0; j < 32; j++)
                {
                    if (CRC_decoded[i][j] != CRC_QRNG[i][j])
                    {
                        CRC_mismatch[i] = 1;
                        break;
                    }
                }

                if ((CRC_mismatch[i] == 0) && (stats.get_frame_error_count(i) == 1)) // It is a codeword and CRC matches but we have an error
                {
                    undetected_error++;
                }
            }
        }
        double No_CRC_mismatch = std::accumulate(CRC_mismatch.begin(), CRC_mismatch.end(), 0);

        if (print_flag == 1)
        {
            std::cout << "Number of decoded codeword CRC mismatches: " << No_CRC_mismatch << std::endl;
            std::cout << "Number of undetected erroneous codewords: " << undetected_error << std::endl;
        }

        return std::make_tuple(No_CRC_mismatch, undetected_error);
    }

    std::vector<std::vector<uint8_t>> calculate_syndrome(decoder LDPC_decoder, std::vector<std::vector<int8_t>> QRNG_output, size_t num_of_decoding_frames, int print_flag)
    {

        std::vector<std::vector<uint8_t>> syndrome(num_of_decoding_frames, std::vector<uint8_t>(LDPC_decoder.get_M(), 0));
#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            std::vector<uint8_t> temp_syndrome = LDPC_decoder.get_syndrome(QRNG_output[i]);

            for (size_t j = 0; j < temp_syndrome.size(); j++)
            {
                syndrome[i][j] = temp_syndrome[j];
            }
        }

        if (print_flag == 1)
            std::cout << "Syndromes are calculated." << std::endl;

        return syndrome;
    }

    utilities::statistics reconcile(const std::vector<double> &alice_states_input, const std::vector<double> &bob_states_input, const double rate, const double noise_variance, const size_t NoI, const int MDR_dim, int layered_flag, int fast_flag, int print_flag, const std::string H_file_name, int lifting_factor)
    {
        // make a copy of alice and bob states
        std::vector<double> alice_states = alice_states_input;
        std::vector<double> bob_states = bob_states_input;

        // Initialize the reconciliation parameters and objects
        size_t num_of_total_states = alice_states.size();

        reconciliation::decoder LDPC_decoder(NoI, layered_flag, fast_flag,  H_file_name, lifting_factor);
        
        LDPC_decoder.set_rate(rate);

        double actual_rate = LDPC_decoder.get_rate();

        size_t N = LDPC_decoder.get_N();

        size_t num_of_decoding_frames = ((size_t)(std::floor((double)(num_of_total_states) / N)));

        // Make sure we don't use more resources than we need

        size_t num_of_threads = omp_get_num_threads();
#pragma omp parallel
        {
            num_of_threads = omp_get_num_threads();
        }

        if (num_of_threads > num_of_decoding_frames)
        {
            num_of_threads = num_of_decoding_frames;
        }

        omp_set_num_threads(num_of_threads);

        utilities::statistics stats(N, num_of_decoding_frames, MDR_dim, layered_flag);

        stats.set_noise_variance(noise_variance);

        if (print_flag == 1)
        {
            std::cout << "Number of quantum states: " << num_of_total_states << std::endl;
            std::cout << "Number of decoding frames: " << num_of_decoding_frames << std::endl;
            std::cout << "Blocklength: " << N << std::endl;
            std::cout << "Number of unused quantum states:" << num_of_total_states - num_of_decoding_frames * N << std::endl;

            std::cout << std::endl;
            std::cout << std::endl;

            std::cout << "Bob :" << std::endl;
        }

        stats.set_rate(actual_rate);
        stats.set_num_of_quantum_states(num_of_total_states);
        stats.set_number_of_unused_states(num_of_total_states - num_of_decoding_frames * N);

        QRNG qrng(N, num_of_decoding_frames);
        MDR mdr(MDR_dim);
        CRC crc;

        // Generate the random raw key material using QRNG and convert it to bipolar form

        auto QRNG_output = qrng.generate_random_sequence();

        auto QRNG_output_bipolar = binary_to_bipolar(QRNG_output, N, num_of_decoding_frames, MDR_dim);

        // Normalize the quantum states and perform the MDR multiplication at Bobs's side

        auto bob_normalization_vector = mdr.normalize(bob_states);
        auto channel_message = mdr.multiplication_Bob(bob_states, QRNG_output_bipolar);

        // Calculate the syndromes of the QRNG output

        auto syndrome = calculate_syndrome(LDPC_decoder, QRNG_output, num_of_decoding_frames);

        // Perform the MDR multiplication at Alice's side
        if (print_flag == 1)
        {
            std::cout << std::endl;
            std::cout << "Alice :" << std::endl;
        }

        auto alice_normalization_vector = mdr.normalize(alice_states);
        auto synthetic_channel_output = mdr.multiplication_Alice(channel_message, alice_states);

        // Calculate the LLR values

        if (print_flag == 1)
        {
            std::cout << "Decoding:" << std::endl;
        }

        auto LLR = get_LLR(synthetic_channel_output, alice_normalization_vector, bob_normalization_vector, noise_variance, MDR_dim, num_of_decoding_frames, N);

        // Perform the LDPC decoding

        std::vector<std::vector<uint8_t>> decoded_frame(num_of_decoding_frames, std::vector<uint8_t>(N, 0));
        std::vector<int> syndrome_flags(num_of_decoding_frames, 0);

        stats.start_timing();

#pragma omp parallel for
        for (size_t i = 0; i < num_of_decoding_frames; i++)
        {
            auto decoding_out = LDPC_decoder.SPA_decoder(LLR[i], syndrome[i]);
            auto temp = std::get<0>(decoding_out);
            std::copy(temp.begin(), temp.end(), decoded_frame[i].begin());

            stats.set_bit_error_count(i, LDPC_decoder.get_decoding_error_count(decoded_frame[i], QRNG_output[i]));
            if (stats.get_bit_error_count(i) > 0)
                stats.set_frame_error_count(i);

            stats.set_num_of_iterations(i, std::get<1>(decoding_out));
            syndrome_flags[i] = std::get<2>(decoding_out);

            if ((syndrome_flags[i] == 1) && (stats.get_frame_error_count(i) == 1))
            {
                stats.mark_codeword_as_wrong(i);
            }
        }

        stats.end_timing();

        // Calculate the CRC of the decoded frames:
        auto CRC_decoded = get_CRC_values(decoded_frame, crc, num_of_decoding_frames);

        // Calculate the CRC of the QRNG output
        auto CRC_QRNG = get_CRC_values(QRNG_output, crc, num_of_decoding_frames);

        // Calculate the CRC mismatchs and compare it to the frame errors
        auto temp = check_CRC(CRC_decoded, CRC_QRNG, syndrome_flags, stats, num_of_decoding_frames);

        stats.set_count_of_undetectable_error_after_CRC(std::get<1>(temp));
        stats.set_count_of_wrong_codewords_detected_by_CRC(std::get<0>(temp));

        if (print_flag == 1)
        {
            utilities::print display;
            display.print_statistics_report(stats);
        }

        return stats;
    }


    template std::vector<double> MDR::multiplication_Alice<std::vector<double>, std::vector<double>>(
        const std::vector<double>&, const std::vector<double>&) const;
    
    template std::vector<double> MDR::multiplication_Bob<std::vector<double>, std::vector<double>>(
        const std::vector<double>&, const std::vector<double>&) const;

    template std::vector<uint8_t> decoder::get_syndrome<std::vector<uint8_t>>(std::vector<uint8_t>) const;

    template std::tuple<std::vector<uint8_t>, int, int> decoder::SPA_decoder<std::vector<double>,std::vector<uint8_t>>(const std::vector<double>&, const std::vector<uint8_t>&) const;    

    template int decoder::get_decoding_error_count<std::vector<uint8_t>,std::vector<uint8_t>>(const std::vector<uint8_t> &, const std::vector<uint8_t> &) const;

    template std::vector<double> MDR::compute_higher_dimension_conjugate<std::vector<double>>(const std::vector<double> &) const;

    template std::vector<int> CRC::get_CRC_checksum<std::vector<uint8_t>>(const std::vector<uint8_t> &);    

    template std::vector<std::vector<uint8_t>> get_CRC_values<uint8_t>(const std::vector<std::vector<uint8_t>> &, CRC &, size_t); 


}
