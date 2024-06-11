/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: use_case.cpp
 * Description: 
 *    This file contains the implemetation of reverse reconciliation. It can be used to read the quantum states of Alice and Bob from files and perform reconciliation. It prints the statistics of the reconciliation process.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 ********************************************************************/

#include <cfenv>
#include <information_reconciliation>
#include <omp.h>
#include "../include/cxx_argp_parser.h"

int main(int argc, char *argv[])
{
    int NoI = 500; // number of decoding iterations
    int fast_decoding_flag = 1; // 1: fast decoding, 0: normal decoding
    size_t MDR_dim = 8;         // {1, 2, 4, 8} for multi dimensional reconcilation
    unsigned int nThreads = std::thread::hardware_concurrency(); // number of decoding threads
    double rate_value = 0.2; //code rate
    std::string Alice_address = "Alice_-4.55_dB_sigma_square_2.851_N_100000000.txt"; //File location of Alice's states
    std::string Bob_address = "Bob_-4.55_dB_sigma_square_2.851_N_100000000.txt"; //File location of Bob's states
    double sigma_square = 2.851; //noise power
    int layered_decoding_flag = 1; // 1: layered decoding, 0: flooding decoding

    // Process the command line arguments 8 arguments

    cxx_argp::parser parser;
    parser.add_option({nullptr, 'A', "Alice", 0, "Alice states",0}, Alice_address);
    parser.add_option({nullptr, 'B', "Bob", 0, "Bob states",0}, Bob_address);
    parser.add_option({"double", 'R', "Rate", 0, "Rate value",0}, rate_value);
    parser.add_option({"size_t", 'd', "MDR", 0, "MDR dimension",0}, MDR_dim);
    parser.add_option({"int", 'f', "Fast decoding", 0, "Fast decoding flag",0}, fast_decoding_flag);
    parser.add_option({"double", 'n', "Noise power", 0, "Noise power",0}, sigma_square);
    parser.add_option({"int", 'i', "NoI", 0, "Number of decoding iterations",0}, NoI);
    parser.add_option({"int", 't', "NoT", 0, "Number of decoding threads",0}, nThreads);
    parser.add_option({"int", 'l', "Layered decoding", 0, "Layered decoding flag",0}, layered_decoding_flag);

    if (!parser.parse(argc, argv))
    {
        std::cerr << "there was parsing error error - exiting\n";
        return 1;
    }

    Alice_address = "../data/" + Alice_address;
    Bob_address = "../data/" + Bob_address;

    omp_set_num_threads(nThreads);

    auto quantum_states = utilities::read_quantum_states(Alice_address);
    auto received_quantum_states = utilities::read_quantum_states(Bob_address);


    auto stats = reconciliation::reconcile(quantum_states, received_quantum_states, rate_value, sigma_square, NoI, MDR_dim, layered_decoding_flag , fast_decoding_flag);

    utilities::print display(1);
    display.print_statistics_report(stats);    

    return 0;
}