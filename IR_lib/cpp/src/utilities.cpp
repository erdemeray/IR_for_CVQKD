/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: utilities.cpp
 * Description:
 *    This file contains the implementation of the classes/methods described in utilities.hpp file.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include "../include/h_files/utilities.hpp"
#include "utilities/statistics.cpp"
#include "utilities/print.cpp"

namespace utilities
{
    double capacity(double snr)
    {
        return 0.5 * std::log2(1 + std::pow(10.0, snr / 10));
    }

    double snr_from_capacity(double capacity)
    {
        return 10 * std::log10(std::pow(2.0, 2 * capacity) - 1);
    }

    std::vector<double> read_quantum_states(const std::string &filename)
    {
        std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open the file: " + filename);
    }

    std::vector<double> states;
    double state;
    while (file >> state) {
        states.push_back(state);
    }

    return states;
    }

} // namespace utilities