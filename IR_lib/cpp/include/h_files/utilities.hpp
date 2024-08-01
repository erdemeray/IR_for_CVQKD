/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: utilities.hpp
 * Description: This file defines the auxiliary classes used for monitoring and displaying the statistics of information reconciliation. It also defines auxilary functions to be used in the simulation.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

#include <iostream>
#include <cstdint>
#include <vector>
#include <stack>
#include <cmath>

#include <fstream>
#include <cmath>   
#include <vector>    
#include <algorithm> 
#include <iostream>
#include <memory>
#include <functional>
#include <exception>
#include <cstdlib>
#include <chrono>
#include <string>
#include <random>
#include <omp.h>
#include <thread>
#include <iomanip>
#include <numeric>
#include <tuple>
#pragma once

#include "utilities/statistics.hpp"
#include "utilities/print.hpp"

namespace utilities
{

  /**
   * @brief This function reads the quantum states from a file.
   *
   * @param filename File name
   * @return std::vector<double> Vector of quantum states
   */
  std::vector<double> read_quantum_states(const std::string &filename);

  /**
   * @brief Calculate the capacity of a quantum channel.
   *
   * @param snr Signal-to-noise ratio (SNR) in dB
   * @return double Capacity of the quantum channel
   */
  double capacity(double snr);

  /**
   * @brief Calculate the signal-to-noise ratio (SNR) of a quantum channel given the capacity.
   *
   * @param capacity Capacity of the quantum channel
   * @return double Signal-to-noise ratio (SNR) in dB
   */
  double snr_from_capacity(double capacity);

}
