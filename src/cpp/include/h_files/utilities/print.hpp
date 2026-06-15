/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: print.hpp
 * Description: This file defines the "print" class to display the the statistics of information reconciliation.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added operator + for statistics class
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/
#include "statistics.hpp"

#pragma once

namespace utilities
{
  
  /**
   * @brief This class is used to configure the printing of the simulation.
   *
   */
  class print
  {
  private:
    bool print_flag; ///> 1: print the statistics report, 0: no printing

    /**
     * @brief Print a line with a label and a value.
     *
     * @tparam T Type of the value
     * @param label Label of the line
     * @param value Value to be printed
     * @param suffix Suffix of the line
     */
    template<typename T>
    void print_line(const std::string& label , const T& value, const std::string& suffix = "") const;

  public:
    /**
     * @brief Construct a new print object.
     *
     * @param flag true: print the statistics report, false: no printing
     */
    explicit print(bool flag = true);

    /**
     * @brief Print the statistics of the simulation.
     *
     * @param stats Statistics object
     */
    void print_statistics_report(statistics stats) const;
  };

}
