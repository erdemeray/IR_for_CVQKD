/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: CRC.cpp
 * Description:
 *    This file contains the implementation of the class "CRC" described in CRC.hpp file.
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

    CRC::CRC(uint32_t polynomial) : poly(polynomial) {}

    template <typename T>
    std::vector<int> CRC::get_CRC_checksum(const T &data)
    {
        uint32_t crc = 0xFFFFFFFF; // Initialize CRC with all bits set

        for (const auto& byte : data) {
            crc ^= byte;
            for (size_t j = 0; j < 8; ++j) {
                crc = (crc & 1) ? (crc >> 1) ^ poly : crc >> 1;
            }
        }

        crc ^= 0xFFFFFFFF; // Finalize the CRC

        std::vector<int> result(32);
        for (int i = 0; i < 32; ++i)
        {
            result[31 - i] = (crc >> i) & 1;
        }

        return result; 
    }

}
