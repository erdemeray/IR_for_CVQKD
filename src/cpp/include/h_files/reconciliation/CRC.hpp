/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: CRC.hpp
 * Description: 
 *    This file contains the class "CRC" to perform CRC-32 check on the converged word.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    31/07/2024 - v1.0 - First stable release
 ********************************************************************/

namespace reconciliation
{

          /**
     * @brief This class is used to calculate the cyclic redundancy check (CRC) of a message.
     *
     */
    class CRC
    {
    private:

        ///> Generator polynomial of the CRC
        uint32_t poly; 

    public:
        /**
         * @brief Construct a new CRC object.
         *
         * @param polynomial Generator polynomial of the CRC, default = 0xEDB88320
         */
        explicit CRC(uint32_t polynomial = 0xEDB88320);

        /**
         * @brief Calculate the CRC checksum of a message.
         *
         * @param data Message
         * @return std::vector<int> CRC checksum of the message
         */
        template <typename T>
        std::vector<int> get_CRC_checksum(const T &data);
    };
}
