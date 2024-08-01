
Examples
##############################

The Information Reconciliation Library for CV-QKD Systems provides example use cases in C++ and Python named `use_case`, demonstrating the library's usage. The test cases are located in the `tests` folders.

C++
---------------------

An example to use the compiled library in C++ is located at the `IR_lib/cpp/tests` directory. The example assumes that the transmitted and received data are saved in a file. Firstly the data is read from the file and then the reconciliation process is performed using the library. The statistics of the reconciliation process are printed to the console. The script that uses the library is provided in the `src/use_case.cpp` file.

In order to compile the example, the library should be compiled first as described in :ref:`installation`. Then navigate to the `IR_lib/cpp/tests` directory and run the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make

The compiled executable is located in `IR_lib/cpp/tests`.

Python
---------------------

Examples to use the compiled library in Python are located at the `IR_lib/tests` directory. The examples firstly sample random data and transmits them to the receiver over an additive white Gaussian noise (AWGN) channel. The receiver receives the data and performs the reconciliation process using the library. The statistics of the reconciliation process are printed to the console. 

In order to compile the library, the library should be compiled first as described in :ref:`installation`. Then, you can use the Jupyter notebook provided in the `IR_lib/tests` directory to run the examples.

The first example named `QCrypt_script.ipynb` simulates the reconciliation for a fixed reconciliation rate for different channel signal-to-noise ratio values and displays the statistics.

By modifying the `slow` and `fast` decoder options and :class:`MDR` dimension, users can compare the performance of the decoders. The results of the test cases are presented in :numref:`result` for the two extreme rate values supported by the library: rate 0.01 and rate 0.2.

.. figure:: figures/results.png
    :width: 95%
    :align: center
    :name: result

    Results of test cases. FER, BER, and NoI represent frame error rate, bit error rate, and number of decoding iterations, respectively. The results are shown for :class:`MDR` dimensions of 8, 4, 2, and 1, from left to right. The `slow` decoder is denoted by points, while the `fast` decoder is represented by lines. Decoding durations are normalized to the number of parallel decoding operations and were obtained using an AMD EPYC :sup:`TM` 7713P 64-Core Processor.

The second example is named `QCrypt_script_fixed_beta.ipynb` and shows the performance of the reconciliation process for different channel signal-to-noise ratio values for a fixed reconciliation efficiency. The user can also change the modulation of the transmitted states between Gaussian and 4-state (BPSK) modulation and compare the efficiency of the multi-dimensional reconciliation.



Examples
##############################

The Information Reconciliation Library for CV-QKD Systems provides example use cases in C++ and Python, demonstrating the library's usage. These examples are located in the `tests` folders.

C++
---------------------

An example demonstrating the use of the compiled library in C++ is located in the `IR_lib/cpp/tests` directory. This example assumes that the transmitted and received data are saved in a file. The data is first read from the file, and then the reconciliation process is performed using the library. The statistics of the reconciliation process are printed to the console. The script that utilizes the library is provided in the `src/use_case.cpp` file.

To compile the example, the library must first be compiled as described in :ref:`installation`. Then, navigate to the `IR_lib/cpp/tests` directory and execute the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make

The compiled executable is located in the `IR_lib/cpp/tests/bin` directory.

Python
---------------------

Examples demonstrating the use of the compiled library in Python are located in the `IR_lib/tests` directory. These examples first sample random data and transmit them to the receiver over an additive white Gaussian noise (AWGN) channel. The receiver then performs the reconciliation process using the library. The statistics of the reconciliation process are printed to the console.

To compile the library, follow the instructions described in :ref:`installation`. Then, you can use the Jupyter notebook provided in the `IR_lib/tests` directory to run the examples.

The first example, named `QCrypt_script.ipynb`, simulates the reconciliation process for a fixed reconciliation rate across different channel signal-to-noise ratio values and displays the statistics.

By modifying the `slow` and `fast` decoder options and the :class:`MDR` dimension, users can compare the performance of the decoders. The results of the test cases are presented in :numref:`result` for the two extreme rate values supported by the library: rate 0.01 and rate 0.2.

.. figure:: figures/results.png
    :width: 95%
    :align: center
    :name: result

    Results of test cases. FER, BER, and NoI represent frame error rate, bit error rate, and number of decoding iterations, respectively. The results are shown for :class:`MDR` dimensions of 8, 4, 2, and 1, from left to right. The `slow` decoder is denoted by points, while the `fast` decoder is represented by lines. Decoding durations are normalized to the number of parallel decoding operations and were obtained using an AMD EPYC™ 7713P 64-Core Processor.

The second example, named `QCrypt_script_fixed_beta.ipynb`, demonstrates the performance of the reconciliation process for different channel signal-to-noise ratio values at a fixed reconciliation efficiency. Users can also change the modulation of the transmitted states between Gaussian and 4-state (BPSK) modulation and compare the efficiency of the multi-dimensional reconciliation.
