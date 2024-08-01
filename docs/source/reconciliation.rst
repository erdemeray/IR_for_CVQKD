
Reconciliation Namespace
##############################

This namespace provides classes and functions for information reconciliation in Continuous Variable Quantum Key Distribution (CV-QKD) systems. The reconciliation process involves comparing correlated observations between two parties and correcting errors.

The reconciliation process consists of the following steps:

1. **Random Number Generation:** The :ref:`QRNG class <QRNG_class>` is simulates the random number generation to be used as the raw key material. 
2. **Multi-dimensional Reconciliation:** The :ref:`MDR class <MDR_class>` performs multi-dimensional reconciliation described in [LAB08]_. 
3. **Error Correction:** The :ref:`decoder class <decoder_class>` corrects the errors in key bits using the syndrome of the  algorithm to correct the errors in the key bits.
4. **Cyclic Redundancy Check(CRC):**  The :ref:`CRC class <CRC_class>` calculates the CRC of the key bits. Parties can verify key bit correctness by comparing CRCs.


The :ref:`reconciliation functions <reconciliation_functions>` are used as helpers in the process. 


Reverse Reconciliation
----------------------

The :func:`reconciliation::reconcile` function implements reverse reconciliation as in :numref:`RR` , allowing the user to perform the reconciliation of correlated observations between parties directly.


.. figure:: figures/RR.png
    :width: 80%
    :align: center
    :name: RR

    The reverse reconciliation algorithm employing multi-dimensional reconciliation (MDR). RNG, LLR and **H** stand for random number generator, log-likelihood ratio and parity check matrix of the code respectively. **x** and **y** are the correlated observations of the parties.




Classes
----------------------

.. _decoder_class:

Decoder Class
============================

The decoder class implements a mechanism for decoding key bits using their corresponding syndrome. The key features of the decoder class are as follows:

* **Forward Error Correction Code**: Uses a rate-adaptive, protograph-based, raptor-like LDPC (Low-Density Parity-Check) code [CS24]_.
* **Decoding Algorithm**: Utilizes the sum-product algorithm for efficient decoding.
* **Update Schedules**: Supports both flooding and layered node update approaches.
* **Rate Range**: Designed to accommodate coding rates between 0.01 and 0.2.
* **Performance**: With the `fast` decoder option that uses look-up table in the decoding, the decoder achieves high performance with reduced decoding time.

.. doxygenclass:: reconciliation::decoder
    :members:
    :protected-members:
    :private-members:    


.. _MDR_class:
MDR Class
============================

The MDR class implements multi-dimensional reconciliation as described in [LAB08]_. The key features of the MDR class are as follows:

* **Dimensionality**: Supports multi-dimensional reconciliation with dimensions ranging from 1 to 8.
* **Implementation**: Uses Cayley-Dickson construction for rotations in higher dimensions.

.. doxygenclass:: reconciliation::MDR
    :members:
    :protected-members:
    :private-members:    

.. _CRC_class:
CRC Class
============================

The CRC class calculates the CRC of the key bits. The key features of the CRC class are as follows:

* **Verification**: Parties can verify the correctness of key bits by comparing CRCs.
* **Implementation**: Uses the CRC-32 algorithm for error detection [MFZ18]_.

.. doxygenclass:: reconciliation::CRC
    :members:
    :protected-members:
    :private-members:    

.. _QRNG_class:
QRNG Class
============================


.. doxygenclass:: reconciliation::QRNG
    :members:
    :protected-members:
    :private-members:    



.. _reconciliation_functions:
Functions
---------------

.. doxygennamespace:: reconciliation

References
----------------

.. [LAB08] A. Leverrier, R. Alléaume, J. Boutros, G. Zémor, and P. Grangier, “Multidimensional reconciliation for a continuous-variable quantum key distribution,” *Phys. Rev. Lett.*, vol. 77, no. 4, Apr. 2008.
.. [CS24] E. E. Cil and L. Schmalen, “Rate-adaptive protograph-based raptor-like LDPC code for continuous-variable quantum key distribution,” in *Proc. Advanced Photonic Congress: Signal Processing in Photonic Communications (SPPCom)*, Quebec City, Canada, Jul. 2024.
.. [MFZ18] M. Milicevic, C. Feng, L. M. Zhang, and P. G. Gulak, “Quasi-cyclic multi-edge LDPC codes for long-distance quantum cryptography,” *npj Quantum Information*, vol. 4, no. 1, Apr. 2018.

