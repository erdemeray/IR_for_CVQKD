
Troubleshooting
=============

- You need to use the compiled Python bindings using the same Python version that was used to compile the library. If you encounter any issues loading the library, make sure that `IR_lib` environment is activated in the terminal and the Python version is the same as the one used to compile the library. You can check the Python version by running the following command:

.. code-block:: bash

    python --version

In case you have multiple Python versions installed in your system and `python` command uses different installation, you can specify the Python installation in `IR_lib` environment by running the following command:

.. code-block:: bash

    ${conda_dir}/envs/IR_lib/bin/python


where `${conda_dir}` is the path to the Conda installation directory.

- If you encounter ".dll not found" error in Python after compilation of the library, you can add the path to the .dll file to the path. You can do this by adding the following lines to your Python script:

.. code-block:: python

    import os
    if os.name == 'nt':
        os.add_dll_directory('C:\\Windows\\SysWOW64')
