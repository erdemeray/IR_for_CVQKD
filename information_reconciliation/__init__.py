import os

# Set the PCM directory to the PCM folder bundled inside this package.
# The C++ extension reads LDPC parity-check matrix files from this path
# at runtime
_pcm_dir = os.path.join(os.path.dirname(__file__), "PCM")
if os.path.isdir(_pcm_dir):
    from ._core import set_pcm_dir
    set_pcm_dir(_pcm_dir)

from ._core import *
