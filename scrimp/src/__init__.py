import getfem as gf
import os


def set_default_path():
    """
    Set the default path folder for outputs to the path of this file + outputs
    """

    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "outputs")


def set_verbose_gf(verbose):
    """
    Set the verbosity level of getfem

    :param verbose: the level
    :type verbose: int
    """

    gf.util_trace_level(verbose)
    gf.util_warning_level(verbose)
