def set_default_path():
    """
    Set the default path folder for outputs to the path of this file + outputs
    """

    return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'outputs')