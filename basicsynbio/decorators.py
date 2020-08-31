"""Module contains decorators."""


def add2docs(*doc_lines):
    """Decorator for adding docstrings to functions.

    :param \*doc_lines: collection of docstrings by lines.
    """
    doc_lines = [
        "\n" + line for line in doc_lines
    ]
    def decor(func):
        for line in doc_lines:
            func.__doc__ += line
        return func
    return decor
