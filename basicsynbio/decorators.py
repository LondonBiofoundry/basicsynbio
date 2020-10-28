"""Module contains decorators."""


def add2docs(*doc_lines, indentation=0):
    """Decorator for adding docstrings to functions.

    :param \*doc_lines: collection of docstrings by lines.
    :param int indentation: white space buffered to each line.
    """
    doc_lines = [
        "\n" + " "*indentation + line for line in doc_lines
    ]
    def decor(func):
        for line in doc_lines:
            func.__doc__ += line
        return func
    return decor
