"""Module contains decorators."""

def add2docs(*doc_lines: str, indentation: int =0):
    """Decorator for adding docstrings to functions.

    Args:
        *doc_lines: collection of docstrings by lines.
        indentation(optional): white space buffered to each line.

    Returns:
        Decorator.
    """
    doc_lines = ["\n" + " " * indentation + line for line in doc_lines]

    def decor(func):
        for line in doc_lines:
            func.__doc__ += line
        return func

    return decor
