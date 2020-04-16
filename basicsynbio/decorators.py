def add2docs(indentation: int, *doc_lines):
    """Decorator for adding docstrings to functions.
    
    Args:
        identation -- white space buffered to each line.
        *doc_lines -- collection of docstrings by lines.
    """
    doc_lines = [
        "\n" + " "*indentation + line for line in doc_lines
    ]
    def decor(func):
        for line in doc_lines:
            func.__doc__ += line
        return func
    return decor