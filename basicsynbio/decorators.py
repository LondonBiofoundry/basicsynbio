"""Module contains decorators."""

from dataclasses import dataclass


@dataclass
class ArgDescription:
    """Class for documenting arguments."""

    argument: str
    description: str


def addargs2docs(*arg_descriptions: ArgDescription):
    """Decorator for adding argument docstrings to functions.

    Note:
        Arguments for documentation must be given in the function docstring with a blank description.

    Args:
        arg_descriptions: The argument to document and its description.

    Returns:
        Decorator.
    """

    def decor(func):
        for arg_description in arg_descriptions:
            if func.__doc__.find(arg_description.argument + ":") == -1:
                raise ValueError(
                    f"string '{arg_description.argument}:' not found in {func.__name__} __doc__"
                )
            description_start = func.__doc__.find(arg_description.argument + ":") + len(
                arg_description.argument + ":"
            )
            func.__doc__ = (
                func.__doc__[:description_start]
                + " "
                + arg_description.description
                + func.__doc__[description_start:]
            )
        return func

    return decor
