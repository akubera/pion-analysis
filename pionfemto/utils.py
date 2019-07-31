
def stringify_arg(arg):
    if isinstance(arg, str):
        return '"%s"' % arg
    else:
        return str(arg)


def stringify_arguments(*args):
    return ', ' .join(map(stringify_arg, args))
