#
# bin/__init__.py
#

def load_python_module(name, path):
    import sys
    from importlib.machinery import SourceFileLoader
    from importlib.util import spec_from_loader, module_from_spec

    loader = SourceFileLoader(name, path)
    spec = spec_from_loader(loader.name, loader)
    mod = module_from_spec(spec)
    loader.exec_module(mod)

    sys.modules[name] = mod

    return mod
