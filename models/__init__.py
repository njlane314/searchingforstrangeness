import importlib

def load_model(name: str):
    """Dynamically import a model architecture from the models package.

    The target module must live under this package and expose a ``build_model``
    function returning an ``nn.Module`` instance.
    """
    module = importlib.import_module(f".{name}", __name__)
    if not hasattr(module, "build_model"):
        raise AttributeError(f"Model module '{name}' lacks a build_model() function")
    return module.build_model()
