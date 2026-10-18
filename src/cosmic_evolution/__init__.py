"""cosmic_evolution: computational investigation of cosmic spacetime evolution (research prototype)."""
from .config import CosmologyParams
from .background.friedmann import FLRWBackground

__version__ = "0.1.0"
__all__ = ["CosmologyParams", "FLRWBackground", "__version__"]
