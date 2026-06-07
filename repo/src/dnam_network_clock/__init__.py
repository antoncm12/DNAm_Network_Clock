"""DNAm Network Clock — apply the trained K=60 PCA epigenetic age clock."""
from .predict import predict_age, clock_info

__version__ = "1.0.0"
__all__ = ["predict_age", "clock_info"]
