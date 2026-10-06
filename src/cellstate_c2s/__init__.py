"""Deviation-from-reference gene ordering for Cell2Sentence cell-state prediction."""

from .builder import DeviationSentenceBuilder
from .c2s_adapter import replace_sentences

__all__ = ["DeviationSentenceBuilder", "replace_sentences"]
__version__ = "0.1.0"
