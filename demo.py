"""Backward-compatible headless entry point.

Prefer ``nontunperc.py`` (GUI by default, or ``python nontunperc.py --cli``).
"""

from nontunperc import run_pipeline, PipelineOptions


if __name__ == "__main__":
    run_pipeline(PipelineOptions())
