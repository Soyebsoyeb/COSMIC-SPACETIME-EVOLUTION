# Contributing

Thank you for helping improve this research code.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

## Standards

1. **Scientific honesty first.** Every new quantity must be labelled established, numerical, approximate, model dependent or speculative, and appear in the relevant file under `docs/`.
2. **Never invent references, DOIs, data or observational constraints.** If a coefficient is recalled and unverified, say so in `docs/REFERENCES.md`.
3. **Test against an exact limit** (Einstein-de Sitter, de Sitter, Milne, closed matter universe, and so on) whenever one exists.
4. **Do not loosen a test tolerance to obtain a pass.** If a test fails, determine whether the code or the expectation is wrong and document the reason.
5. **Figures** are registered with the `@figure` decorator, saved as PNG and PDF, accompany a CSV dataset, and carry a caption.
6. No emoji in code, documentation or commit messages.
7. Run `pytest`, `python scripts/validate_numerics.py` and, when physics changes, `python scripts/generate_all_figures.py` before opening a pull request.

## Pull requests

Describe the physics change, the equations used, the reference, and the validation performed. Update `CHANGELOG.md`.
