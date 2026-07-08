"""``python -m parametric_ibp_lf_reducer`` entry point (pure Python, no Wolfram runtime)."""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
