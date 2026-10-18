.PHONY: install test validate figures sweep reproduce notebooks clean

install:
	pip install -e ".[dev]"

test:
	pytest

validate:
	python scripts/validate_numerics.py

figures:
	python scripts/generate_all_figures.py

sweep:
	python scripts/run_parameter_sweep.py

reproduce:
	python scripts/reproduce_results.py

notebooks:
	python scripts/execute_notebooks.py

clean:
	rm -rf build dist .pytest_cache src/*.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
