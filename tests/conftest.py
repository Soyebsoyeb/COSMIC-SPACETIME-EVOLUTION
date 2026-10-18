import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cosmic_evolution.background.friedmann import FLRWBackground  # noqa: E402
from cosmic_evolution.config import CosmologyParams  # noqa: E402

# T_cmb = 0.01 K and N_eff = 0 make the radiation density ~1e-15, so analytic matter/Lambda/curvature formulae apply
# without a radiation correction.
NO_RAD = dict(T_cmb=1e-2, N_eff=0.0, omega_b_h2=0.0)


def no_rad(**kw) -> FLRWBackground:
    return FLRWBackground(CosmologyParams(**{**NO_RAD, **kw}))


@pytest.fixture(scope="session")
def bg():
    return FLRWBackground(CosmologyParams())


@pytest.fixture(scope="session")
def params():
    return CosmologyParams()
