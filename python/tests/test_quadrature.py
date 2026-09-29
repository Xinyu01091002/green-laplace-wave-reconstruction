import math

import numpy as np

from green_laplace.quadrature import gauss_laguerre_rule


def test_gauss_laguerre_moments():
    nodes, weights = gauss_laguerre_rule(8)
    assert np.all(nodes > 0)
    assert np.all(weights > 0)
    for degree in range(16):
        assert np.isclose(
            np.sum(weights * nodes**degree), math.factorial(degree), rtol=3e-11
        )

