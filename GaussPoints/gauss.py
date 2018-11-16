
if __name__ == "__main__" and __package__ is None:
    __package__ = "scilib.math.integrate"
    import scilib.math.integrate

from numpy import zeros, cos, pi, ndarray, array, concatenate, where, linspace
from numpy import abs as np_abs
from ..special_functions import jacobi_poly
from ..._exception_handling import MaxIterError

def _generate_gauss_raw(n_order):
    """
    This function generates a module with the gauss points for the
    Legendre, Radau and  Lobatto quadratures. In tha way it is possible
    to use a gauss points for integration in a faster way

    Parameters
    ----------
    n_order: int
        Desired top order of integration.
    """

    with open("gauss_points_raw.py", "w") as file:
        # Import modules
        file.write("\n")
        file.write("from numpy import array\n")
        file.write("\n")

        # Write Gauss-legendre points
        file.write("def get_gauss_legendre_points(n_order):\n")
        file.write("\n")
        file.write("    if n_order < 1:\n")
        file.write("        raise ValueError('n_order must be equal or higher to 1.')\n")
        file.write("\n")
        file.write("    if n_order == 1:\n")
        count = 1
        while True:
            roots, weights = gauss_points(count)
            numbers_string = "%0.16f"
            for j in range(1, roots.shape[0]):
                numbers_string += ", %0.16f"
            roots_string = "        roots = array([" + numbers_string + "])\n"
            weights_string = "        weights = array([" + numbers_string + "])\n"
            file.write(roots_string % tuple(roots))
            file.write(weights_string % tuple(weights))

            if count >= n_order:
                break
            count += 1
            file.write("    elif n_order == %d:\n" % count)
        file.write("    else:\n")
        file.write("        raise ValueError('n_order specified is not available.')\n")
        file.write("\n")
        file.write("    return roots, weights\n")
        file.write("\n")
        file.write("\n")

        # Write Gauss Lobatto Points
        file.write("def get_gauss_lobatto_points(n_order):\n")
        file.write("\n")
        file.write("    if n_order < 2:\n")
        file.write("        raise ValueError('n_order must be equal or higher to 2.')\n")
        file.write("\n")
        file.write("    if n_order == 2:\n")
        count = 2
        while True:
            roots, weights = gauss_lobatto_points(count)
            numbers_string = "%0.16f"
            for j in range(1, roots.shape[0]):
                numbers_string += ", %0.16f"
            roots_string = "        roots = array([" + numbers_string + "])\n"
            weights_string = "        weights = array([" + numbers_string + "])\n"
            file.write(roots_string % tuple(roots))
            file.write(weights_string % tuple(weights))

            if count >= n_order:
                break
            count += 1
            file.write("    elif n_order == %d:\n" % count)
        file.write("    else:\n")
        file.write("        raise ValueError('n_order specified is not available.')\n")
        file.write("\n")
        file.write("    return roots, weights\n")
        file.write("\n")
        file.write("\n")

        # Write Gauss Radau Points
        file.write("def get_gauss_radau_points(n_order):\n")
        file.write("\n")
        file.write("    if n_order < 2:\n")
        file.write("        raise ValueError('n_order must be equal or higher to 2.')\n")
        file.write("\n")
        file.write("    if n_order == 2:\n")
        count = 2
        while True:
            roots, weights = gauss_radau_points(count)
            numbers_string = "%0.16f"
            for j in range(1, roots.shape[0]):
                numbers_string += ", %0.16f"
            roots_string = "        roots = array([" + numbers_string + "])\n"
            weights_string = "        weights = array([" + numbers_string + "])\n"
            file.write(roots_string % tuple(roots))
            file.write(weights_string % tuple(weights))

            if count >= n_order:
                break
            count += 1
            file.write("    elif n_order == %d:\n" % count)
        file.write("    else:\n")
        file.write("        raise ValueError('n_order specified is not available.')\n")
        file.write("\n")
        file.write("    return roots, weights\n")


def gauss_points(n_order):
    """
    Calculate Gauss-Legendre integration points.

    Parameters
    ----------
    n_order: int
        Number of gauss integrations points.

    Returns
    -------
    roots: ndarray
        Roots for Gauss-Legendre integration.
    weights: ndarray
        Weights for the corresponding root of integration.
    """

    # Set up shape parameters to Lagrange polynomial
    alpha = 0
    beta = 0

    # Calculate roots
    roots = jacobi_poly_roots(alpha, beta, n_order)

    # Calculate weights
    [pm, dpm] = jacobi_poly(roots, alpha, beta, n_order, 'derivative')
    weights = 2 / ((1 - roots ** 2) * dpm ** 2)

    return roots, weights


def gauss_radau_points(n_order):
    """
    Calculate Gauss-Radau integration points.

    Parameters
    ----------
    n_order: int
        Number of gauss integrations points.

    Returns
    -------
    roots: ndarray
        Roots for Gauss-Radau integration.
    weights: ndarray
        Weights for the corresponding root of integration.
    """

    # Set up shape parameters to Lagrange polynomial
    alpha = 0
    beta = 0

    # Calculate roots
    roots = concatenate((array([-1.0]), jacobi_poly_roots(alpha, beta + 1, n_order - 1)))

    # Calculate weights
    pm = jacobi_poly(roots, alpha, beta, n_order - 1)
    weights = (1 - roots) / (n_order **  2.0 * pm ** 2.0)

    return roots, weights


def gauss_lobatto_points(n_order):
    """
    Calculate Gauss-Lobatto integration points.

    Parameters
    ----------
    n_order: int
        Number of gauss integrations points.

    Returns
    -------
    roots: ndarray
        Roots for Gauss-Lobatto integration.
    weights: ndarray
        Weights for the corresponding root of integration.
    """

    # Set up shape parameters to Lagrange polynomial
    alpha = 0
    beta = 0

    # Calculate roots
    roots = concatenate((array([-1.0]), jacobi_poly_roots(alpha+1, beta+1, n_order-2), array([1.0])))

    # Calculate weights
    pm = jacobi_poly(roots, alpha, beta, n_order - 1)
    weights = 2.0 / (n_order * (n_order - 1.0) * pm ** 2.0)

    return roots, weights


def jacobi_poly_roots(alpha, beta, n_order, MAX_ITER=100, TOL=1e-16):
    """
    This function calculates the jacobi polynomial roots. It uses an iterative
    method based on the Newton-Rapshon algorithm. The complete description of
    the algorithm used for the calculation of the roots can be found on:
    Spectral/hp Element Methods for Computational Fluid Dynamics
    G.E. Karniadakis and S. Sherwin. Appendix B. pag 598.

    Parameters
    ----------
    alpha: int
        First shape parameter of the jacobi polynomial.
    beta: int
        Second shape parameter of the jacobi polynomial.
    n_order: int
        Order for the Jacobi polynomial to calculate its roots.
    MAX_ITER: int, optional
        Maximum number of iterations for the convergence of the iterative
        process. See reference above.
    TOL: float, optional
        Error tolerance for the convergence problem. By default machine precision
        is set.
    """

    # Allocate roots vector
    roots = zeros((n_order, ))

    # Find first roots estimation
    x = linspace(-1.0, 1.0, 1000)
    y = jacobi_poly(x, alpha, beta, n_order)
    pos1 = where(((y[1:] > 0.0) * (y[0:-1] < 0.0)) + ((y[1:] < 0.0) * (y[0:-1] > 0.0)))
    roots_est = x[pos1]

    # Calculate roots
    for k in range(n_order):
        # Calculate initial guess(Average based on Chebyshev polys)
        # r = -cos((2 * k + 1) * pi / 2.0 / n_order)
        r = roots_est[k]

        # Newton - Raphson iteration
        count = 0
        while True:
            s = 0.0
            [pm, dpm] = jacobi_poly(r, alpha, beta, n_order, derivative=True)
            delta = -pm / (dpm - pm * s)
            r = r + delta
            if abs(delta) < TOL:
                break

            if count > MAX_ITER:
                raise MaxIterError('Maximum number of iterations reached.')
            count = count + 1

        roots[k] = r

    # Check in case of asymmetric modes (alpha != 0 or beta != 0) if the are repeated roots (usually the first one)


    return roots

if __name__ == "__main__":
    _generate_gauss_raw(20)
