
if __name__ == "__main__" and __package__ is None:
    __package__ = "scilib.math"
    import scilib.math

from copy import deepcopy
from numpy import ones_like, zeros_like, where, ndarray, zeros, array, ones

def _jacobi_poly_root(x, alpha, beta, i, pn, pn1):
    an1 = 2 * (i + 1) * (i + alpha + beta + 1) * (2 * i + alpha + beta)
    an2 = (2 * i + alpha + beta + 1) * (alpha ** 2 - beta ** 2)
    an3 = (2 * i + alpha + beta) * (2 * i + alpha + beta + 1) * (2 * i + alpha + beta + 2)
    an4 = 2 * (i + alpha) * (i + beta) * (2 * i + alpha + beta + 2)
    return ((an2 + an3 * x) * pn - an4 * pn1) / an1

def _jacobi_poly_der_root(x, alpha, beta, i, pn, pn1):
    bn1 = (2 * (i + 1) + alpha + beta) * (1 - x ** 2)
    bn2 = (i + 1) * (alpha - beta - (2 * (i + 1) + alpha + beta) * x)
    bn3 = 2 * (i + 1 + alpha) * (i + 1 + beta)
    return (bn2 * pn + bn3 * pn1) / bn1

def jacobi_poly(x, alpha, beta, n_order, derivative=False, recursive=False, acumulative=False):
    """
    This function calculates the Jacobi polynomial using recursion. The
    recursion algorithm can be found on: Spectral/hp Element Methods for
    Computational Fluid Dynamics, G.E. Karniadakis & S. Sherwin. (Appendix A. pag 585)

    Parameters
    ----------
    x: float, ndarray
        This is the abscissas where it is desired to evaluate the Jacobi polynomial. Input valid
        values are in the region: (-1 1).
    alpha: float
        Polynomial shape parameter, see reference above. It cannot be lower or equal to -1.
    beta: float
        Polynomial shape parameter, see reference above. It cannot be lower or equal to -1.
    n_order: integer
        Order of the polynomial required. In case the derivative is required it is
        also the order of the derivative.
    derivative: bool, optional
        This flags control if the first derivative of the polynomial is required at the same
        time as the polynomial.
    recursive: bool, optional
        This flags allows to return the jacobi polynomial evaluation at each step up to the
        specified n_order.
    acumulative: bool, optional
        This flags allows to return the acumulative jacobi polynomial evaluation at each step up to the
        specified n_order. It cannot be combined with recursive flag.

    Returns
    -------
    poly: float, ndarray
        Polynomial values for the given abscissas.
    poly_der: float, ndarray
        Polynomial derivative for the given abscissas.
    """

    # Check input arguments
    if n_order < 0:
        ValueError('n_order argument cannot be lower than 0. See book of reference.')

    if alpha <= -1:
        raise ValueError('Alpha value cannot be lower or equal to -1.')

    if beta <= -1:
        raise ValueError('Beta value cannot be lower or equal to -1.')
    if (any(where(x <=-1.0)[0]) or any(where(x >= 1.0)[0])) and derivative:
        raise ValueError('Input value is out of bound. Calculation domain is: (-1, 1) not [-1, 1].')

    if recursive and acumulative:
        raise ValueError('Recursive and Acumulative flags cannot be applied at the same time.')

    # Calculate only the Jacobi polynomial
    if not derivative:
        pn1 = 1.0*ones_like(x)
        pn = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x)

        if n_order < 1:
            return pn1

        if not recursive and not acumulative:
            for i in range(1, n_order):
                pn_temp = pn
                pn = _jacobi_poly_root(x, alpha, beta, i, pn, pn1)
                pn1 = pn_temp
        elif acumulative:
            pn_i = deepcopy(pn)
            pn = pn1 + pn
            for i in range(1, n_order):
                pn_temp = pn_i
                pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                pn1 = pn_temp
                pn += pn_i
        else:
            pn_i = deepcopy(pn)
            if type(x) is ndarray:
                if len(x.shape) == 1:
                    pn = zeros((n_order+1, x.shape[0]))
                elif len(x.shape) == 2:
                    pn = zeros((n_order+1, x.shape[0], x.shape[1]))
                else:
                    raise ValueError("Only 1D and 2D matrixes are allowed for the calculation.")
                pn[0, :] = pn1
                pn[1, :] = pn_i
                for i in range(1, n_order):
                    pn_temp = pn_i
                    pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                    pn1 = pn_temp
                    pn[i + 1, :] = deepcopy(pn_i)
            else:
                pn = zeros((n_order+1, ))
                pn[0] = pn1
                pn[1] = pn_i
                for i in range(1, n_order):
                    pn_temp = pn_i
                    pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                    pn1 = pn_temp
                    pn[i+1] = deepcopy(pn_i)
        return pn

    # Evaluate the polynomial an its derivative
    elif derivative:
        pn1 = 1.0*ones_like(x)
        dpn1 = zeros_like(x)
        pn = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x)
        dpn = 0.5 * (alpha + beta + 2.0)*ones_like(x)

        if n_order < 1:
            return pn1, dpn1

        if not recursive and not acumulative:
            for i in range(1, n_order):
                # Polynomial
                pn_temp = pn
                pn = _jacobi_poly_root(x, alpha, beta, i, pn, pn1)
                pn1 = pn_temp

                # Derivative
                dpn = _jacobi_poly_der_root(x, alpha, beta, i, pn, pn1)

        elif acumulative:
            pn_i = deepcopy(pn)
            pn = pn1 + pn
            dpn = deepcopy(dpn) * ones_like(x)
            for i in range(1, n_order):
                # Polynomial
                pn_temp = pn_i
                pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                pn1 = pn_temp
                pn += pn_i

                # Derivative
                dpn_i = _jacobi_poly_der_root(x, alpha, beta, i, pn_i, pn1)
                dpn += dpn_i

        else:
            pn_i = deepcopy(pn)
            dpn2 = dpn
            if type(x) is ndarray:
                if len(x.shape) == 1:
                    pn = zeros((n_order + 1, x.shape[0]))
                    dpn = zeros((n_order + 1, x.shape[0]))
                elif len(x.shape) == 2:
                    pn = zeros((n_order + 1, x.shape[0], x.shape[1]))
                    dpn = zeros((n_order + 1, x.shape[0], x.shape[1]))
                else:
                    raise ValueError("Only 1D and 2D matrixes are allowed for the calculation.")
                pn[0, :] = pn1
                pn[1, :] = pn_i
                dpn[0, :] = dpn1
                dpn[1, :] = dpn2
                for i in range(1, n_order):
                    pn_temp = pn_i
                    pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                    pn1 = pn_temp
                    pn[i + 1, :] = deepcopy(pn_i)
                    dpn_i = _jacobi_poly_der_root(x, alpha, beta, i, pn_i, pn1)
                    dpn[i + 1, :] = dpn_i
            else:
                pn = zeros((n_order + 1,))
                dpn = zeros((n_order + 1,))
                pn[0] = pn1
                pn[1] = pn_i
                dpn[0] = dpn1
                dpn[1] = dpn
                for i in range(1, n_order):
                    pn_temp = pn_i
                    pn_i = _jacobi_poly_root(x, alpha, beta, i, pn_i, pn1)
                    pn1 = pn_temp
                    pn[i + 1] = deepcopy(pn_i)
                    dpn_i = _jacobi_poly_der_root(x, alpha, beta, i, pn_i, pn1)
                    dpn[i+1] = dpn_i

        return pn, dpn
