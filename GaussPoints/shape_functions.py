
if __name__ == "__main__" and __package__ is None:
    __package__ = "scilib.math.pde.fem"
    import scilib.math.pde.fem

from numpy import zeros, linspace, array, matmul, transpose, meshgrid, newaxis
from numpy.matlib import repmat
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from ...special_functions import jacobi_poly

def modal_functions_1d(epsilon, n_order):
    """
    Calculate 1D modal functions. It is calculated using the normal
    tensor product expansions. More description about the
    modal set and the collapsed coordinate system can be found on chapter 3 of
    Spectral/hp Element Methods for Computational Fluid Dynamics
    G.E. Karniadakis and S. Sherwin.

    Parameters
    ----------
    epsilon: float, ndarray
        Values of the coordinate 1.
    n_order: int
        Order of the expansion set in the first direction.

    Returns
    -------
    shapes: ndarray
        Values at the specified points (epsilon, eta) of the expansion set up to the orders
        p_order and q_order.
    """

    if n_order < 1:
        raise ValueError("Required order for the tensor expansion must be higher or equal to 1."
                         " See specifications.")
    if len(epsilon.shape) > 1:
        raise ValueError("The shape of epsilon argument must be: (x,).")

    # Calculate initial node mode
    shapes = zeros((n_order+1, epsilon.shape[0]))
    shapes[0, :] = (1 - epsilon) / 2.0

    # Calculate intermediate values if any
    if n_order > 1:
        jacobi_poly_rec = jacobi_poly(epsilon, 1.0, 1.0, n_order - 2, recursive=True)
        if n_order == 2:
            shapes[1,:] = (1-epsilon)*(1+epsilon)/4.0*jacobi_poly_rec
        else:
            shapes[1:-1,:] = repmat((1-epsilon)*(1+epsilon)/4.0, jacobi_poly_rec.shape[0], 1)*jacobi_poly_rec

    # Calculate end mode
    shapes[-1, :] = (1 + epsilon) / 2.0

    return shapes

def modal_functions_2d(epsilon, p_order, eta=None, q_order=None):
    """
    Calculate 2D modal functions. It is calculated using the normal
    tensor product expansions. More description about the
    modal set and the collapsed coordinate system can be found on chapter 3 of
    Spectral/hp Element Methods for Computational Fluid Dynamics
    G.E. Karniadakis and S. Sherwin.

    Parameters
    ----------
    epsilon: float, ndarray
        Values of the coordinate 1.
    p_order: int
        Order of the expansion set in the first direction.
    eta: float, ndarray, optional
        Values of the coordinate 2.
    q_order: int, optional
        Order of the expansion set in the second direction.

    Returns
    -------
    shapes: ndarray
        Values at the specified points (epsilon, eta) of the expansion set up to the orders
        p_order and q_order.
    """
    if q_order is None:
        shape1d = modal_functions_1d(epsilon, p_order)
        if eta is None:
            shapes = zeros((shape1d.shape[0], shape1d.shape[0], epsilon.shape[0], epsilon.shape[0]))
            for i in range(shape1d.shape[0]):
                for j in range(shape1d.shape[0]):
                    shapes[i, j, :, :] = matmul(shape1d[i, :].reshape(-1, 1), shape1d[j, :][newaxis])
        else:
            shape1d_eta = modal_functions_1d(eta, p_order)
            shapes = zeros((shape1d.shape[0], shape1d_eta.shape[0], epsilon.shape[0], eta.shape[0]))
            for i in range(shape1d.shape[0]):
                for j in range(shape1d_eta.shape[0]):
                    shapes[i, j, :, :] = matmul(shape1d[i, :].reshape(-1, 1), shape1d_eta[j, :][newaxis])
    else:
        shape1d = modal_functions_1d(epsilon, p_order)
        if eta is None:
            eta = epsilon
        shape1d_eta = modal_functions_1d(eta, q_order)
        shapes = zeros((shape1d.shape[0], shape1d_eta.shape[0], epsilon.shape[0], eta.shape[0]))
        for i in range(shape1d.shape[0]):
            for j in range(shape1d_eta.shape[0]):
                shapes[i, j, :, :] = matmul(shape1d[i, :].reshape(-1, 1), shape1d_eta[j, :][newaxis])
    return shapes

if __name__ == "__main__":
    # n_order = 5
    # x = linspace(-1.0, 1.0)
    # y = modal_functions_1d(x, n_order)
    # for i in range(n_order+1):
    #     plt.plot(x, y[i, :])
    # plt.show()
    epsilon = linspace(-1.0, 1.0)
    shape_function = modal_functions_2d(epsilon, 4)
    X, Y = meshgrid(epsilon, epsilon)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, shape_function[2, 2, :, :])
    plt.show()
