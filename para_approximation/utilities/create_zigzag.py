import numpy as np
import matplotlib.pyplot as plt

from definitions import *
from zigzag_utilities import *

np.random.seed(2)


def sample_zigzag_points(lb, ub, L, min_dist=0.01, max_dist=1.0):
    """
        create points for a zigzagging function from lb to ub with maximal Lipschitz constant of L,
        two consecutive points have a minimal distance of min_dist and maximal distance of max_dist
    """
    assert lb <= ub
    assert L > 0, "Lipschitz constant must be positive"
    assert max_dist > 0
    assert min_dist > 0
    assert min_dist <= max_dist

    # sample the first y point to the first x point at lb
    x_points = [lb]
    y_points = [np.random.uniform(-L, L)]

    # sample x and y points until the entire interval [lb, ub] is covered
    while x_points[-1] < ub:
        # sample another x point uniformly in min/max distance, not exceeding the upper bound
        x_next = min(x_points[-1] + np.random.uniform(min_dist, max_dist), ub)

        # compute the y interval to sample from and do so
        x_diff = x_next - x_points[-1]
        y_next = y_points[-1] + np.random.uniform( - L * x_diff, L * x_diff)

        # add both to the corresponding arrays
        x_points.append(x_next)
        y_points.append(y_next)

    return x_points, y_points


def plot_zigzag(x_points, y_points):
    plt.plot(x_points, y_points)
    plt.plot(0, zigzag(0), c='r', marker='x')
    plt.show()


if __name__ == "__main__":
    # sample points
    x_samples, y_samples = sample_zigzag_points(-5, 5, 1)
    # plot results
    plot_zigzag(x_samples, y_samples)
    # write parameter
    zigzag_path = os.path.join(ROOT_DIR, "para_approximation", "utilities", "zigzag_parameters.csv")
    with open(zigzag_path, "w") as file:
        file.write(",".join(str(x) for x in x_samples) + "\n")
        file.write(",".join(str(y) for y in y_samples))