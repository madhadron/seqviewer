from numpy import *

def boundaries(xs, centers):
    """Calculate boundaries between bases.

    Given an array *xs* of x positions at which each sample of this
    chromatogram will be displayed, and an array *centers* of indices
    into *xs* telling where the center of each base is called to be,
    return a list of boundaries between each base, starting with 0 to
    the left of the first base, and ending with the last element of
    *xs* to the right of the last base.

    The algorithm is braindead: split the space between two centers in
    half.  This is what programs like 4Peaks do as well.
    """
    bounds = [0.0]
    for cnext, cprev in zip(centers[1:], centers[:-1]):
        bounds.append((xs[cnext]+xs[cprev]) / 2.0)
    bounds.append(xs[-1])


def scale_to(these_centers, other_centers, this_n_samples, other_n_samples):
    """Produce x positions that will align *these_centers* to *other_centers*.

    For example,::

        n_samples = 20
        centers = [3, 7, 12, 17]
        other_centers = [0.5, 1.5, 2.5, 3.5]
        other_n_samples = 4
        print scale_to(centers, other_centers, n_samples, other_n_samples)

    evaluates to [0.0, 0.17, 0.33, 0.5, 0.75, 1.0, 1.25, 1.5, 1.7,
    1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.83, 4.17].
    """
    assert len(these_centers) == len(other_centers)
    assert this_n_samples > these_centers[-1]
    # Fill in the first half space
    xs = [other_centers[0]*i/float(these_centers[0]) for i in range(these_centers[0])]
    # Next fill a single space for each step
    for i in range(1, len(these_centers)):
        spacing = other_centers[i] - other_centers[i-1]
        offset = other_centers[i-1]
        n = these_centers[i] - these_centers[i-1]
        xs.extend([offset + spacing*i/float(n) for i in range(n)])
    # And add the last half spacing
    spacing = other_n_samples - other_centers[-1]
    offset = other_centers[-1]
    n = this_n_samples - these_centers[-1]
    xs.extend([offset + spacing*i/float(n) for i in range(n)])
    return array(xs)
