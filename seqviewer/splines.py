import numpy as np

def several(x,n):
    return [x for i in range(n)]

def solve_tridiag(a, b, c, d):
    assert len(a) == len(b)-1
    assert len(c) == len(b)-1
    assert len(b) == len(d)
    N = len(b)
    bprime = np.zeros(len(b), dtype=b.dtype)
    bprime[0] = b[0]
    bprime[1:] = b[1:] - a*c/b[:-1]
    dprime = np.zeros(len(d), dtype=d.dtype)
    dprime[0] = d[0]
    dprime[1:] = d[1:] - d[:-1]*a/b[:-1]
    xs = np.zeros(len(dprime), dtype=d.dtype)
    xs[N-1] = dprime[N-1] / bprime[N-1]
    for i in range(N-2, -1, -1):
        xs[i] = (dprime[i] - xs[i+1]*c[i-1])/bprime[i]
    return xs

class Point(object):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)
    def __repr__(self):
        return "Point(%s,%s)" % (self.x,self.y)
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)
    def __mul__(self, other):
        return Point(self.x*other, self.y*other)
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    def __div__(self, other):
        return Point(self.x/other, self.y/other)
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)


def interpolate(points):
    N = len(points)-1

    a = np.array(several(1,N-2) + [2.])
    b = np.array([2] + several(4,N-2) + [7.])
    c = np.array(several(1,N-1))
    d = np.array([points[1]*2 + points[0]] + \
                     [points[i]*4 + points[i+1]*2 for i in range(1,N-1)] + \
                     [points[N-1]*8 + points[N]])

    left_control_points = solve_tridiag(a,b,c,d)

    ar = np.array(several(1.0, N-1))
    br = np.array([7.0] + several(4.0, N-2) + [2.0])
    cr = np.array([2.0] + several(1.0, N-2))
    dr = np.array([points[0] + points[1]*8] + \
                      [points[i]*2 + points[i+1]*4 for i in range(1,N-1)] + \
                      [points[N-1]*2 + points[N]])
    right_control_points = solve_tridiag(ar,br,cr,dr)

    return {'left': left_control_points,
            'right': right_control_points}
