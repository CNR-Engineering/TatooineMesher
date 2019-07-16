import math
import numpy as np


class CubicHermiteSpline:
    """
    CubicHermiteSpline: Generate a cubic Hermite spline from a key points.
    Key points: [[x0, y0], [x1, y1], [x2, y2], ...].
    """
    FINITE_DIFF = 0  # Tangent method: finite difference method
    CARDINAL = 1  # Tangent method: Cardinal spline (c is used)
    ZERO = 0  # End tangent: zero
    GRAD = 1  # End tangent: gradient (m is used)

    class TKeyPoint:
        X = 0.0  # Input
        Y = 0.0  # Output
        M = 0.0  # Gradient

        def __str__(self):
            return '[' + str(self.X) + ', ' + str(self.Y) + ', ' + str(self.M) + ']'

    def __init__(self):
        self.idx_prev = 0

    def find_idx(self, x, idx_prev=0):
        idx = idx_prev
        if idx >= len(self.KeyPts):
            idx = len(self.KeyPts) - 1
        while idx + 1 < len(self.KeyPts) and x > self.KeyPts[idx + 1].X:
            idx += 1
        while idx >= 0 and x < self.KeyPts[idx].X:
            idx -= 1
        return idx

    # Return interpolated value at t
    def evaluate(self, x):
        idx = self.find_idx(x, self.idx_prev)
        if abs(x - self.KeyPts[-1].X) < 1.0e-6:
            idx = len(self.KeyPts) - 2
        if idx < 0 or idx >= len(self.KeyPts) - 1:
            print('WARNING: Given t= %f is out of the key points (index: %i)' % (x, idx))
            if idx < 0:
                idx = 0
                x = self.KeyPts[0].X
            else:
                idx = len(self.KeyPts) - 2
                x = self.KeyPts[-1].X

        h00 = lambda t: t * t * (2.0 * t - 3.0) + 1.0
        h10 = lambda t: t * (t * (t - 2.0) + 1.0)
        h01 = lambda t: t * t * (-2.0 * t + 3.0)
        h11 = lambda t: t * t * (t - 1.0)

        self.idx_prev = idx
        p0 = self.KeyPts[idx]
        p1 = self.KeyPts[idx+1]
        xr = (x - p0.X) / (p1.X - p0.X)
        return h00(xr) * p0.Y + h10(xr) * (p1.X - p0.X) * p0.M + h01(xr) * p1.Y + h11(xr) * (p1.X - p0.X) * p1.M

    def Initialize(self, data, tan_method=CARDINAL, end_tan=GRAD, c=0.0, m=1.0):
        self.KeyPts = [self.TKeyPoint() for i in range(len(data))]
        for idx in range(len(data)):
            self.KeyPts[idx].X = data[idx][0]
            self.KeyPts[idx].Y = data[idx][1]

        grad = lambda idx1, idx2: (self.KeyPts[idx2].Y - self.KeyPts[idx1].Y) / (
                self.KeyPts[idx2].X - self.KeyPts[idx1].X)

        for idx in range(1, len(self.KeyPts) - 1):
            self.KeyPts[idx].M = 0.0
        if tan_method == self.FINITE_DIFF:
            for idx in range(1, len(self.KeyPts) - 1):
                self.KeyPts[idx].M = 0.5 * grad(idx, idx + 1) + 0.5 * grad(idx - 1, idx)
        elif tan_method == self.CARDINAL:
            for idx in range(1, len(self.KeyPts) - 1):
                self.KeyPts[idx].M = (1.0 - c) * grad(idx - 1, idx + 1)
        else:
            raise NotImplementedError

        if end_tan == self.ZERO:
            self.KeyPts[0].M = 0.0
            self.KeyPts[-1].M = 0.0
        elif end_tan == self.GRAD:
            self.KeyPts[0].M = m * grad(0, 1)
            self.KeyPts[-1].M = m * grad(-2, -1)
        else:
            raise NotImplementedError


# def get_interp_cubic_hermite_spline(distances, x_array, y_array, distances_new, tan_method, c=None):
#     spline_x = CubicHermiteSpline()  # x = spline_x(dist)
#     spline_y = CubicHermiteSpline()  # y = spline_y(dist)
#
#     spline_x.Initialize(np.vstack((distances, x_array)).T, tan_method=tan_method, c=c)
#     spline_y.Initialize(np.vstack((distances, y_array)).T, tan_method=tan_method, c=c)
#
#     xx = []
#     yy = []
#     for dist in distances_new:
#         x = spline_x.evaluate(dist)
#         y = spline_y.evaluate(dist)
#         xx.append(x)
#         yy.append(y)
#     return xx, yy
#
#
# def export_points(i, liste_profiles, type_tang):
#     x_array = np.zeros(len(liste_profiles))
#     y_array = np.zeros(len(liste_profiles))
#     distances = np.zeros(len(liste_profiles))
#
#     # Calcul des distances curvilignes
#     x_old, y_old = None, None
#     for j, profile in enumerate(liste_profiles):
#         x = profile.coord.array['X'][i]
#         y = profile.coord.array['Y'][i]
#         x_array[j] = x
#         y_array[j] = y
#         if j != 0:
#             distances[j] = distances[j - 1] + math.sqrt((x - x_old) ** 2 + (y - y_old) ** 2)
#         x_old, y_old = x, y
#
#     distances_new = np.arange(distances.min(), distances.max(), 0.1)
#     distances_new = np.unique(np.concatenate((distances_new, distances)))
#
#     X2, Y2 = get_interp_cubic_hermite_spline(distances, x_array, y_array, distances_new, type_tang)
#     with open('points[%i]_hcs_(%s).i2s' % (i, type_tang), 'w') as fileout:
#         for valx, valy in zip(X2, Y2):
#             fileout.write('%f %f \n' % (valx, valy))
