

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
