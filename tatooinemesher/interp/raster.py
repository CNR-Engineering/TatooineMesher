from scipy import interpolate
import numpy as np


def interp_raster(raster):
    # Get raster geometry
    transform = raster.GetGeoTransform()  # to calculate the pixel location of a geospatial coordinate
    x_left, pixelWidth, _, y_top, _, pixelHeight = transform  # x_left : x_origin  and  y_top : y_origin
    # x_left = transform[0], pixelWidth = transform[1], y_top = transform[3], pixelHeight = transform[5]

    # RasterXSize and RasterYSize: they are as a fraction of the input image size
    nb_cols = raster.RasterXSize  # number of pixels _ column
    nb_rows = raster.RasterYSize  # number of pixels _ row

    x_right = x_left + nb_cols * pixelWidth
    y_bottom = y_top + nb_rows * pixelHeight

    x = np.linspace(x_left, x_right, num=nb_cols)
    y = np.linspace(y_bottom, y_top, num=nb_rows)

    # Get elevation and build interpolation
    z = raster.ReadAsArray()
    z[z < -1e38] = np.nan
    new_z = np.flip(z, 0).T

    return interpolate.RegularGridInterpolator((x, y), new_z, bounds_error=False)
