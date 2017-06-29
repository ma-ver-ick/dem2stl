import numpy as np
import rasterio
import math
import time
import h5py
import os
import sys
# from rasterio.warp import reproject

from affine import Affine
from pyproj import Proj, transform

GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif'
OUTFILE = './condensed.hdf5'

# If you want a model with 300m diameter (150 radius) you need about the following angle for detail
# 0.5 = 1.3mm
# 0.2 = 0.52mm
# 0.1 = 0.26mm


def advance_til(projSrc, projDest, start_point, diff, T0):
    start_point_t = (start_point[1], start_point[0]) * T0

    start_longitude, start_latitude = transform(projSrc, projDest, start_point_t[0], start_point_t[1])

    end_x = start_point[0]
    if diff[0] > 0:
        for x in range(0, 5000):
            start_point_t = (start_point[1] + x, start_point[0]) * T0
            temp_long, temp = transform(projSrc, projDest, start_point_t[0], start_point_t[1])
            if abs(temp_long - start_longitude) >= diff[0]:
                end_x = x + start_point[0]
                break

    end_y = start_point[1]
    if diff[1] > 0:
        for y in range(0, 5000):
            start_point_t = (start_point[1], start_point[0] + y) * T0
            temp, temp_lat = transform(projSrc, projDest, start_point_t[0], start_point_t[1])
            if abs(temp_lat - start_latitude) >= diff[1]:
                end_y = y + start_point[1]
                break

    return end_x, end_y


def main():
    if os.path.exists(OUTFILE):
        os.remove(OUTFILE)

    with h5py.File(OUTFILE) as out:
        h5_grp_data = out.create_group("data")

        with rasterio.open(GEOTIFF) as src:
            T0 = src.affine  # upper-left pixel corner affine transform
            T0 = T0 * Affine.translation(0.5, 0.5)  # Get affine transform for pixel centres
            projSrc = Proj(src.crs)
            # Function to convert pixel row/column index (from 0) to easting/northing at centre
            # rc2en = lambda r, c: (c, r) * T0
            projDest = Proj(proj='latlong', datum='WGS84')

            last_image = None
            last_image_window = None

            step_size = 0.1
            longitude_min = int(-180 * 1.0 / step_size)
            longitude_max = int(-longitude_min + 1.0 / step_size)
            latitude_min = int(90 * 1.0 / step_size)  # we iterate from + to - as the GEOTIFF files seems to be mapped that way
            latitude_max = int(-latitude_min - 1)

            h5_datatype = np.dtype('(3,)i')
            h5_dataset = h5_grp_data.create_dataset("height_map",
                                                    (abs(longitude_min) + abs(longitude_max),
                                                     abs(latitude_min) + abs(latitude_max)),
                                                    dtype=h5_datatype)

            h5_x, h5_y = 0, 0
            src_x = 0
            for longitude_i in range(longitude_min, longitude_max):  # - to +
                src_y = 0
                h5_y = 0

                for latitude_i in range(latitude_min, latitude_max, -1):  # + to -
                    longitude = longitude_i / 10.0
                    latitude = latitude_i / 10.0

                    sphere = (longitude, latitude)

                    till_point = advance_til(projSrc, projDest, (src_x, src_y), (-1, step_size), T0)

                    # till_point[1] + 1, because x is not advanced (see -1 above) and read will return all height
                    area = src.read(1, window=((src_y, till_point[1] + 1), (src_x, till_point[0] + 1)))

                    if area.size > 0:
                        m = area.mean()
                        h5_dataset[h5_x, h5_y] = [longitude_i, latitude_i, m]
                        # print longitude, latitude, m, (src_x, src_y), till_point
                    else:
                        print longitude, latitude, "miss", (src_x, src_y), till_point

                    src_x, src_y = till_point
                    h5_y += 1

                src_x, src_y = advance_til(projSrc, projDest, (src_x, src_y), (step_size, -1), T0)
                h5_x += 1

                    #         pixel_t = pixel * T0
                    #
                    #         longitude, latitude = transform(projSrc, projDest, pixel_t[0], pixel_t[1])
                    # for col in range(0, src.height, 50):
                    #     for row in range(0, src.width, 50):
                    #     #     area = src.read(1, window=((col, col + 1000), (row, row + 1000)))
                    #     #
                    #     #     min_height = min(min_height, area.min())
                    #     #     max_height = max(max_height, area.max())
                    #     #
                    #     #     print min_height, max_height
                    #
                    #
                    #         print longitude, latitude
                    #     # break

if __name__ == '__main__':
    start = time.time()

    main()

    end = time.time()
    print("Took: %g" % (end - start))