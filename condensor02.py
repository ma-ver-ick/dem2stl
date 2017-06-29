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

# GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LrocKaguya_DEMmerge_60N60S_512ppd.tif'
GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif'
GEOTIFF_CHUCK_SIZE = 50000
OUTFILE = './condensed-22gb.hdf5'


# If you want a model with 300m diameter (150 radius) you need about the following angle for detail
# 0.5 = 1.3mm
# 0.2 = 0.52mm
# 0.1 = 0.26mm


def advance_til(projSrc, projDest, start_point, diff, T0, sphere):
    start_point_t = (start_point[1], start_point[0]) * T0

    start_longitude, start_latitude = transform(projSrc, projDest, start_point_t[0], start_point_t[1])

    # check / validate that the start point is in the same ballpark as we are at the sphere
    if not sphere is None and (abs(start_latitude - sphere[1]) > 0.1):
        return None

    SEARCH_RANGE = 5000

    end_x = start_point[0]
    if diff[0] > 0:
        for x in range(0, SEARCH_RANGE):
            start_point_t = (start_point[1] + x, start_point[0]) * T0
            temp_long, temp = transform(projSrc, projDest, start_point_t[0], start_point_t[1])
            if abs(temp_long - start_longitude) >= diff[0]:
                end_x = x + start_point[0]
                break

    end_y = start_point[1]
    if diff[1] > 0:
        for y in range(0, SEARCH_RANGE):
            start_point_t = (start_point[1], start_point[0] + y) * T0
            temp, temp_lat = transform(projSrc, projDest, start_point_t[0], start_point_t[1])
            if abs(temp_lat - start_latitude) >= diff[1]:
                end_y = y + start_point[1]
                break

    return end_x, end_y


def create_last_image_coordinates(last_image_window, src_x, src_y, till_point, delta_x):
    start_x = src_x - last_image_window[1][0]
    end_x = till_point[0] - last_image_window[1][0] + delta_x
    start_y = src_y - last_image_window[0][0]
    end_y = till_point[1] - last_image_window[0][0] + 1
    return end_x, end_y, start_x, start_y


def point_inside(last_image_window, start_point, till_point, delta_x):
    end_x, end_y, start_x, start_y = create_last_image_coordinates(last_image_window, start_point[0], start_point[1],
                                                                   till_point, delta_x)

    if not last_image_window[1][0] <= start_point[0] <= last_image_window[1][1]:
        return False
    if not last_image_window[0][0] <= start_point[1] <= last_image_window[0][1]:
        return False
    if not 0 <= end_x <= last_image_window[1][1] - last_image_window[1][0]:
        return False
    if not 0 <= end_y <= last_image_window[0][1] - last_image_window[0][0]:
        return False
    return True


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

            # we iterate from + to - as the GEOTIFF files seems to be mapped that way
            latitude_min = int(90 * 1.0 / step_size)
            latitude_max = int(-latitude_min - 1)
            delta_x = advance_til(projSrc, projDest, (0, 0), (step_size, -1), T0, None)[0]  # 26px!

            h5_datatype = np.dtype('(3,)i')
            h5_dataset = h5_grp_data.create_dataset("height_map",
                                                    (abs(longitude_min) + abs(longitude_max),
                                                     # 900+900),
                                                     abs(latitude_min) + abs(latitude_max)),
                                                    dtype=h5_datatype)

            GEOTIFF_CHUCK_SIZE = src.height

            h5_x, h5_y = 0, 0
            src_x = 0
            for longitude_i in range(longitude_min, longitude_max):  # - to +
                src_y = 0
                h5_y = 0
                print "Longitude ", (longitude_i / 10.0)

                for latitude_i in range(latitude_min, latitude_max, -1):  # + to -
                    longitude = longitude_i / 10.0
                    latitude = latitude_i / 10.0

                    sphere = (longitude, latitude)

                    till_point = advance_til(projSrc, projDest, (src_x, src_y), (-1, step_size), T0, None)
                    # if till_point is None:
                    #     print " ! Not in src file", (longitude, latitude)
                    #     h5_y += 1
                    #     continue

                    till_point = (src_x + delta_x, till_point[1])

                    till_point = (min(till_point[0], src.width - 1), min(till_point[1], src.height - 1))
                    if till_point[0] - src_x < 0 or till_point[1] - src_y < 0:
                        break  # continue with next iteration

                    # till_point[1] + 1, because x is not advanced (see -1 above) and read will return all height
                    # area = src.read(1, window=((src_y, till_point[1] + 1), (src_x, till_point[0] + 1)))
                    if last_image is None or not point_inside(last_image_window, (src_x, src_y), till_point, 0):
                        range_y = (src_y, GEOTIFF_CHUCK_SIZE)
                        range_x = (src_x, src_x + delta_x)
                        last_image_window = (range_y, range_x)
                        print " * Loading from geotiff", range_x, range_y
                        last_image = src.read(1, window=last_image_window)

                    end_x, end_y, start_x, start_y = create_last_image_coordinates(last_image_window, src_x, src_y,
                                                                                   till_point, 0)
                    area = last_image[start_y:end_y, start_x:end_x]

                    if area.size > 0:
                        m = area.mean()
                        h5_dataset[h5_x, h5_y] = [longitude, latitude, m]
                        # print longitude, latitude, m, (src_x, src_y), till_point
                    else:
                        print longitude, latitude, "miss", (src_x, src_y), till_point

                    src_y = till_point[1]
                    h5_y += 1

                src_x, src_y = advance_til(projSrc, projDest, (src_x, src_y), (step_size, -1), T0, None)
                h5_x += 1
                out.flush()


if __name__ == '__main__':
    start = time.time()

    main()

    end = time.time()
    print("Took: %g" % (end - start))
