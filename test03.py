import stl
import numpy as np
import rasterio
import math
import time
import h5py
import random
# from rasterio.warp import reproject

from affine import Affine
from pyproj import Proj, transform

# https://github.com/jswhit/pyproj
# https://gis.stackexchange.com/questions/129847/obtain-coordinates-and-corresponding-pixel-values-from-geotiff-using-python-gdal#129857

# https://astrogeology.usgs.gov/search/map/Moon/LRO/LOLA/Lunar_LRO_LrocKaguya_DEMmerge_60N60S_512ppd
# GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LrocKaguya_DEMmerge_60N60S_512ppd.tif'

# https://astrogeology.usgs.gov/search/map/Moon/LRO/LOLA/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014#open
GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014.tif'
# MOON_RADIUS = 1737400
# MAX_MOON_RADIUS = 1745085
# GEOTIFF_SCALE = 0.5

OUTFILE = './condensed.hdf5'

# https://en.wikipedia.org/wiki/Spherical_coordinate_system
# lonlat23d = lambda r, lo, la: (r * math.sin(la) * math.cos(lo), r * math.sin(la) * math.sin(lo), r * math.cos(lo))


class test03:
    infile = None
    outfile = None
    step_size = (10, 10)  # longitude, latitude
    srcProj = None
    destProj = None
    last_range_window = None
    last_range_image = None
    last_range_window_size = (20000, 20000)
    t0 = None
    t0_inv = None
    moon_radius = 150      # mm
    moon_radius_scale = 2.0  # half the values from the input file

    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

        if self.infile:
            self.srcProj = Proj(self.infile.crs)
            self.t0 = self.infile.affine
            self.t0_inv = ~self.t0
            self.destProj = Proj(proj='latlong', datum='WGS84')

    def sphere2cart(self, sphere, radius):
        la, lo = sphere  # SWAPPED! ATTENTION
        la_rad = la / 180.0 * math.pi
        lo_rad = lo / 180.0 * math.pi

        ret = (radius * math.cos(la_rad) * math.cos(lo_rad),
               radius * math.cos(la_rad) * math.sin(lo_rad),
               radius * math.sin(la_rad))  # more z-stable ?
        return ret

    def window_contains_point(self, point):
        w = self.last_range_window
        if self.last_range_window is None:
            return False
        return w[0][0] <= point[0] < w[0][1] and w[1][0] <= point[1] < w[1][1]

    def radius_for(self, dataset, sphere):
        if dataset is None:
            return self.moon_radius

        entry = dataset[sphere]
        longitude = entry[0]
        latitude = entry[1]
        height = entry[2]

        range = 17787.0+21281.0

        # if longitude != sphere[0] or latitude != sphere[1]:
        #    print "Data mismatch ", longitude, latitude, sphere

        return height/range * self.moon_radius_scale + self.moon_radius

    def convert(self):
        with h5py.File(OUTFILE) as out:
            h5_grp_data = out["data"]
            stl.write_solid_start(self.outfile)

            step_size = 0.1
            longitude_min = int(-180 * 1.0 / step_size)
            longitude_max = int(-longitude_min + 1.0 / step_size)

            # we iterate from + to - as the GEOTIFF files seems to be mapped that way
            latitude_min = int(90 * 1.0 / step_size)
            latitude_max = int(-latitude_min - 1)

            h5_dataset = h5_grp_data["height_map"]
            swap = lambda a, b: (b, a)

            # min_height = 100000000
            # max_height = -min_height
            # for row in h5_dataset:
            #     for cell in row:
            #         min_height = min(min_height, cell[2])
            #         max_height = max(max_height, cell[2])
            #
            # print min_height, max_height

            h5_x, h5_y = 0, 0
            src_x = 0
            for longitude_i in range(longitude_min, longitude_max):  # - to +
                src_y = 0
                h5_y = 0
                print "Longitude ", (longitude_i / 10.0)

                for latitude_i in range(latitude_min, latitude_max, -1):  # + to -
                    longitude = longitude_i / 10.0
                    latitude = latitude_i / 10.0

                    top_left      = (longitude              , latitude)
                    top_right     = (longitude + step_size  , latitude)
                    bottom_left   = (longitude              , latitude + step_size)
                    bottom_right  = (longitude + step_size  , latitude + step_size)

                    top_left_i =      (longitude_i    , latitude_i)
                    top_right_i =     (longitude_i + 1, latitude_i)
                    bottom_left_i =   (longitude_i    , latitude_i + 1)
                    bottom_right_i =  (longitude_i + 1, latitude_i + 1)

                    top_left = self.sphere2cart(top_left, self.radius_for(h5_dataset, top_left_i))
                    top_right = self.sphere2cart(top_right, self.radius_for(h5_dataset, top_right_i))
                    bottom_left = self.sphere2cart(bottom_left, self.radius_for(h5_dataset, bottom_left_i))
                    bottom_right = self.sphere2cart(bottom_right, self.radius_for(h5_dataset, bottom_right_i))

                    # Write a patch:
                    # +---------+
                    # |        /|
                    # |   1  /  |
                    # |    /    |
                    # |  /   2  |
                    # |/        |
                    # +---------+

                    stl.write_3_faces(self.outfile, [top_left, top_right, bottom_left])
                    stl.write_3_faces(self.outfile, [bottom_left, top_right, bottom_right])

            stl.write_solid_end(self.outfile)


if __name__ == '__main__':
    start = time.time()

    with open("test.stl", "w+") as dest:
        test03(None, dest).convert()

    end = time.time()
    print("Took: %g" % (end - start))
