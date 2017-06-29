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

# Used to print an example and kept for historic reasons.



#
# Assumes that the geotiff file has been resized in a way, that every pixel is a vertex "on the moon"
#

# https://github.com/jswhit/pyproj
# https://gis.stackexchange.com/questions/129847/obtain-coordinates-and-corresponding-pixel-values-from-geotiff-using-python-gdal#129857

# https://astrogeology.usgs.gov/search/map/Moon/LRO/LOLA/Lunar_LRO_LrocKaguya_DEMmerge_60N60S_512ppd
# GEOTIFF = '/Users/msei/Downloads/Lunar_LRO_LrocKaguya_DEMmerge_60N60S_512ppd.tif'

# https://astrogeology.usgs.gov/search/map/Moon/LRO/LOLA/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014#open
# https://astrogeology.usgs.gov/search/map/Moon/LRO/LOLA/Lunar_LRO_LOLA_Global_LDEM_118m_Mar2014
GEOTIFF = 'lunar_large.tiff'
# MOON_RADIUS = 1737400
# MAX_MOON_RADIUS = 1745085
# GEOTIFF_SCALE = 0.5

# https://en.wikipedia.org/wiki/Spherical_coordinate_system
# lonlat23d = lambda r, lo, la: (r * math.sin(la) * math.cos(lo), r * math.sin(la) * math.sin(lo), r * math.cos(lo))

to_rad = lambda t: t/180.0 * math.pi


class test04:
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
    moon_radius_scale = 4.0  # half the values from the input file
    inner_moon_radius = 145
    min_radius = 100000000
    max_radius = -100000000

    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

        if self.infile:
            self.srcProj = Proj(self.infile.crs)
            self.t0 = self.infile.affine
            self.t0_inv = ~self.t0
            self.destProj = Proj(proj='latlong', datum='WGS84')

    def sphere2cart(self, sphere, radius):
        lo_rad, la_rad = sphere

        ret = (radius * math.cos(la_rad) * math.cos(lo_rad),
               radius * math.cos(la_rad) * math.sin(lo_rad),
               radius * math.sin(la_rad))  # more z-stable ?
        return ret

    def transform_to_latlong(self, point):
        p = (point[0], point[1]) * self.t0

        lo, la = transform(self.srcProj, self.destProj, p[1], p[0])
        la_rad = to_rad(la)
        lo_rad = to_rad(lo)

        return la_rad, lo_rad

    def radius_for(self, height):
        range = 17787.0 + 21281.0
        ret = self.moon_radius + height/range * self.moon_radius_scale

        self.min_radius = min(self.min_radius, ret)
        self.max_radius = max(self.max_radius, ret)
        return ret

    def write_sphere_corners(self, top_point, top_height, bottom_point, bottom_height):
        outer_top_point = self.sphere2cart(top_point, self.radius_for(top_height))
        inner_top_point = self.sphere2cart(top_point, self.inner_moon_radius)
        outer_bottom_point = self.sphere2cart(bottom_point, self.radius_for(bottom_height))
        inner_bottom_point = self.sphere2cart(bottom_point, self.inner_moon_radius)

        # inner_top_point
        # |
        # +---------+ (outer_top_point)
        # |        /|
        # |   1  /  |
        # |    /    |
        # |  /   2  |
        # |/        |
        # +---------+ (outer_bottom_point)

        stl.write_3_faces(self.outfile, [outer_top_point, outer_bottom_point, inner_bottom_point])
        stl.write_3_faces(self.outfile, [inner_bottom_point, inner_top_point, outer_top_point])

    def convert(self):
        stl.write_solid_start(self.outfile)

        height_map = self.infile.read(1)
        height_for = lambda p: height_map[p[1] % src.height, p[0] % src.width]
        length_vector = lambda v: math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        diff_vector = lambda v0, v1: (v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2])

        step = 1
        min_length_x = 10000000
        max_length_x = -min_length_x
        min_length_y = min_length_x
        max_length_y = -min_length_y

        min_x_degree = to_rad(-180)
        max_x_degree = to_rad(-110)
        min_y_degree = to_rad(-90)
        max_y_degree = to_rad(50)

        for x in range(0, self.infile.width, step):
            print "x ", x, "/", self.infile.width

            for y in range(0, self.infile.height, step):
                top_left      = (x       , y)
                top_right     = (x + step, y)
                bottom_left   = (x       , y + step)
                bottom_right  = (x + step, y + step)

                outer_top_left_r = height_for(top_left)
                top_left = self.transform_to_latlong(top_left)

                outer_top_right_r = height_for(top_right)
                top_right = self.transform_to_latlong(top_right)

                outer_bottom_left_r = height_for(bottom_left)
                bottom_left = self.transform_to_latlong(bottom_left)

                outer_bottom_right_r = height_for(bottom_right)
                bottom_right = self.transform_to_latlong(bottom_right)

                if top_left[0] >= max_x_degree:
                    break

                if bottom_left[1] <= max_y_degree:
                    self.write_sphere_corners(top_left, outer_top_left_r, top_right, outer_top_right_r)
                    break

                if x == 0:
                    self.write_sphere_corners(top_left, outer_top_left_r, bottom_left, outer_bottom_left_r)

                # if x == 0:  # and x > 0:  # top_left[1] >= min_y_degree:
                if self.transform_to_latlong((x + step, y))[0] >= max_x_degree:
                    self.write_sphere_corners(top_right, outer_top_right_r, bottom_right, outer_bottom_right_r)

                outer_top_left = self.sphere2cart(top_left, self.radius_for(outer_top_left_r))
                outer_top_right = self.sphere2cart(top_right, self.radius_for(outer_top_right_r))
                outer_bottom_left = self.sphere2cart(bottom_left, self.radius_for(outer_bottom_left_r))
                outer_bottom_right = self.sphere2cart(bottom_right, self.radius_for(outer_bottom_right_r))

                inner_top_left = self.sphere2cart(top_left, self.inner_moon_radius)
                inner_top_right = self.sphere2cart(top_right, self.inner_moon_radius)
                inner_bottom_left = self.sphere2cart(bottom_left, self.inner_moon_radius)
                inner_bottom_right = self.sphere2cart(bottom_right, self.inner_moon_radius)

                # Write a patch:
                # +---------+
                # |        /|
                # |   1  /  |
                # |    /    |
                # |  /   2  |
                # |/        |
                # +---------+

                # length x
                length = length_vector(diff_vector(outer_top_left, outer_top_right))
                min_length_x = min(min_length_x, length)
                max_length_x = max(max_length_x, length)

                # length y
                length = length_vector(diff_vector(outer_top_right, outer_bottom_right))
                min_length_y = min(min_length_x, length)
                max_length_y = max(max_length_x, length)

                # cw
                # v1 = [top_left, top_right, bottom_left]
                # v2 = [bottom_left, top_right, bottom_right]

                # ccw
                ov1 = [outer_top_right, outer_top_left, outer_bottom_left]
                ov2 = [outer_bottom_right, outer_top_right, outer_bottom_left]
                iv1 = [inner_top_left, inner_top_right, inner_bottom_left]
                iv2 = [inner_bottom_left, inner_top_right, inner_bottom_right]

                stl.write_3_faces(self.outfile, ov1)
                stl.write_3_faces(self.outfile, ov2)
                stl.write_3_faces(self.outfile, iv1)
                stl.write_3_faces(self.outfile, iv2)

        stl.write_solid_end(self.outfile)

        print "min/max vertex length, x: ", min_length_x, max_length_x, " y:", min_length_y, max_length_y
        print "min/max radius ", self.min_radius, self.max_radius, self.max_radius - self.min_radius


if __name__ == '__main__':
    start = time.time()

    with rasterio.open(GEOTIFF) as src:
        with open("test.stl", "w+") as dest:
            test04(src, dest).convert()

    end = time.time()
    print("Took: %g" % (end - start))
