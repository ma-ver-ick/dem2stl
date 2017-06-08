import stl
import numpy as np
import rasterio
import math
import time
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
MOON_RADIUS = 1737400
MAX_MOON_RADIUS = 1745085
GEOTIFF_SCALE = 0.5

# https://en.wikipedia.org/wiki/Spherical_coordinate_system
# lonlat23d = lambda r, lo, la: (r * math.sin(la) * math.cos(lo), r * math.sin(la) * math.sin(lo), r * math.cos(lo))


class test02:
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
    moon_radius = 1737400  # meters
    moon_radius_scale = 5  # half the values from the input file

    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile

        if self.infile:
            self.srcProj = Proj(self.infile.crs)
            self.t0 = self.infile.affine
            self.t0_inv = ~self.t0
            self.destProj = Proj(proj='latlong', datum='WGS84')

    def sphere2cart(self, sphere, radius):
        lo, la = sphere
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

    def radius_for(self, sphere):
        try:
            x, y = transform(self.destProj, self.srcProj, sphere[0], sphere[1])
            x, y = self.t0_inv * (x, y)

            if not self.window_contains_point((x, y)) or self.last_range_image is None:
                s = self.last_range_window_size
                self.last_range_window = ((x - s[0], y - s[0]), (x + s[0], y + s[1]))
                self.last_range_image = self.infile.read(1, window=self.last_range_window)

            w = self.last_range_window[0]
            pos = (x - w[0], y - w[1])
            v = self.last_range_image[pos]
            return v * self.moon_radius_scale + self.moon_radius
        except:
            print "Error ", sphere
            pass
        return self.moon_radius

    def convert(self):
        stl.write_solid_start(self.outfile)

        for longitude in range(-180, 180 + self.step_size[0], self.step_size[0]):
            for latitude in range(-90, 90 + self.step_size[1], self.step_size[1]):
                # for longitude in np.linspace(-180, 180 - self.step_size[0], self.step_size[0]):
                #    for latitude in np.linspace(-90, 90 - self.step_size[1], self.step_size[1]):
                # print latitude
                # sphere = (longitude, latitude)
                # point3d = self.sphere2cart(sphere, self.radius_for(sphere))
                # print point3d

                top_left =      (longitude                       , latitude)
                top_right =     (longitude + self.step_size[0]   , latitude)
                bottom_left =   (longitude                       , latitude + self.step_size[1])
                bottom_right =  (longitude + self.step_size[0]   , latitude + self.step_size[1])

                top_left = self.sphere2cart(top_left, self.radius_for(top_left))
                top_right = self.sphere2cart(top_right, self.radius_for(top_right))
                bottom_left = self.sphere2cart(bottom_left, self.radius_for(bottom_left))
                bottom_right = self.sphere2cart(bottom_right, self.radius_for(bottom_right))

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

        # for col in range(0, infile.height - PATCH_SIZE, PATCH_SIZE):
        #     for row in range(0, infile.width - PATCH_SIZE, PATCH_SIZE):
        #         if col + PATCH_SIZE >= infile.height or row + PATCH_SIZE >= infile.width:
        #             continue
        #
        #         # transformTo3D(transformToLatLong(projSrc, projDest, swap((row, col)), MOON_RADIUS, T0))
        #         # continue
        #         part_shape = (min(infile.width, row + PATCH_SIZE), min(infile.height, col + PATCH_SIZE))
        #
        #         # _f file/window coordinates
        #         top_left_f = swap((0, 0))  # (row, col)
        #         top_right_f = swap((0, part_shape[1]))  # (row + PATCH_SIZE, col)
        #         bottom_left_f = swap((part_shape[0], 0))  # (row, col + PATCH_SIZE)
        #         bottom_right_f = swap((part_shape[0], part_shape[1]))  # (row + PATCH_SIZE, col + PATCH_SIZE)
        #
        #         # _r real world PIXEL - coordinates
        #         top_left_r = swap((row, col))
        #         top_right_r = swap((row + part_shape[0], col))
        #         bottom_left_r = swap((row, col + part_shape[1]))
        #         bottom_right_r = swap((row + part_shape[0], col + part_shape[1]))
        #
        #         radius_top_left = pixel2m(get_point_from(infile, top_left_r))
        #         radius_top_right = pixel2m(get_point_from(infile, top_right_r))
        #         radius_bottom_left = pixel2m(get_point_from(infile, bottom_left_r))
        #         radius_bottom_right = pixel2m(get_point_from(infile, bottom_right_r))
        #
        #         min_height = min(min_height, radius_top_left, radius_bottom_left, radius_top_right, radius_bottom_right)
        #         max_height = max(max_height, radius_top_left, radius_bottom_left, radius_top_right, radius_bottom_right)
        #
        #         # Write a patch:
        #         # +---------+
        #         # |        /|
        #         # |      /  |
        #         # |    /    |
        #         # |  /      |
        #         # |/        |
        #         # +---------+
        #
        #         temp = transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))
        #
        #         if False:  # clockwise
        #             # top left triangle
        #             writeFaces(outfile, [
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))),
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0)))]
        #                        )
        #
        #             # bottom right triangle
        #             writeFaces(outfile, [
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0))),
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_right_r, radius_bottom_right, T0)))]
        #                        )
        #         else:  # counter-clockwise
        #             # top left triangle
        #             writeFaces(outfile, [
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))),
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0)))]
        #                        )
        #
        #             # bottom right triangle
        #             writeFaces(outfile, [
        #                 translate(
        #                     transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0))),
        #                 translate(transformTo3D(
        #                     transformToLatLong(projSrc, projDest, bottom_right_r, radius_bottom_right, T0)))]
        #                        )
        #     break


if __name__ == '__main__':
    start = time.time()

    with open("test.stl", "w+") as dest:
        with rasterio.open(GEOTIFF) as src:
            test02(src, dest).convert()

    end = time.time()
    print("Took: %g" % (end - start))
