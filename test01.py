import numpy as np
import rasterio
import math
import time
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

pixel2m = lambda v: v * GEOTIFF_SCALE + MOON_RADIUS
# https://en.wikipedia.org/wiki/Spherical_coordinate_system
# lonlat23d = lambda r, lo, la: (r * math.sin(la) * math.cos(lo), r * math.sin(la) * math.sin(lo), r * math.cos(lo))
PATCH_SIZE = 2000


def writeOLD(solid, file):
    name = solid.name
    if name is None:
        name = "unnamed"

    file.write(("solid %s\n" % name).encode())
    for facet in solid.facets:
        file.write(("  facet normal %g %g %g\n" % facet.normal).encode())
        file.write(b"    outer loop\n")
        for vertex in facet.vertices:
            file.write(("      vertex %g %g %g\n" % vertex).encode())
        file.write(b"    endloop\n")
        file.write(b"  endfacet\n")
    file.write(("endsolid %s\n" % name).encode())


def writeSolidStart(file, name="test"):
    file.write(("solid %s\n" % name).encode())


def tuple_3_diff(a, b):
    return (b[0] - a[0], b[1] - a[1], b[2] - a[2])


def writeFaces(file, vertices):
    v1 = tuple_3_diff(vertices[0], vertices[1])
    v2 = tuple_3_diff(vertices[0], vertices[2])

    # doesn't seem correct
    normal = (
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    )
    normal_length = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]
    normal_length = math.sqrt(normal_length)
    normal = (normal[0] / normal_length, normal[1] / normal_length, normal[2] / normal_length)

    file.write(("  facet normal %g %g %g\n" % normal).encode())
    file.write(b"    outer loop\n")
    for vertex in vertices:
        file.write(("      vertex %g %g %g\n" % vertex).encode())
    file.write(b"    endloop\n")
    file.write(b"  endfacet\n")


def writeSolidEnd(file, name="test"):
    file.write(("endsolid %s\n" % name).encode())


def transformToLatLong(projSrc, projDest, point, r, T0):
    p = (point[1], point[0]) * T0

    lo, la = transform(projSrc, projDest, p[0], p[1])
    la_rad = p[0] / 180.0 * math.pi
    lo_rad = p[1] / 180.0 * math.pi
    # ret = (r * math.sin(la_rad) * math.cos(lo_rad), r * math.sin(la_rad) * math.sin(lo_rad), r * math.cos(la_rad))

    # more z-stable ?
    ret = (r * math.cos(la_rad) * math.cos(lo_rad), r * math.cos(la_rad) * math.sin(lo_rad), r * math.sin(la_rad))
    # print lo, la, ret
    return ret


def transformTo3D(p):
    return p


def get_point_from(infile, point):
    return 0.0

    x = point[0]
    y = point[1]
    area = infile.read(1, window=((x, y), (x, y)))

    # why does infile.read sometimes return an empty area?
    if len(area.shape) == 0 or area.shape[0] == 0 or area.shape[1] == 0:
        return MOON_RADIUS

    return area[(0, 0)]


def convert(infile, outfile):
    print("Width %g, height %g, crs %s, transform %s" % (infile.width, infile.height, infile.crs, infile.transform))

    #
    # SETUP Transformation
    #
    T0 = src.affine  # upper-left pixel corner affine transform
    T0 = T0 * Affine.translation(0.5, 0.5)  # Get affine transform for pixel centres
    projSrc = Proj(src.crs)
    # Function to convert pixel row/column index (from 0) to easting/northing at centre
    # rc2en = lambda r, c: (c, r) * T0
    projDest = Proj(proj='latlong', datum='WGS84')
    pixel2m = lambda v: MOON_RADIUS + v * GEOTIFF_SCALE
    translate = lambda p: (p[0] + MAX_MOON_RADIUS, p[1] + MAX_MOON_RADIUS, p[2] + MAX_MOON_RADIUS)
    swap = lambda p: p #(p[1], p[0])

    writeSolidStart(outfile)

    min_height = MOON_RADIUS * 2
    max_height = -(MOON_RADIUS * 2)

    for col in range(0, infile.height - PATCH_SIZE, PATCH_SIZE):
        for row in range(0, infile.width - PATCH_SIZE, PATCH_SIZE):
            if col + PATCH_SIZE >= infile.height or row + PATCH_SIZE >= infile.width:
                continue

            # transformTo3D(transformToLatLong(projSrc, projDest, swap((row, col)), MOON_RADIUS, T0))
            # continue
            part_shape = (min(infile.width, row + PATCH_SIZE), min(infile.height, col + PATCH_SIZE))

            # _f file/window coordinates
            top_left_f = swap((0, 0))  # (row, col)
            top_right_f = swap((0, part_shape[1]))  # (row + PATCH_SIZE, col)
            bottom_left_f = swap((part_shape[0], 0))  # (row, col + PATCH_SIZE)
            bottom_right_f = swap((part_shape[0], part_shape[1]))  # (row + PATCH_SIZE, col + PATCH_SIZE)

            # _r real world PIXEL - coordinates
            top_left_r = swap((row, col))
            top_right_r = swap((row + part_shape[0], col))
            bottom_left_r = swap((row, col + part_shape[1]))
            bottom_right_r = swap((row + part_shape[0], col + part_shape[1]))

            radius_top_left = pixel2m(get_point_from(infile, top_left_r))
            radius_top_right = pixel2m(get_point_from(infile, top_right_r))
            radius_bottom_left = pixel2m(get_point_from(infile, bottom_left_r))
            radius_bottom_right = pixel2m(get_point_from(infile, bottom_right_r))

            min_height = min(min_height, radius_top_left, radius_bottom_left, radius_top_right, radius_bottom_right)
            max_height = max(max_height, radius_top_left, radius_bottom_left, radius_top_right, radius_bottom_right)

            # Write a patch:
            # +---------+
            # |        /|
            # |      /  |
            # |    /    |
            # |  /      |
            # |/        |
            # +---------+

            temp = transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))

            if False:  # clockwise
                # top left triangle
                writeFaces(outfile, [
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0)))]
                )

                # bottom right triangle
                writeFaces(outfile, [
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_right_r, radius_bottom_right, T0)))]
                )
            else:  # counter-clockwise
                # top left triangle
                writeFaces(outfile, [
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_left_r, radius_top_left, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0)))]
                )

                # bottom right triangle
                writeFaces(outfile, [
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, top_right_r, radius_top_right, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_left_r, radius_bottom_left, T0))),
                    translate(transformTo3D(transformToLatLong(projSrc, projDest, bottom_right_r, radius_bottom_right, T0)))]
                )
        break

    writeSolidEnd(outfile)
    print "max/min height", max_height, min_height


#     part = src.read(1, window=((0, 1000), (0, 1000)))
#     print pixel2m(part[(100, 100)])
#
#     # Function to convert pixel row/column index (from 0) to easting/northing at centre
#     rc2en = lambda r, c: (c, r) * T1
#
#     # All eastings and northings (there is probably a faster way to do this)
#     eastlings, northings = rc2en(100, 100)
#
#     # Project all longitudes, latitudes
#
#     print(long, lat)
#     print(lonlat23d(MOON_RADIUS, long, lat))


start = time.time()

with open("test.stl", "w+") as outfile:
    with rasterio.open(GEOTIFF) as src:
        convert(src, outfile)

end = time.time()
print("Took: %g" % (end - start))

# writeSolidStart(f)
# writeFaces(f, [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
# writeSolidEnd(f)
# with rasterio.open(GEOTIFF) as src:
#     print(src.width, src.height)
#     print(src.crs)
#     print(src.transform)
#     print(src.count)
#     print(src.indexes)
#
#     pixel2m = lambda v: v * GEOTIFF_SCALE + MOON_RADIUS
#
#     part = src.read(1, window=((0, 1000), (0, 1000)))
#     print pixel2m(part[(100, 100)])
#
#     T0 = src.affine  # upper-left pixel corner affine transform
#     p1 = Proj(src.crs)
#
#     # Get affine transform for pixel centres
#     T1 = T0 * Affine.translation(0.5, 0.5)
#     # Function to convert pixel row/column index (from 0) to easting/northing at centre
#     rc2en = lambda r, c: (c, r) * T1
#
#     # All eastings and northings (there is probably a faster way to do this)
#     eastlings, northings = rc2en(100, 100)
#
#     # Project all longitudes, latitudes
#     p2 = Proj(proj='latlong', datum='WGS84')
#     long, lat = transform(p1, p2, eastlings, northings)
#
#     print(long, lat)
#
#     # https://en.wikipedia.org/wiki/Spherical_coordinate_system
#     lonlat23d = lambda r, lo, la: (r * math.sin(la) * math.cos(lo), r * math.sin(la) * math.sin(lo), r * math.cos(lo))
#
#     print(lonlat23d(MOON_RADIUS, long, lat))
#
# # temp = rasterio.band(src, 1)
# # print src.affine
# # part = src.read(1, window=((0, 10), (0, 10)))
# # print part[(0, 0)]
# #  src.build_overviews([4096], Resampling.average)
