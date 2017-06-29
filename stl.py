import math

# stl2bin test.stl test_bin.stl


def tuple_3_diff(a, b):
    return b[0] - a[0], b[1] - a[1], b[2] - a[2]


def write_solid_end(f, name="test"):
    f.write(("endsolid %s\n" % name).encode())


def write_solid_start(f, name="test"):
    f.write(("solid %s\n" % name).encode())


def write_3_faces(f, vertices):
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

    f.write(("  facet normal %g %g %g\n" % normal).encode())
    f.write(b"    outer loop\n")
    for vertex in vertices:
        f.write(("      vertex %g %g %g\n" % vertex).encode())
    f.write(b"    endloop\n")
    f.write(b"  endfacet\n")
