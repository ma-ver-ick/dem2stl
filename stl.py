import math


def tuple_3_diff(a, b):
    return b[0] - a[0], b[1] - a[1], b[2] - a[2]


def write_solid_end(file, name="test"):
    file.write(("endsolid %s\n" % name).encode())


def write_solid_start(file, name="test"):
    file.write(("solid %s\n" % name).encode())


def write_3_faces(file, vertices):
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
