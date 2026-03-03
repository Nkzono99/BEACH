from beach import make_box, make_cylinder, make_plane, make_sphere


def main() -> None:
    templates = {
        "plane": make_plane(size_x=1.0, size_y=1.0, nx=8, ny=8),
        "box": make_box(size=(1.0, 0.8, 0.6), nx=4, ny=4, nz=3),
        "cylinder": make_cylinder(radius=0.4, height=1.2, n_theta=32, n_z=8, cap=True),
        "sphere": make_sphere(radius=0.5, n_lon=32, n_lat=16),
    }

    for name, mesh in templates.items():
        print(f"{name:10s} triangles={mesh.nelem:6d} area_sum={mesh.areas.sum():.6e}")


if __name__ == "__main__":
    main()
