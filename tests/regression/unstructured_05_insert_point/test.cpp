#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/delaunay.hpp"
#include "../../../src/refine.hpp"


/*
 * This tests point insertion and triangulation afterwards.
 */

int main(int argc, char* argv[]) {
    double hmin = 0.1;
    umr::mesh::MeshUnstructuredConstrained mesh;
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(1, 1);
    umr::mesh::Point p3(0, 1);
    umr::mesh::refine::uniform_grid(mesh, p0, p2, hmin, 3);

    bool is_valid = validate_mesh(mesh);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    umr::mesh::Point* pa = new umr::mesh::Point(0.39, 0.5);
    umr::mesh::Edge* e = mesh.get_initial_edge();
    e = &umr::mesh::delaunay::locate(*pa, *e);
    e = &umr::mesh::delaunay::make_cavity(mesh, *pa, *e);

    is_valid = validate_mesh(mesh, 2);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    std::vector<umr::mesh::Point*> new_points;
    new_points.push_back(pa);
    umr::mesh::delaunay::split_cavity(mesh, *e, new_points);

    is_valid = validate_mesh(mesh, 3);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    umr::mesh::Point* pb = new umr::mesh::Point(0.75, 0.5);
    e = umr::mesh::refine::insert_point(mesh, *pb, *e, hmin);

    is_valid = validate_mesh(mesh, 4);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    e = umr::mesh::refine::refine_triangles(mesh, *e, hmin);

    is_valid = validate_mesh(mesh, 5);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
