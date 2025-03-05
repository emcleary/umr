#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/refine.hpp"


/*
 * This tests the giftwrapping algorithm for inserting edges shorter
 * than and longer than neighboring cells. The initial mesh is
 * generated with a uniform mesh.
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

    // Inserting short edge
    umr::mesh::Point& pa = mesh.make_point(0.5, 0.5);
    umr::mesh::Point& pb = mesh.make_point(pa.x, pa.y+hmin);
    umr::mesh::Edge* e = &umr::mesh::quadedge::make_edge(pa, pb);
    umr::mesh::refine::insert_edge(mesh, *e, hmin);

    is_valid = validate_mesh(mesh, 2);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    // Inserting long edge
    umr::mesh::Point& pc = mesh.make_point(0.3, 0.0);
    umr::mesh::Point& pd = mesh.make_point(0.7, 0.0);
    e = &umr::mesh::quadedge::make_edge(pc, pd);
    umr::mesh::refine::insert_edge(mesh, *e, hmin);

    is_valid = validate_mesh(mesh, 3);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
