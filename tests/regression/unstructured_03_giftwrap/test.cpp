#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/delaunay.hpp"

/*
 * This tests the giftwrapping algorithm for inserting edges.  The
 * mesh used is a square mesh with segments split into subsegments and
 * trinagulated with the Bowyer-Watson algorithm. No refinement is
 * done. Inserting the edge merely flips edges connecting segment
 * points back and forth.
 */

int main(int argc, char* argv[]) {
    umr::mesh::MeshUnstructured m;
    const size_t n = 10;
    for (size_t i = 0; i < n; i++) {
        m.make_point(i, 0);
        m.make_point(n, i);
        m.make_point(n-i, n);
        m.make_point(0, n-i);
    }
    umr::mesh::delaunay::bowyer_watson(m);

    bool is_valid = validate_mesh(m);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }
    
    umr::mesh::Point& p0 = m.make_point(0, 0);
    umr::mesh::Point& p1 = m.make_point(n, n);

    umr::mesh::Edge& e0 = *m.get_initial_edge();
    umr::mesh::Edge& e1 = umr::mesh::quadedge::make_edge(p0, p1);
    umr::mesh::delaunay::insert_edge(m, e1, e0);

    is_valid = validate_mesh(m, 2);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }
    
    umr::mesh::Point& p2 = m.make_point(0, n);
    umr::mesh::Point& p3 = m.make_point(n, 0);
    umr::mesh::Edge& e2 = umr::mesh::quadedge::make_edge(p2, p3);
    umr::mesh::delaunay::insert_edge(m, e2, e0);

    m.remove_external_edges();

    is_valid = validate_mesh(m, 3);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }
    

    return EXIT::SUCCESS;
}
