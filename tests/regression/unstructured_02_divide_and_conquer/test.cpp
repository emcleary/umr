#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/delaunay.hpp"


/*
 * A simple test of the Divide and Conquer algorithm for triangulating
 * a set of points.
 */

int main(int argc, char* argv[]) {
    std::vector<std::vector<double>> point_data = {
        {0, 0},
        {5, 0},
        {5, 6},
        {0, 5},
        {2, 2},
        {0, 3},
        {4, 4},
        {3, 4},
        {1, 2},
        {4, 1},
    };

    umr::mesh::MeshUnstructured m;
    for (std::vector<double>& p : point_data)
        m.make_point(p[0], p[1]);

    umr::mesh::delaunay::divide_and_conquer(m);

    bool is_valid = validate_mesh(m);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }
    
    return EXIT::SUCCESS;
}
