#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/refine.hpp"


/*
 * This tests generation of a uniform mesh given points
 * representing the lower and upper bounds.
 */

int main(int argc, char* argv[]) {
    umr::mesh::MeshUnstructuredConstrained mesh;
    umr::mesh::Point pmin(0, 0);
    umr::mesh::Point pmax(1, 1);
    umr::mesh::refine::uniform_grid(mesh, pmin, pmax, 0.1, 3);

    bool is_valid = validate_mesh(mesh);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
