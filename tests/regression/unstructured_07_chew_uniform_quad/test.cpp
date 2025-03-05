#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/sources.hpp"
#include "../../../src/refine.hpp"
#include "../../../src/builder.hpp"

/*
  Simple test refining a quadralateral and the Chew uniform algorithm.
  Tests that density affects the initial split of source segments.
*/
int main(int argc, char* argv[]) {
    umr::mesh::Point p0(0.0, 0.0);
    umr::mesh::Point p1(0.0, 1.0);
    umr::mesh::Point p2(1.0, 1.0);
    umr::mesh::Point p3(1.0, 0.0);

    auto rect = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);
    rect->set_number_subsegments({1, 2, 3, 4});

    umr::source::SourceManager sources;
    sources.add_shape(rect);

    umr::density::DensityManager density;

    umr::Builder builder;
    builder.set_refinement_chew_uniform();
    builder.set_sources(&sources);
    builder.set_density(&density);

    density.set_minimum_density(1);
    umr::MeshGenerator gen0 = builder.build();
    const umr::mesh::MeshUnstructuredConstrained& m0 = gen0.triangulate();
    bool is_valid = validate_mesh(m0);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    density.set_minimum_density(0.2);
    umr::MeshGenerator gen1 = builder.build();
    const umr::mesh::MeshUnstructuredConstrained& m1 = gen1.triangulate();
    is_valid = validate_mesh(m1, 2);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    density.set_no_minimum_density();
    umr::MeshGenerator gen2 = builder.build();
    const umr::mesh::MeshUnstructuredConstrained& m2 = gen1.triangulate();
    is_valid = validate_mesh(m2, 3);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
