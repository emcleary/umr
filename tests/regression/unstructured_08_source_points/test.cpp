#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/sources.hpp"
#include "../../../src/refine.hpp"
#include "../../../src/builder.hpp"

/*
  Simple test refining a quadralateral with the Chew uniform
  algorithm. It validates that source points can be inserted,
  kept, and that they affect how segments get initially split.
*/

int main(int argc, char* argv[]) {
    umr::mesh::Point p0(0.0, 0.0);
    umr::mesh::Point p1(0.0, 1.0);
    umr::mesh::Point p2(1.0, 1.0);
    umr::mesh::Point p3(1.0, 0.0);
    umr::mesh::Point ps0(0.25, 0.5);
    umr::mesh::Point ps1(0.50, 0.5);

    double hmin = 1;

    auto rect = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);

    umr::source::SourceManager sources0;
    sources0.add_shape(rect);
    sources0.add_point(ps0);

    umr::source::SourceManager sources1;
    sources1.add_shape(rect);
    sources1.add_point(ps1);

    umr::density::DensityManager density;
    density.set_minimum_density(hmin);

    umr::Builder builder;
    builder.set_refinement_chew_uniform();
    builder.set_density(&density);

    builder.set_sources(&sources0);
    umr::MeshGenerator gen0 = builder.build();
    const umr::mesh::MeshUnstructuredConstrained& m0 = gen0.triangulate();
    bool is_valid = validate_mesh(m0);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    builder.set_sources(&sources1);
    umr::MeshGenerator gen1 = builder.build();
    const umr::mesh::MeshUnstructuredConstrained& m1 = gen1.triangulate();
    is_valid = validate_mesh(m1, 2);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
