#include <iostream>
#include <fstream>
#include <json/json.h>

#include "../src/common.hpp"
#include "../src/validate.hpp"
#include "../../../src/sources.hpp"
#include "../../../src/builder.hpp"

/*
  This test checks that the Terminator algorithms prevents convergence
  issues due to small input geometry angles. Both the Ruppert and
  Chew-Nonuniform refinements algorithms are tested with and without
  the Terminator algorithm.
*/


static const int MAX_ITERATIONS = 25;


umr::mesh::Point make_point(double r, double theta) {
    double x = r * cos(theta * std::numbers::pi / 180);
    double y = r * sin(theta * std::numbers::pi / 180);
    return umr::mesh::Point(x, y);
}


umr::MeshGenerator run_test(umr::Builder& builder) {
    umr::MeshGenerator generator = builder.build();
    generator.triangulate(MAX_ITERATIONS);
    return generator;
}


int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    builder.set_min_angle(20);

    umr::density::DensityManager density;
    density.set_no_minimum_density();
    builder.set_density(&density);

    umr::source::SourceManager sources;
    builder.set_sources(&sources);

    // External boundary of mesh is a triangle
    // with angles of 10, 28.9, and 141.1
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1 = make_point(1.3, 0);
    umr::mesh::Point p2 = make_point(1, 10);

    auto s01 = std::make_shared<umr::parametric::Line>(p0, p1);
    auto s12 = std::make_shared<umr::parametric::Line>(p1, p2);
    auto s20 = std::make_shared<umr::parametric::Line>(p2, p0);

    sources.add_segment(s01);
    sources.add_segment(s12);
    sources.add_segment(s20);

    builder.set_refinement_chew_nonuniform();
    builder.set_terminator_on();
    umr::MeshGenerator g0 = run_test(builder);
    if (!g0.is_converged()) {
        std::cout << "Chew nonuniform with Terminator did not converge\n";
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    bool is_valid = validate_mesh(g0.finalize());
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    builder.set_refinement_chew_nonuniform();
    builder.set_terminator_off();
    umr::MeshGenerator g1 = run_test(builder);
    if (g1.is_converged()) {
        std::cout << "Chew nonuniform without Terminator did converge\n";
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }
    
    builder.set_refinement_ruppert();
    builder.set_terminator_on();
    umr::MeshGenerator g2 = run_test(builder);
    if (!g2.is_converged()) {
        std::cout << "Ruppert with Terminator did not converge\n";
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    is_valid = validate_mesh(g2.finalize(), 1);
    if (!is_valid) {
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    builder.set_refinement_ruppert();
    builder.set_terminator_off();
    umr::MeshGenerator g3 = run_test(builder);
    if (g3.is_converged()) {
        std::cout << "Ruppert without Terminator did converge\n";
        std::cout << "Aborting\n";
        return EXIT::FAIL;
    }

    return EXIT::SUCCESS;
}
