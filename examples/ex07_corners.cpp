#include "../src/builder.hpp"
#include "../src/command_line_interface.hpp"


/*
 * Shewchuck's Terminator algorithm was implemented in both Ruppert's
 * algorithm and Chew's Nonuniform algorithm to tackle convergence
 * issues in geometries with small input segment angles. This example
 * demonstrates the use of this algorithm.
 *
 * Running this program with the Terminator algorithm turned off will
 * not converge, but running it like that is a good demonstration of
 * the convergence issue. This can be done with the following, which
 * dumps mesh data at each iteration.
 *
 * ./ex07_corners --no-terminator --data-frequency 1
 */

umr::mesh::Point make_point(double r, double theta) {
    double x = r * cos(theta * std::numbers::pi / 180);
    double y = r * sin(theta * std::numbers::pi / 180);
    return umr::mesh::Point(x, y);
}

int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    umr::source::SourceManager sources;
    umr::density::DensityManager density;
    builder.set_sources(&sources);
    builder.set_density(&density);

    /*
     * Set the exterior and interior segments
     */

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1 = make_point(1.5, 0);
    umr::mesh::Point p2 = make_point(2.4, 10);
    umr::mesh::Point p3 = make_point(1.0, 20);
    umr::mesh::Point p4 = make_point(3.0, 30);
    umr::mesh::Point p5 = make_point(1.5, 40);

    auto s01 = std::make_shared<umr::parametric::Line>(p0, p1);
    auto s02 = std::make_shared<umr::parametric::Line>(p0, p2);
    auto s03 = std::make_shared<umr::parametric::Line>(p0, p3);
    auto s12 = std::make_shared<umr::parametric::Line>(p1, p2);
    auto s23 = std::make_shared<umr::parametric::Line>(p2, p3);

    sources.add_segment(s01);
    sources.add_segment(s02);
    sources.add_segment(s03);
    sources.add_segment(s12);
    sources.add_segment(s23);

    auto s04 = std::make_shared<umr::parametric::Line>(p0, p4);
    auto s34 = std::make_shared<umr::parametric::Line>(p3, p4);
    sources.add_segment(s04);
    sources.add_segment(s34);

    auto s05 = std::make_shared<umr::parametric::Line>(p0, p5);
    auto s45 = std::make_shared<umr::parametric::Line>(p4, p5);
    sources.add_segment(s05);
    sources.add_segment(s45);

    /*
     * Set default convergence criteria.
     */

    builder.set_min_angle(10); // degrees

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_chew_nonuniform();

    /*
     * Build the mesh generator.
     */

    umr::MeshGenerator generator = builder.build();

    /*
     * Triangulate and refine the meshes.
     */

    generator.triangulate();
    generator.dump("ex07_corners");

    /*
     * Display mesh quality. 
     */

    generator.quality();

    return 0;
}
