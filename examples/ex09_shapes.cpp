#include "../src/builder.hpp"


/*
 * As with example 6, this is another example with shapes.
 */

int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    umr::density::DensityManager density;
    umr::source::SourceManager sources;
    builder.set_density(&density);
    builder.set_sources(&sources);

    /*
     * Set the exterior and interior segments
     */

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(3, 0);
    umr::mesh::Point p2(2, 1);
    umr::mesh::Point p3(0, 1);
    auto quad = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);
    sources.add_shape(quad);

    umr::mesh::Point p4(0.1, 0.1);
    umr::mesh::Point p5(0.9, 0.1);
    umr::mesh::Point p6(0.9, 0.9);
    umr::mesh::Point p7(0.1, 0.9);
    auto square = std::make_shared<umr::source::ShapeQuadrilateral>(p4, p5, p6, p7);
    sources.add_shape(square);

    umr::mesh::Point p8(2.00, 0.2);
    umr::mesh::Point p9(2.60, 0.2);
    umr::mesh::Point p10(2.00, 0.8);
    auto triangle = std::make_shared<umr::source::ShapeTriangle>(p8, p9, p10);
    sources.add_shape(triangle);

    umr::mesh::Point center(1.5, 0.5);
    double radius = 0.4;
    auto circle = std::make_shared<umr::source::ShapeCircle>(center, radius);
    sources.add_shape(circle);

    /*
     * Set default convergence criteria.
     */

    density.set_minimum_density(0.3);
    builder.set_min_angle(25); // degrees

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_ruppert();

    /*
     * Build the mesh generators.
     */

    umr::MeshGenerator generator = builder.build();

    /*
     * Triangulate and refine the meshes.
     */
    
    generator.triangulate();
    generator.dump("ex09_shapes");

    /*
     * Display mesh quality. 
     */

    generator.quality();

    return 0;
}
