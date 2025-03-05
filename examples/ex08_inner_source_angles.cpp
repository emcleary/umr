#include "../src/builder.hpp"


/*
  This is another example of Shewchuck's Terminator algorithm.  This
  example demonstrates its use when interior segments form the narrow
  angles, rather than exterior segments.
*/

umr::mesh::Point make_point(umr::mesh::Point& center, double radius, double angle) {
    angle *= std::numbers::pi / 180;
    double x = center.x + radius * cos(angle);
    double y = center.y + radius * sin(angle);
    return umr::mesh::Point(x, y);
}

int main(int argc, char** argv) {

    umr::Builder builder(argc, argv);
    umr::source::SourceManager sources;
    umr::density::DensityManager density;
    builder.set_sources(&sources);
    builder.set_density(&density);

    /*
     * Set the exterior and interior segments
     */

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(1, 1);
    umr::mesh::Point p3(0, 1);

    umr::mesh::Point pc(0.5, 0.5);
    std::vector<umr::mesh::Point> points = {
        make_point(pc, 0.2, 0),
        make_point(pc, 0.2, 15),
        make_point(pc, 0.2, 120),
        make_point(pc, 0.15, 140),
        make_point(pc, 0.25, 175),
        make_point(pc, 0.2, 185),
        make_point(pc, 0.2, 195),
        make_point(pc, 0.3, 265),
        make_point(pc, 0.1, 275),
        make_point(pc, 0.3, 290),
    };

    auto square = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);
    sources.add_shape(square);

    for (auto point : points) {
        auto line = std::make_shared<umr::parametric::Line>(pc, point);
        sources.add_segment(line);
    }

    /*
     * Set default convergence criteria.
     */

    density.set_no_minimum_density();
    builder.set_min_angle(25);

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_ruppert();

    /*
     * Build the mesh generator.
     */

    umr::MeshGenerator generator = builder.build();


    /*
     * Triangulate and refine the meshes.
     */

    generator.triangulate();
    generator.dump("ex08_inner_source_angles.vtk");

    /*
     * Display mesh quality. 
     */

    generator.quality();

    return 0;
}
