#include "../src/builder.hpp"


/*
 * As with density, sources can be used to specify mesh
 * resolution. Thi can be done through external segment lengths,
 * interior segments lengths, a set of interior source points, or a
 * distance between any combination of them.
 */

umr::source::SourceManager create_sources_with_points(umr::source::SourceManager& exterior) {
    umr::source::SourceManager manager(exterior);

    double x0 = 0.3;
    double y0 = 0.3;
    double hx = 0.03;
    double hy = 0.03;
    int n = 5;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            manager.add_point(x0 + i * hx, y0 + j * hy);
    
    return manager;
}


umr::source::SourceManager create_sources_with_segments(umr::source::SourceManager& exterior) {
    umr::source::SourceManager manager(exterior);

    using namespace umr::mesh;

    auto add_segments = [&](Point&& p0, Point&& p1, double dx, double dy) {
        Point p2(p0.x + dx, p0.y + dy);
        Point p3(p1.x + dx, p1.y + dy);
        auto line01 = std::make_shared<umr::parametric::Line>(p0, p1);
        auto line23 = std::make_shared<umr::parametric::Line>(p2, p3);
        manager.add_segment(line01);
        manager.add_segment(line23);
    };

    add_segments(Point(0.2, 0.2), Point(0.2, 0.4), 0.02, 0);
    add_segments(Point(0.5, 0.5), Point(0.5, 0.8), 0.04, 0);
    add_segments(Point(0.5, 0.1), Point(0.8, 0.1), 0, 0.03);

    return manager;
}


int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    umr::density::DensityManager density;
    umr::source::SourceManager exterior;
    builder.set_density(&density);

    /*
     * Set the exterior source segments
     */

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(1, 1);
    umr::mesh::Point p3(0, 1);

    auto square = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);
    exterior.add_shape(square);

    /*
     * Create managers with various interior sources
     */

    umr::source::SourceManager sources_with_points = create_sources_with_points(exterior);
    umr::source::SourceManager sources_with_segments = create_sources_with_segments(exterior);

    /*
     * Set default convergence criteria.
     */

    builder.set_min_angle(25); // degrees

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_ruppert();

    /*
     * Build the mesh generators.
     */

    builder.set_sources(&sources_with_points);
    umr::MeshGenerator generator_points = builder.build();

    builder.set_sources(&sources_with_segments);
    umr::MeshGenerator generator_segments = builder.build();

    /*
     * Triangulate and refine the meshes.
     */
    
    generator_points.triangulate();
    generator_points.dump("ex05_mesh_source_points");

    generator_segments.triangulate();
    generator_segments.dump("ex05_mesh_source_segments");

    /*
     * Display mesh quality. 
     */

    std::cout << "\nInterior source points\n\n";
    generator_points.quality();

    std::cout << "\nInterior source segments\n\n";
    generator_segments.quality();

    return 0;
}
