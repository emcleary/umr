#include "../src/builder.hpp"


/*
 * The library contains tools for creating shapes of source segments.
 * They include Quadrilateral, Triangle, and Circle, along with Shape
 * which can be constructed with any combination of IParametric
 * objects: Line, Arc, or Parametric (functions).
 *
 * Each shape comes with a method for setting the initial resolution
 * by specifying a number of subsegments of each segment. See the
 * documentation and the examples below. For any nonlinear functions
 * used to specify a segment in an shape, it is important that it
 * starts with a good set of subsegments to resolve any of its
 * features.
 */

umr::source::SourceManager create_circular_shapes() {

    umr::mesh::Point p0(0, 0);
    double r0 = 2;
    auto exterior = std::make_shared<umr::source::ShapeCircle>(p0, r0);
    exterior->set_number_subsegments({10});

    umr::mesh::Point p1(1.5, 0);
    double r1 = 0.3;
    auto interior = std::make_shared<umr::source::ShapeCircle>(p1, r1);
    interior->set_number_subsegments({10});

    umr::source::SourceManager manager;
    manager.add_shape(exterior);
    manager.add_shape(interior);

    return manager;
}


umr::source::SourceManager create_quad_shapes() {

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(3, 0);
    umr::mesh::Point p2(2, 1);
    umr::mesh::Point p3(1, 1);
    auto exterior = std::make_shared<umr::source::ShapeQuadrilateral>(p0, p1, p2, p3);

    umr::mesh::Point p4(1.3, 0.1);
    umr::mesh::Point p5(1.7, 0.1);
    umr::mesh::Point p6(1.7, 0.5);
    umr::mesh::Point p7(1.3, 0.5);
    auto interior = std::make_shared<umr::source::ShapeQuadrilateral>(p4, p5, p6, p7);
    interior->set_number_subsegments({5, 5, 5, 5});

    umr::source::SourceManager manager;
    manager.add_shape(exterior);
    manager.add_shape(interior);

    return manager;
}


umr::source::SourceManager create_triangle_shapes() {

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(3, 0);
    umr::mesh::Point p2(0, 3);
    auto exterior = std::make_shared<umr::source::ShapeTriangle>(p0, p1, p2);

    double x = 0.05;
    double y = 0.4;
    double h = 0.5;
    umr::mesh::Point p3(x, y);
    umr::mesh::Point p4(x + h * sqrt(3) / 2, y + h / 2);
    umr::mesh::Point p5(x, y + h);
    auto interior = std::make_shared<umr::source::ShapeTriangle>(p3, p4, p5);

    umr::source::SourceManager manager;
    manager.add_shape(exterior);
    manager.add_shape(interior);

    return manager;
}


umr::source::SourceManager create_arbitrary_shapes() {

    auto exterior = std::make_shared<umr::source::Shape>();

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    int n_line_subseg = 3;
    exterior->add_line(p0, p1, n_line_subseg);

    umr::mesh::Point center(1, 1);
    double radius = 1;
    double t_start_arc = - std::numbers::pi / 2;
    double t_end_arc = std::numbers::pi / 2;
    int n_arc_subseg = 5;
    exterior->add_arc(center, radius, t_start_arc, t_end_arc, n_arc_subseg);

    double t_start = 0;
    double t_end = 1;
    umr::parametric::ParametricFunction f = [t_end](double t) -> std::pair<double, double> {
        double x = t_end - t;
        double y = 2 + 0.1 * sin(x * 2 * std::numbers::pi);
        return std::make_pair(x, y);
    };
    int n_param_subseg = 5;
    exterior->add_parametric(f, t_start, t_end, n_param_subseg);

    umr::mesh::Point p2(f(t_end));
    exterior->add_line(p2, p0);

    umr::source::SourceManager manager;
    manager.add_shape(exterior);

    return manager;
}


int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    umr::density::DensityManager density;
    umr::source::SourceManager exterior;
    builder.set_density(&density);

    /*
     * Create managers with various sources
     */

    umr::source::SourceManager sources_circles = create_circular_shapes();
    umr::source::SourceManager sources_quads = create_quad_shapes();
    umr::source::SourceManager sources_triangles = create_triangle_shapes();
    umr::source::SourceManager sources_arbitrary = create_arbitrary_shapes();

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

    builder.set_sources(&sources_circles);
    umr::MeshGenerator generator_circles = builder.build();

    builder.set_sources(&sources_quads);
    umr::MeshGenerator generator_quads = builder.build();

    builder.set_sources(&sources_triangles);
    umr::MeshGenerator generator_triangles = builder.build();

    builder.set_sources(&sources_arbitrary);
    umr::MeshGenerator generator_arbitrary = builder.build();

    /*
     * Triangulate and refine the meshes.
     */
    
    generator_circles.triangulate();
    generator_circles.dump("ex06_shape_circles");

    generator_quads.triangulate();
    generator_quads.dump("ex06_shape_quads");

    generator_triangles.triangulate();
    generator_triangles.dump("ex06_shape_triangles");

    generator_arbitrary.triangulate();
    generator_arbitrary.dump("ex06_shape_arbitrary");

    /*
     * Display mesh quality. 
     */

    std::cout << "\nCircular shapes\n\n";
    generator_circles.quality();

    std::cout << "\nQuadrilateral shapes\n\n";
    generator_quads.quality();

    std::cout << "\nTriangular shapes\n\n";
    generator_triangles.quality();

    std::cout << "\nArbitrary shape\n\n";
    generator_arbitrary.quality();

    return 0;
}
