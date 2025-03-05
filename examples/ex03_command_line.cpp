#include "../src/umr.hpp"

/*
 * Command line arguments are a very handy means of testing scripts.
 * The most notable features implemented here are setting the minimum
 * angle, a minimum density, and the algorithm.
 *
 * This example is hard coded to use the Chew Uniform algorithm. It
 * can be run with other algorithms from the command line, e.g.
 *
 * ./ex03_command_line --refinement ruppert
 * ./ex03_command_line --refinement chew-nonuniform
 *
 * Convergence criteria can be also be set, e.g.
 *
 * ./ex03_command_line --refinement ruppert --angle 15
 * ./ex03_command_line --refinement chew-nonuniform --density 0.03
 *
 * For a full list of valid arguments, run
 * 
 * ./ex03_command_line --help
 */

int main(int argc, char* argv[]) {

    /*
     * To allow for command line argument overriding, construct
     * the builder with the main function's input arguments.
     */

    umr::Builder builder(argc, argv);
    umr::source::SourceManager sources;
    umr::density::DensityManager density;
    builder.set_sources(&sources);
    builder.set_density(&density);

    /*
     * Add a source to the source manager.
     */

    umr::mesh::Point center(0, 0);
    double radius = 1;

    auto circle = std::make_shared<umr::source::ShapeCircle>(center, radius);
    sources.add_shape(circle);

    /*
     * Set default convergence criteria.
     */

    density.set_minimum_density(0.2);
    builder.set_min_angle(20); // degrees

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_chew_uniform();

    /*
     * The builder will override the above settings with command line
     * arguments, if any, while constructing the generator.
     */

    umr::MeshGenerator generator = builder.build();

    /*
     * Lastly triangulate and refine the mesh, and print out mesh
     * quality.
     */

    generator.triangulate();
    generator.dump("ex03_command_line");
    generator.quality();

    return 0;
}
