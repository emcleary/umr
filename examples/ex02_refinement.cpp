#include "../src/umr.hpp"

/*
 * The UMR framework has three algorithms implemented for mesh
 * refinement:
 *
 *   1) Chew Uniform
 *   2) Chew Nonuniform
 *   3) Ruppert
 *
 * The objective of these algorithms is to refine a mesh with the
 * purpose of increasing the minimum angle of its triangles.
 *
 * The Chew Uniform algorithm is guaranteed to have a minimum angle of
 * 30 degrees, but results in a very high density mesh.
 *
 * The Chew Nonuniform algorithm returns a coarse mesh and is
 * guaranteed to converge for input minimum angles of up to 26.5
 * degrees. The major con of this algorithm is that it is not strictly
 * meet the Delaunay criteria.
 *
 * The Ruppert algorithm is similar to the Chew Nonuniform algorithm
 * in that it returns a coarser mesh than the Chew Uniform algorithm.
 * Its result is guaranteed to meet the Delaunay criteria, but is only
 * guaranteed to converge for minimum angles up to 20.7 degrees.
 *
 * More often than not, the Chew Nonuniform and Ruppert algorithms
 * will converge for minimum angles up to 30 degrees. Just be aware
 * that there is no guarantee.
 *
 * For references, I recommend the following papers:
 *
 * Chew, L. Paul. Guaranteed-quality triangular meshes. Cornell University, 1989.
 *
 * Chew, L. Paul. "Guaranteed-quality mesh generation for curved
 * surfaces." Proceedings of the ninth annual symposium on
 * Computational geometry. 1993.
 *
 * Ruppert, Jim. "A new and simple algorithm for quality 2-dimensional
 * mesh generation." Proceedings of the fourth annual ACM-SIAM
 * Symposium on Discrete algorithms. 1993.
 *
 * Shewchuk, J. R. "Mesh generation for domains with small angles."
 * Proceedings of the sixteenth annual Symposium on Computational
 * Geometry. 2000.
 */

int main(int argc, char* argv[]) {

    umr::Builder builder;
    umr::source::SourceManager sources;
    umr::density::DensityManager density;
    builder.set_sources(&sources);
    builder.set_density(&density);

    /*
     * Refined meshes require source edges, or segments, to define a
     * closed shape. Below is an example of a circle shape. More
     * shapes are available and are demonstrated in other examples.
     *
     * As with ex01, these sources must be added to the source manager.
     */

    umr::mesh::Point center(0, 0);
    double radius = 1;

    auto circle = std::make_shared<umr::source::ShapeCircle>(center, radius);
    sources.add_shape(circle);

    /*
     * There are 2 main parameters implemented for convergence with
     * these algorithms: the minimum angle, and the density. Both are
     * set here.
     */

    density.set_minimum_density(0.2);
    builder.set_min_angle(20); // degrees

    /*
     * Now the builder is ready to build some mesh generators.
     * Specify the algorithm to be used, then run the builder.
     */

    builder.set_refinement_chew_uniform();
    umr::MeshGenerator generator_chew_uniform = builder.build();

    builder.set_refinement_chew_nonuniform();
    umr::MeshGenerator generator_chew_nonuniform = builder.build();

    builder.set_refinement_ruppert();
    umr::MeshGenerator generator_ruppert = builder.build();

    /*
     * Finally, each mesh can be trianugulated and refined.
     */

    generator_chew_uniform.triangulate();
    generator_chew_uniform.dump("ex02_refinement_chew_uniform");

    generator_chew_nonuniform.triangulate();
    generator_chew_nonuniform.dump("ex02_refinement_chew_nonuniform");

    generator_ruppert.triangulate();
    generator_ruppert.dump("ex02_refinement_ruppert");

    /*
     * Lastly, the mesh generator can be used to display quality of
     * the meshes.
     */

    std::cout << "\nChew Uniform algorithm\n\n";
    generator_chew_uniform.quality();

    std::cout << "\nChew Nonuniform algorithm\n\n";
    generator_chew_nonuniform.quality();

    std::cout << "\nRuppert algorithm\n\n";
    generator_ruppert.quality();

    return 0;
}
