#include "../src/umr.hpp"


/*
 * The UMR framework is centered around using a Builder. Once
 * parameters are set, the builder can be used to build a mesh
 * generator, which is then run to triangulate the mesh and write out
 * the results.
 *
 * This example demonstrates Delaunay triangulation of a set of
 * points. It runs the Divide and Conquer algorithm as implemented by
 * Guibas and Stolfi. See "Guibas, Leonidas, and Jorge
 * Stolfi. "Primitives for the manipulation of general subdivisions
 * and the computation of Voronoi." ACM transactions on graphics (TOG)
 * 4.2 (1985): 74-123."
 */
int main(int argc, char* argv[]) {

    // Define a set of points to use in the triangulation.
    std::vector<std::pair<double, double>> points = {
        {0.0, 0.5}, {0.1, 0.9}, {0.2, 0.1}, {0.3, 0.6},
        {0.4, 0.2}, {0.5, 1.0}, {0.6, 0.3}, {0.7, 0.8},
        {0.8, 0.4}, {0.9, 0.7}, {1.0, 0.0},
    };

    // Input points are set as source points. They must be added to a
    // SourceManager.
    umr::source::SourceManager sources;
    for (auto [x, y] : points)
        sources.add_point(x, y);

    // Add the sources to the Builder.
    umr::Builder builder;
    builder.set_sources(&sources);

    // Build the mesh generator, run the Divide and Conquer algorithm,
    // then dump the mesh to a VTK file.
    umr::MeshGenerator generator = builder.build();
    generator.triangulate();
    generator.dump("ex01_delaunay_triangulation");

    return 0;
}
