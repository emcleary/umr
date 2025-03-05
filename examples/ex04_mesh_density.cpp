#include "../src/builder.hpp"


/*
 * Density is an important criterion in the refinement algorithms, as
 * it specifies the largest triangle allowed at each point in the
 * mesh's domain. There are 3 types of density profiles than can be
 * used:
 *
 * 1) A fixed density
 * 2) A density function
 * 3) An interpolator
 *
 * Any number of these can be added to a single DensityManager.  It
 * will evaluate each and return the minimum value. Errors will be
 * thrown if the minimum value is negative. To avoid these errors, it
 * is recommended to set a fixed minimum density which will override
 * smaller positive and all negative values.
 *
 * If no densites are set, then refinement algorithms will only use
 * the minimum angle criterion for convergence.
 */

umr::density::DensityManager make_density_uniform() {
    umr::density::DensityManager density;
    density.set_minimum_density(0.1);
    return density;
}

umr::density::DensityManager make_density_interpolator() {
    umr::density::DensityManager density;
    density.add_interpolator("ex04.csv");
    return density;
}

umr::density::DensityManager make_density_function() {
    umr::density::DensityFunction function = [](double x, double y) {
        return 0.2 + 0.12 * sin((x - y) * std::numbers::pi * 3 / 4);
    };
    umr::density::DensityManager density;
    density.add_functor(function);
    return density;
}


int main(int argc, char* argv[]) {

    umr::Builder builder(argc, argv);
    umr::source::SourceManager sources;
    builder.set_sources(&sources);

    /*
     * Set the external boundary of the mesh
     */

    double xmin = 0;
    double xmax = 4;
    double yavg = 2;
    double amp = 0.1;
    umr::mesh::Point p0(xmin, yavg);
    umr::mesh::Point p1(xmin, 0);
    umr::mesh::Point p2(xmax, 0);
    umr::mesh::Point p3(xmax, yavg);

    double tstart = xmin;
    double tend = xmax;
    umr::parametric::ParametricFunction f = [tend, yavg, amp](double t) -> std::pair<double, double> {
        double x = tend - t;
        double y = yavg + amp * sin(x * 2 * std::numbers::pi);
        return std::make_pair(x, y);
    };

    auto exterior = std::make_shared<umr::source::Shape>();

    exterior->add_line(p0, p1);
    exterior->add_line(p1, p2);
    exterior->add_line(p2, p3);

    int n_param_subseg = 24;
    exterior->add_parametric(f, tstart, tend, n_param_subseg);

    sources.add_shape(exterior);

    /*
     * Set default convergence criteria.
     */

    umr::density::DensityManager density_uniform = make_density_uniform();
    umr::density::DensityManager density_interpolator = make_density_interpolator();
    umr::density::DensityManager density_function = make_density_function();
    builder.set_min_angle(25); // degrees

    /*
     * Specify a default refinement algorithm.
     */

    builder.set_refinement_ruppert();

    /*
     * Build the mesh generators.
     */

    builder.set_density(&density_uniform);
    umr::MeshGenerator generator_uniform = builder.build();

    builder.set_density(&density_interpolator);
    umr::MeshGenerator generator_interpolator = builder.build();

    builder.set_density(&density_function);
    umr::MeshGenerator generator_function = builder.build();

    /*
     * Triangulate and refine the meshes.
     */
    
    generator_uniform.triangulate();
    generator_uniform.dump("ex04_density_uniform");
    generator_interpolator.triangulate();
    generator_interpolator.dump("ex04_density_interpolate");
    generator_function.triangulate();
    generator_function.dump("ex04_density_function");

    /*
     * Display mesh quality. 
     */

    std::cout << "\nUniform density\n\n";
    generator_uniform.quality();

    std::cout << "\nInterpolated density\n\n";
    generator_interpolator.quality();

    std::cout << "\nFunction density\n\n";
    generator_function.quality();

    return 0;
}
