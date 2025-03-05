#include <cassert>
#include <memory>

#include "common.hpp"
#include "../src/integrators.hpp"
#include "../src/parametrics.hpp"
#include "../src/inequalities.hpp"


BOOST_AUTO_TEST_SUITE(test_integrators)


BOOST_AUTO_TEST_CASE(test_integrate_line) {
    const double tol = 1e-8;
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(0, 1);
    umr::mesh::Point p2(1, 1);
    umr::parametric::Line line01(p0, p1);
    umr::parametric::Line line02(p0, p2);

    umr::math::integrate::SegmentLength integrator;
    
    double d = integrator.run(&line01);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = integrator.run(&line02);
    BOOST_CHECK_CLOSE_FRACTION(d, sqrt(2), tol);
}


BOOST_AUTO_TEST_CASE(test_integrate_arc) {
    const double tol = 1e-8;

    umr::mesh::Point pc(0, 0);
    double radius0 = 1;
    double radius1 = 3;

    umr::parametric::Arc arc0(pc, radius0, 0, std::numbers::pi);
    umr::parametric::Arc arc1(pc, radius1, 0, 2 * std::numbers::pi);

    umr::math::integrate::SegmentLength integrator;
    
    double d = integrator.run(&arc0);
    BOOST_CHECK_CLOSE_FRACTION(d, radius0 * std::numbers::pi, tol);

    d = integrator.run(&arc1);
    BOOST_CHECK_CLOSE_FRACTION(d, 2 * radius1 * std::numbers::pi, tol);
}


BOOST_AUTO_TEST_CASE(test_integrate_parametrics) {
    const double tol = 1e-4;

    auto circle = [](double t) -> std::pair<double, double> {
        double x = t;
        double y = sqrt(1 - t * t);
        return {x, y};
    };

    auto differentiate = [](double t) -> std::pair<double, double> {
        bool at_singularity = (umr::inequalities::is_close(t, 1, 1e-12, 0)
                || umr::inequalities::is_close(t, -1, 1e-12, 0));
        double x = 1;
        double y = at_singularity ? 0 : - t / sqrt(1 - t * t);
        return {x, y};
    };

    umr::parametric::Parametric param(circle, -1, 1);
    umr::parametric::Parametric param_deriv(circle, differentiate, -1, 1);

    umr::math::integrate::SegmentLength integrator;
    double d = integrator.run(&param);
    BOOST_CHECK_CLOSE_FRACTION(d, std::numbers::pi, tol);

    d = integrator.run(&param_deriv);
    BOOST_CHECK_CLOSE_FRACTION(d, std::numbers::pi, tol);
}


BOOST_AUTO_TEST_SUITE_END()
