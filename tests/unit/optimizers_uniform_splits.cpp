#include <cassert>

#include "common.hpp"
#include "../src/optimizers.hpp"
#include "../src/point.hpp"


BOOST_AUTO_TEST_SUITE(test_optimizers_uniform_splits)


BOOST_AUTO_TEST_CASE(test_split_line) {
    const double tol = 1e-8;
    const unsigned int nseg = 2;
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(10, 10);

    umr::parametric::Line line(p0, p1);
    umr::math::optimize::SplitSegmentUniformly opt;

    opt.set_num_subsegments(nseg);
    std::vector<double> parameters = opt.run(&line);
    assert(parameters.size() == 3);

    BOOST_CHECK_CLOSE_FRACTION(parameters[0], 0, tol);
    BOOST_CHECK_CLOSE_FRACTION(parameters[1], 0.5, tol);
    BOOST_CHECK_CLOSE_FRACTION(parameters[2], 1, tol);
}


BOOST_AUTO_TEST_CASE(test_split_line_many) {
    const double tol = 1e-7;
    const unsigned int nseg_max = 10;
    umr::mesh::Point p0(-42, 6);
    umr::mesh::Point p1(24, 98);

    umr::parametric::Line line(p0, p1);
    umr::math::optimize::SplitSegmentUniformly opt;

    for (unsigned int nseg = 3; nseg < nseg_max; ++nseg) {
        opt.set_num_subsegments(nseg);
        std::vector<double> parameters = opt.run(&line);
        assert(parameters.size() == nseg + 1);

        for (unsigned int i = 0; i <= nseg; i++)
            BOOST_CHECK_CLOSE_FRACTION(parameters[i], (double)i/nseg, tol);
    }

}


BOOST_AUTO_TEST_CASE(test_split_arc) {
    const double tol = 1e-6;
    const double two_pi = 2 * std::numbers::pi;
    const unsigned int nseg = 3;

    umr::mesh::Point pc(0, 0);
    double radius = 1;
    double ti = 0;
    double tf = two_pi;
    umr::parametric::Arc arc(pc, radius, ti, tf);

    umr::math::optimize::SplitSegmentUniformly opt;
    opt.set_num_subsegments(nseg);
    std::vector<double> parameters = opt.run(&arc);
    assert(parameters.size() == nseg+1);
    BOOST_CHECK_SMALL(parameters[0], tol);
    BOOST_CHECK_CLOSE_FRACTION(parameters[1], two_pi * 1 / 3, tol);
    BOOST_CHECK_CLOSE_FRACTION(parameters[2], two_pi * 2 / 3, tol);
    BOOST_CHECK_CLOSE_FRACTION(parameters[3], two_pi, tol);
}


BOOST_AUTO_TEST_CASE(test_split_arc_many) {
    const double tol = 1e-6;
    const double two_pi = 2 * std::numbers::pi;
    const unsigned int nseg_max = 30;

    umr::mesh::Point pc(1, 2);
    double radius = 3;
    double ti = 0;
    double tf = two_pi;
    umr::parametric::Arc arc(pc, radius, ti, tf);

    umr::math::optimize::SplitSegmentUniformly opt;
    for (unsigned int nseg = 3; nseg < nseg_max; ++nseg) {
        opt.set_num_subsegments(nseg);
        std::vector<double> parameters = opt.run(&arc);
        assert(parameters.size() == nseg+1);

        BOOST_CHECK_SMALL(parameters[0], tol);
        for (unsigned int i = 1; i <= nseg; ++i)
            BOOST_CHECK_CLOSE_FRACTION(parameters[i], two_pi * i / nseg, tol);
    }
}


BOOST_AUTO_TEST_CASE(test_split_trig) {
    const double tol = 1e-6;
    const double two_pi = 2 * std::numbers::pi;
    const double n_periods = 4;
    double tmin = 0;
    double tmax = n_periods * two_pi;
    const unsigned int nseg = 5 * n_periods;

    auto function = [&](double t) -> std::pair<double, double> {
        assert(tmin - tol < t && t < tmax + tol);
        double x = t;
        double y = sin(t);
        return {x, y};
    };

    umr::parametric::Parametric parametric(function, tmin, tmax);

    umr::math::optimize::SplitSegmentUniformly opt;
    opt.set_num_subsegments(nseg);
    std::vector<double> parameters = opt.run(&parametric);
    assert(parameters.size() == nseg+1);

    std::vector<umr::mesh::Point> points;
    for (unsigned int i = 0; i < parameters.size(); i++) {
        auto [x, y] = function(parameters[i]);
        points.push_back(umr::mesh::Point(x, y));
    }

    // compare equal lengths between points since
    // parameters themselves won't be uniformly spaced
    // in nonlinear functions
    double lmin = std::numeric_limits<double>::max();
    double lmax = std::numeric_limits<double>::min();
    for (unsigned int i = 1; i < points.size(); i++) {
        double length = points[i].distance_to(points[i-1]);
        lmin = std::min(lmin, length);
        lmax = std::max(lmax, length);
    }
    BOOST_CHECK_CLOSE_FRACTION(lmin, lmax, tol);
}


BOOST_AUTO_TEST_SUITE_END()
