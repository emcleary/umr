#include <cassert>
#include <numbers>
#include <cmath>

#include "common.hpp"
#include "../src/parametrics.hpp"


BOOST_AUTO_TEST_SUITE(test_parametrics)


BOOST_AUTO_TEST_CASE(test_line) {
    const double tol = 1e-8;
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(10, 10);
    umr::parametric::Line line(p0, p1);

    auto [x0, y0] = line.evaluate(0);
    BOOST_CHECK_EQUAL(x0, p0.x);
    BOOST_CHECK_EQUAL(y0, p0.y);

    auto [x1, y1] = line.evaluate(0.3);
    BOOST_CHECK_CLOSE_FRACTION(x1, 3, tol);
    BOOST_CHECK_CLOSE_FRACTION(y1, 3, tol);

    auto [x2, y2] = line.evaluate(0.5);
    BOOST_CHECK_CLOSE_FRACTION(x2, 5, tol);
    BOOST_CHECK_CLOSE_FRACTION(y2, 5, tol);

    auto [x3, y3] = line.evaluate(1);
    BOOST_CHECK_EQUAL(x3, p1.x);
    BOOST_CHECK_EQUAL(y3, p1.y);
}


BOOST_AUTO_TEST_CASE(test_line_many) {
    const double tol = 1e-7;
    umr::mesh::Point p0(-42, 6);
    umr::mesh::Point p1(24, 98);
    umr::parametric::Line line(p0, p1);
    std::vector<double> t = {0.1, 0.234, 0.65, 0.999};

    auto [xmin, ymin] = line.evaluate(0);
    BOOST_CHECK_EQUAL(xmin, p0.x);
    BOOST_CHECK_EQUAL(ymin, p0.y);

    auto [xmax, ymax] = line.evaluate(1);
    BOOST_CHECK_EQUAL(xmax, p1.x);
    BOOST_CHECK_EQUAL(ymax, p1.y);

    double dx = p1.x - p0.x;
    double dy = p1.y - p0.y;
    for (auto ti : t) {
        auto [x, y] = line.evaluate(ti);
        BOOST_CHECK_CLOSE_FRACTION(x, p0.x + dx*ti, tol);
        BOOST_CHECK_CLOSE_FRACTION(y, p0.y + dy*ti, tol);
    }
}


BOOST_AUTO_TEST_CASE(test_arc) {
    const double tol = 1e-8;

    umr::mesh::Point pc(0, 0);
    double radius = 1;
    double ti = 0;
    double tf = std::numbers::pi;
    umr::parametric::Arc arc(pc, radius, ti, tf);

    auto [x0, y0] = arc.evaluate(0);
    BOOST_CHECK_CLOSE_FRACTION(x0, 1, tol);
    BOOST_CHECK_SMALL(y0, tol);
    
    auto [x1, y1] = arc.evaluate(1./4 * std::numbers::pi);
    BOOST_CHECK_CLOSE_FRACTION(x1, 1 / sqrt(2), tol);
    BOOST_CHECK_CLOSE_FRACTION(y1, 1 / sqrt(2), tol);
    
    auto [x2, y2] = arc.evaluate(1./2 * std::numbers::pi);
    BOOST_CHECK_SMALL(x2, tol);
    BOOST_CHECK_CLOSE_FRACTION(y2, 1, tol);

    auto [x3, y3] = arc.evaluate(3./4 * std::numbers::pi);
    BOOST_CHECK_CLOSE_FRACTION(x3, -1 / sqrt(2), tol);
    BOOST_CHECK_CLOSE_FRACTION(y3, 1 / sqrt(2), tol);

    auto [x4, y4] = arc.evaluate(std::numbers::pi);
    BOOST_CHECK_CLOSE_FRACTION(x4, -1, tol);
    BOOST_CHECK_SMALL(y4, tol);
}


BOOST_AUTO_TEST_CASE(test_arc_many) {
    const double tol = 1e-8;
    const double two_pi = 2 * std::numbers::pi;

    umr::mesh::Point pc(1, 2);
    double radius = 3;
    double ti = 1./8 * two_pi;
    double tf = 7./8 * two_pi;
    // scaled between [0, 1] for simplicity
    std::vector<double> parameters = {0, 0.1, 0.234, 0.5, 0.654, 0.999, 1};

    umr::parametric::Arc arc(pc, radius, ti, tf);

    auto parametric = [&](double s) -> std::pair<double, double> {
        assert(-tol < s && s < 1 + tol);
        double t = ti * (1-s) + tf * s;
        return arc.evaluate(t);
    };

    auto [x0, y0] = arc.evaluate(ti);
    auto analytical = [&](double s) -> std::pair<double, double> {
        assert(-tol < s && s < 1 + tol);
        double theta = s * (tf - ti);
        double x = x0 - pc.x;
        double y = y0 - pc.y;
        double cosine = cos(theta);
        double sine = sin(theta);
        double xs = x * cosine - y * sine;
        double ys = x * sine + y * cosine;
        return {xs + pc.x, ys + pc.y};
    };

    for (double si : parameters) {
        auto [xa, ya] = analytical(si);
        auto [xp, yp] = parametric(si);
        BOOST_CHECK_CLOSE_FRACTION(xa, xp, tol);
        BOOST_CHECK_CLOSE_FRACTION(ya, yp, tol);
    }
}


BOOST_AUTO_TEST_CASE(test_polynomial) {
    const double tol = 1e-8;
    double tmin = 3;
    double tmax = 10;

    auto polynomial = [&](double t) -> std::pair<double, double> {
        assert(tmin - tol < t && t < tmax + tol);
        double x = t + 1;
        double y = t * t + t;
        return {x, y};
    };

    umr::parametric::Parametric parametric(polynomial, tmin, tmax);
    
    auto [x0, y0] = parametric.evaluate(tmin);
    auto [xp0, yp0] = polynomial(tmin);
    BOOST_CHECK_EQUAL(x0, xp0);
    BOOST_CHECK_EQUAL(y0, yp0);

    auto [x1, y1] = parametric.evaluate(tmax);
    auto [xp1, yp1] = polynomial(tmax);
    BOOST_CHECK_EQUAL(x1, xp1);
    BOOST_CHECK_EQUAL(y1, yp1);

    unsigned int nmax = 10;
    for (unsigned int i = 1; i < nmax; i++) {
        double t = tmin + (tmax - tmin) * i / nmax;
        std::cout << "t = " << t << '\n';
        auto [x, y] = parametric.evaluate(t);
        auto [xp, yp] = polynomial(t);
        BOOST_CHECK_CLOSE_FRACTION(x, xp, tol);
        BOOST_CHECK_CLOSE_FRACTION(y, yp, tol);
    }
}


BOOST_AUTO_TEST_CASE(test_periodic) {
    const double tol = 1e-8;
    const double two_pi = 2 * std::numbers::pi;
    double tmin = 0;
    double tmax = two_pi;
    double radius = 1;

    auto function = [&](double t) -> std::pair<double, double> {
        assert(tmin - tol < t && t < tmax + tol);
        double theta = t;
        double x = radius * cos(theta);
        double y = radius * sin(theta);
        return {x, y};
    };

    umr::parametric::Parametric parametric(function, tmin, tmax);
    
    auto [x0, y0] = parametric.evaluate(tmin);
    BOOST_CHECK_SMALL(x0 - 1, tol);
    BOOST_CHECK_SMALL(y0 - 0, tol);

    auto [x1, y1] = parametric.evaluate(tmax);
    BOOST_CHECK_SMALL(x1 - 1, tol);
    BOOST_CHECK_SMALL(y1 - 0, tol);

    // Since periodic, endpoints are expected to be EXACTLY equal
    BOOST_CHECK_EQUAL(x1, x0);
    BOOST_CHECK_EQUAL(y1, y0);

    unsigned int nmax = 10;
    for (unsigned int i = 1; i < nmax; i++) {
        double t = tmin + (tmax - tmin) * i / nmax;
        std::cout << "t = " << t << '\n';
        auto [x, y] = parametric.evaluate(t);
        auto [xp, yp] = function(t);
        BOOST_CHECK_CLOSE_FRACTION(x, xp, tol);
        BOOST_CHECK_CLOSE_FRACTION(y, yp, tol);
    }
}


BOOST_AUTO_TEST_SUITE_END()
