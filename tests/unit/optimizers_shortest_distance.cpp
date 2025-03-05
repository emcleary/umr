#include <cassert>
#include <memory>

#include "common.hpp"
#include "../src/optimizers.hpp"
#include "../src/point.hpp"
#include "../src/sources.hpp"

BOOST_AUTO_TEST_SUITE(test_optimizers_shortest_distance)


double calculate_distance(std::shared_ptr<umr::source::IShape> s0, std::shared_ptr<umr::source::IShape> s1) {
    umr::source::SourceManager sm0;
    sm0.add_shape(s0);
    sm0.finalize();

    umr::source::SourceManager sm1;
    sm1.add_shape(s1);
    sm1.finalize();

    umr::math::optimize::ShortestDistanceParametricParametric opt;
    opt.set_num_iter(1000);
    opt.set_epsilon(1e-14);
    double dmin = std::numeric_limits<double>::max();
    for (umr::parametric::IParametric* p0 : sm0.get_segments()) {
        opt.set_parametric_0(p0);

        for (umr::parametric::IParametric* p1 : sm1.get_segments()) {
            opt.set_parametric_1(p1);

            std::vector<double> parameters = opt.run(&opt);
            auto [x0, y0] = p0->evaluate(parameters[0]);
            auto [x1, y1] = p1->evaluate(parameters[1]);
            double dx = x1 - x0;
            double dy = y1 - y0;
            double d = sqrt(dx * dx + dy * dy);
            dmin = std::min(dmin, d);
        }
    }

    return dmin;
}

double calculate_distance(umr::parametric::IParametric* parametric0, umr::parametric::IParametric* parametric1) {
    umr::math::optimize::ShortestDistanceParametricParametric opt;
    opt.set_parametric_0(parametric0);
    opt.set_parametric_1(parametric1);

    std::vector<double> parameters = opt.run(&opt);
    assert(parameters.size() == 2);
    
    auto [x0, y0] = parametric0->evaluate(parameters[0]);
    auto [x1, y1] = parametric1->evaluate(parameters[1]);
    double dx = x1 - x0;
    double dy = y1 - y0;
    return sqrt(dx * dx + dy * dy);
}

double calculate_distance(umr::parametric::IParametric* parametric, umr::mesh::Point* point) {
    umr::math::optimize::ShortestDistanceParametricPoint opt;
    opt.set_parametric(parametric);
    opt.set_point(point);

    std::vector<double> parameters = opt.run(&opt);
    assert(parameters.size() == 1);
    
    auto [x, y] = parametric->evaluate(parameters[0]);
    double dx = x - point->x;
    double dy = y - point->y;
    return sqrt(dx * dx + dy * dy);
}


BOOST_AUTO_TEST_CASE(test_func_func) {
    const double tol = 1e-5;
    const unsigned int n_periods = 4;
    const double tmin = 0;
    const double tmax = n_periods * 2 * std::numbers::pi;

    auto func0 = [](double t) -> std::pair<double, double> {
        double x = t;
        double y = sin(t);
        return {x, y};
    };
    
    auto func1 = [&](double t) -> std::pair<double, double> {
        return {t, 5};
    };

    umr::mesh::Point p0(0, 5);
    umr::mesh::Point p1(2 * std::numbers::pi, 5);
    umr::parametric::Line line01(p0, p1);

    umr::parametric::Parametric parametric0(func0, tmin, tmax);
    umr::parametric::Parametric parametric1(func1, tmin, tmax);

    umr::mesh::Point p2(0, 0);
    umr::mesh::Point p3(0, -2);
    umr::mesh::Point p4(tmax, -2);
    umr::mesh::Point p5(tmax, 0);
    auto s_wave = std::make_shared<umr::source::Shape>();
    s_wave->add_parametric(func0, tmin, tmax, n_periods);
    s_wave->add_line(p5, p4);
    s_wave->add_line(p4, p3);
    s_wave->add_line(p3, p2);

    umr::mesh::Point p6(0, 6);
    auto s_tri = std::make_shared<umr::source::ShapeTriangle>(p0, p1, p6);

    // nonintersecting
    double d = calculate_distance(&parametric0, &parametric1);
    BOOST_CHECK_CLOSE_FRACTION(d, 4, tol);

    // stuck at local minimum
    d = calculate_distance(&parametric0, &line01);
    BOOST_TEST(d != 4, tt::tolerance(tol));

    // finds global minimum with sufficient splitting
    d = calculate_distance(s_wave, s_tri);
    BOOST_CHECK_CLOSE_FRACTION(d, 4, tol);
}


BOOST_AUTO_TEST_CASE(test_line_point) {
    const double tol = 1e-5;
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(2, 0);
    umr::mesh::Point p2(1, -1);
    umr::mesh::Point p3(1, 0);
    umr::mesh::Point p4(1, 1);

    umr::parametric::Line line(p0, p1);

    // nonintersecting
    double d = calculate_distance(&line, &p2);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = calculate_distance(&line, &p4);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    // intersecting
    d = calculate_distance(&line, &p0);
    BOOST_CHECK_SMALL(d, tol);

    d = calculate_distance(&line, &p1);
    BOOST_CHECK_SMALL(d, tol);

    d = calculate_distance(&line, &p3);
    BOOST_CHECK_SMALL(d, tol);
}


BOOST_AUTO_TEST_CASE(test_arc_point) {
    const double tol = 1e-5;
    umr::mesh::Point p0(0, 0);

    umr::mesh::Point pc0(0, 1);
    umr::mesh::Point pc1(0, 2);
    double radius = 1;
    umr::parametric::Arc arc0(pc0, radius, 0, std::numbers::pi); 
    umr::parametric::Arc arc1(pc1, radius, std::numbers::pi, 2*std::numbers::pi); 

    // nonintersecting -- nonunique
    double d = calculate_distance(&arc0, &p0);
    BOOST_CHECK_CLOSE_FRACTION(d, sqrt(2), tol);

    // nonintersecting -- unique
    d = calculate_distance(&arc1, &p0);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    // intersecting
    d = calculate_distance(&arc1, &pc1);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);
}

/*
 * NB: These arc-arc tests are a bit inconsistent. I believe the main
 * reason is that solutions can be non-unique, as tested a bit more
 * carefully in other tests. The important thing is that the
 * circle-circle cases seems to work fine since they are split into
 * subsegments. Only optimizing with shape subsegments is done in the
 * code; optimizing parametrics themselves is not.
 */
BOOST_AUTO_TEST_CASE(test_arc_arc) {
    const double tol = 1e-8;
    const double tmin = 0;
    const double tmax = 2 * std::numbers::pi;
    
    umr::mesh::Point pc0(0, 0);
    umr::mesh::Point pc1(1, 0);
    umr::mesh::Point pc2(2, 0);

    umr::parametric::Arc arc0(pc0, 1, tmin, tmax);
    umr::parametric::Arc arc1(pc0, 2, tmin, tmax);
    umr::parametric::Arc arc2(pc1, 2, tmin, tmax);
    umr::parametric::Arc arc3(pc1, 3, tmin, tmax);

    auto circ0 = std::make_shared<umr::source::ShapeCircle>(pc0, 1);
    auto circ1 = std::make_shared<umr::source::ShapeCircle>(pc1, 1);
    
    // contained circles
    double d = calculate_distance(&arc0, &arc1);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = calculate_distance(&arc0, &arc3);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    // tangent
    d = calculate_distance(&arc0, &arc2);
    BOOST_CHECK_SMALL(d, tol);

    // this test is notably difficult; adding a phase shift to the
    // parameter helps
    double phase = std::numbers::pi / 2;
    umr::parametric::Arc arcint(pc1, 1, tmin + phase, tmax + phase);
    d = calculate_distance(&arc0, &arcint);
    BOOST_CHECK_SMALL(d, tol);

    // intersecting -- WORKS when using sources to split into subsegments!
    circ0->set_number_subsegments({3});
    circ1->set_number_subsegments({3});
    d = calculate_distance(circ0, circ1);
    BOOST_CHECK_SMALL(d, tol);
}

// lines with unique solutions of closest distance
BOOST_AUTO_TEST_CASE(test_lines) {
    const double tol = 1e-5;

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(5, 0);
    umr::mesh::Point p2(6, -1);
    umr::mesh::Point p3(6, 1);
    umr::mesh::Point p4(3, 1);
    umr::mesh::Point p5(5, 6);
    umr::mesh::Point p6(0, -1);
    umr::mesh::Point p7(5, 1);
    umr::mesh::Point p8(0, 1);

    umr::parametric::Line line01 = umr::parametric::Line(p0, p1);
    umr::parametric::Line line23 = umr::parametric::Line(p2, p3);
    umr::parametric::Line line45 = umr::parametric::Line(p4, p5);
    umr::parametric::Line line67 = umr::parametric::Line(p6, p7);
    umr::parametric::Line line68 = umr::parametric::Line(p6, p8);

    // nonintersecting
    double d = calculate_distance(&line01, &line23);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    d = calculate_distance(&line01, &line45);
    BOOST_CHECK_CLOSE_FRACTION(d, 1, tol);

    // intersecting
    d = calculate_distance(&line01, &line67);
    BOOST_CHECK_SMALL(d, tol);

    d = calculate_distance(&line01, &line68);
    BOOST_CHECK_SMALL(d, tol);
}

// separate tests since solutions aren't (always) unique
BOOST_AUTO_TEST_CASE(test_parallel_lines) {
    const double tol = 1e-4;
    const double y = 1;

    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(-2, y);
    umr::mesh::Point p3(-1, y);
    umr::mesh::Point p4(0, y);
    umr::mesh::Point p5(1, y);
    umr::mesh::Point p6(2, y);
    umr::mesh::Point p7(3, y);

    umr::parametric::Line line01 = umr::parametric::Line(p0, p1);
    umr::parametric::Line line23 = umr::parametric::Line(p2, p3);
    umr::parametric::Line line34 = umr::parametric::Line(p3, p4);
    umr::parametric::Line line45 = umr::parametric::Line(p4, p5);
    umr::parametric::Line line56 = umr::parametric::Line(p5, p6);
    umr::parametric::Line line67 = umr::parametric::Line(p6, p7);
    umr::parametric::Line line46 = umr::parametric::Line(p4, p6);
    umr::parametric::Line line47 = umr::parametric::Line(p4, p7);
    umr::parametric::Line line57 = umr::parametric::Line(p5, p7);

    // neighboring lines
    double d = calculate_distance(&line01, &line23);
    BOOST_CHECK_CLOSE_FRACTION(d, p3.distance_to(p0), tol);

    d = calculate_distance(&line01, &line34);
    BOOST_CHECK_CLOSE_FRACTION(d, y, tol);

    d = calculate_distance(&line01, &line45);
    BOOST_CHECK_CLOSE_FRACTION(d, y, tol);

    d = calculate_distance(&line01, &line56);
    BOOST_CHECK_CLOSE_FRACTION(d, y, tol);

    d = calculate_distance(&line01, &line67);
    BOOST_CHECK_CLOSE_FRACTION(d, p6.distance_to(p1), tol);

    // consecutive lines
    d = calculate_distance(&line34, &line56);
    BOOST_CHECK_CLOSE_FRACTION(d, p4.distance_to(p5), tol);

    // overlapping lines
    d = calculate_distance(&line56, &line47);
    BOOST_CHECK_SMALL(d, tol);

    d = calculate_distance(&line46, &line57);
    BOOST_CHECK_SMALL(d, tol);

    d = calculate_distance(&line45, &line56);
    BOOST_CHECK_SMALL(d, tol);
}


BOOST_AUTO_TEST_SUITE_END()
