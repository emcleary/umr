#include <cmath>

#include "common.hpp"
#include "../src/quadedge.hpp"
#include "../src/math.hpp"


BOOST_AUTO_TEST_SUITE(test_quadedge_math)


BOOST_AUTO_TEST_CASE(test_triangle_determinant) {
    umr::mesh::Point p0 = umr::mesh::Point(0, 0);
    umr::mesh::Point p1 = umr::mesh::Point(1, 0);
    umr::mesh::Point p2 = umr::mesh::Point(0, 1);

    double det = umr::math::triangle_determinant(p0, p1, p2);
    BOOST_CHECK(det > 0);

    det = umr::math::triangle_determinant(p0, p2, p1);
    BOOST_CHECK(det < 0);

    umr::mesh::Point p3 = umr::mesh::Point(0, 1);
    umr::mesh::Point p4 = umr::mesh::Point(0, 2);
    umr::mesh::Point p5 = umr::mesh::Point(0, 3);
    umr::mesh::Point p6 = umr::mesh::Point(1, 0);
    umr::mesh::Point p7 = umr::mesh::Point(2, 0);
    umr::mesh::Point p8 = umr::mesh::Point(3, 0);

    det = umr::math::triangle_determinant(p3, p4, p5);
    BOOST_CHECK(det == 0);
    det = umr::math::triangle_determinant(p4, p5, p3);
    BOOST_CHECK(det == 0);
    det = umr::math::triangle_determinant(p3, p5, p4);
    BOOST_CHECK(det == 0);

    det = umr::math::triangle_determinant(p6, p7, p8);
    BOOST_CHECK(det == 0);
    det = umr::math::triangle_determinant(p7, p8, p6);
    BOOST_CHECK(det == 0);
    det = umr::math::triangle_determinant(p6, p8, p7);
    BOOST_CHECK(det == 0);
}


BOOST_AUTO_TEST_CASE(test_counter_clockwise_points) {
    umr::mesh::Point p0 = umr::mesh::Point(0, 0);
    umr::mesh::Point p1 = umr::mesh::Point(1, 0);
    umr::mesh::Point p2 = umr::mesh::Point(0, 1);

    bool ccw = umr::math::is_ccw(p0, p1, p2);
    BOOST_CHECK(ccw);

    ccw = umr::math::is_ccw(p0, p2, p1);
    BOOST_CHECK(!ccw);

    ccw = umr::math::is_ccw(p0, p0, p0);
    BOOST_CHECK(!ccw);
}


BOOST_AUTO_TEST_CASE(test_left_right_on_line) {
    umr::mesh::Point p0 = umr::mesh::Point(0, 0);
    umr::mesh::Point p1 = umr::mesh::Point(10, 0);
    umr::mesh::Point p2 = umr::mesh::Point(5, 1);
    umr::mesh::Point p3 = umr::mesh::Point(5, 0);
    umr::mesh::Edge e01 = umr::mesh::quadedge::make_edge(p0, p1);

    BOOST_CHECK(umr::math::left_of(p2, p0, p1));
    BOOST_CHECK(umr::math::left_of(p2, e01));
    BOOST_CHECK(!umr::math::left_of(p2, e01.sym()));

    BOOST_CHECK(umr::math::right_of(p2, p1, p0));
    BOOST_CHECK(!umr::math::right_of(p2, e01));
    BOOST_CHECK(umr::math::right_of(p2, e01.sym()));

    BOOST_CHECK(!umr::math::on_line(p2, p0, p1));
    BOOST_CHECK(!umr::math::on_line(p2, p1, p0));
    BOOST_CHECK(!umr::math::on_line(p2, e01));
    BOOST_CHECK(!umr::math::on_line(p2, e01.sym()));

    BOOST_CHECK(umr::math::on_line(p3, p0, p1));
    BOOST_CHECK(umr::math::on_line(p3, p1, p0));
    BOOST_CHECK(umr::math::on_line(p3, e01));
    BOOST_CHECK(umr::math::on_line(p3, e01.sym()));
}


BOOST_AUTO_TEST_CASE(test_circumcenter) {
    umr::mesh::Point p0 = umr::mesh::Point(1, 0);
    umr::mesh::Point p1 = umr::mesh::Point(0, 1);
    umr::mesh::Point p2 = umr::mesh::Point(-1, 0);
    auto [x0, y0] = umr::math::calculate_circumcenter(p0, p1, p2);
    BOOST_CHECK_CLOSE(x0, 0, 1e-30);
    BOOST_CHECK_CLOSE(y0, 0, 1e-30);

    // Calculation involves division by a determinant which will be 0
    // if the 3 points are on the same line.
    umr::mesh::Point p3 = umr::mesh::Point(1, 1);
    umr::mesh::Point p4 = umr::mesh::Point(2, 2);
    umr::mesh::Point p5 = umr::mesh::Point(3, 3);
    auto [x1, y1] = umr::math::calculate_circumcenter(p3, p4, p5);
    BOOST_TEST_CHECK(std::isinf(x1));
    BOOST_TEST_CHECK(std::isinf(y1));
}


BOOST_AUTO_TEST_CASE(test_in_circle) {
    umr::mesh::Point p0 = umr::mesh::Point(1, 0);
    umr::mesh::Point p1 = umr::mesh::Point(0, 1);
    umr::mesh::Point p2 = umr::mesh::Point(-1, 0);

    double tol = 1e-10;
    umr::mesh::Point q0 = umr::mesh::Point(0, 0);
    umr::mesh::Point q1 = umr::mesh::Point(1-tol, 0);
    umr::mesh::Point q2 = umr::mesh::Point(1, 0);
    umr::mesh::Point q3 = umr::mesh::Point(1+tol, 0);
    umr::mesh::Point q4 = umr::mesh::Point(2, 0);
    BOOST_CHECK(umr::math::in_circle(p0, p1, p2, q0));
    BOOST_CHECK(umr::math::in_circle(p0, p1, p2, q1));
    BOOST_CHECK(!umr::math::in_circle(p0, p1, p2, q2));
    BOOST_CHECK(!umr::math::in_circle(p0, p1, p2, q3));
    BOOST_CHECK(!umr::math::in_circle(p0, p1, p2, q4));
}


BOOST_AUTO_TEST_CASE(test_intersect) {
    umr::mesh::Point p0 = umr::mesh::Point(-1, 0);
    umr::mesh::Point p1 = umr::mesh::Point(1, 0);
    umr::mesh::Edge e01 = umr::mesh::quadedge::make_edge(p0, p1);

    umr::mesh::Point p2 = umr::mesh::Point(0, -1);
    umr::mesh::Point p3 = umr::mesh::Point(0, 1);
    umr::mesh::Edge e23 = umr::mesh::quadedge::make_edge(p2, p3);

    umr::mesh::Point p4 = umr::mesh::Point(1, 1);
    umr::mesh::Edge e14 = umr::mesh::quadedge::make_edge(p1, p4);

    umr::mesh::Point p5 = umr::mesh::Point(10, -1);
    umr::mesh::Point p6 = umr::mesh::Point(-1, 10);
    umr::mesh::Edge e56 = umr::mesh::quadedge::make_edge(p5, p6);

    umr::mesh::Point p7 = umr::mesh::Point(2, 0);
    umr::mesh::Point p8 = umr::mesh::Point(3, 0);
    umr::mesh::Edge e17 = umr::mesh::quadedge::make_edge(p1, p7);
    umr::mesh::Edge e78 = umr::mesh::quadedge::make_edge(p7, p8);

    umr::mesh::Point p9 = umr::mesh::Point(0, 0);
    umr::mesh::Edge e97 = umr::mesh::quadedge::make_edge(p9, p7);

    BOOST_CHECK( umr::math::intersect(e01, e23)); // edges cross at the origin
    BOOST_CHECK( umr::math::intersect(e01, e14)); // edges share point p1
    BOOST_CHECK(!umr::math::intersect(e01, e56)); // lines don't cross; e56 is on both sides of e01
    BOOST_CHECK(!umr::math::intersect(e23, e56)); // lines don't cross; e23 is on both sides of e56
    BOOST_CHECK( umr::math::intersect(e01, e17)); // 2 line segments on the same line with a common point
    BOOST_CHECK(!umr::math::intersect(e01, e78)); // 2 non-overlapping line segments on the same line
    BOOST_CHECK(!umr::math::intersect(e01, e97)); // 2 overlapping line segments on the same line;
                                                  // NB: should be TRUE, but edges should never be like this
                                                  //     in this application
}


BOOST_AUTO_TEST_CASE(test_edges) {
    umr::mesh::Point p0 = umr::mesh::Point(1, 0);
    umr::mesh::Point p1 = umr::mesh::Point(0, 1);
    umr::mesh::Point p2 = umr::mesh::Point(-1, 0);
    umr::mesh::Edge e0 = umr::mesh::quadedge::make_edge(p0, p1);
    BOOST_CHECK(!umr::mesh::quadedge::part_of_triangle(e0));

    umr::mesh::Edge e1 = umr::mesh::quadedge::extend_edge(e0, p2);
    BOOST_CHECK(!umr::mesh::quadedge::part_of_triangle(e0));
    BOOST_CHECK(!umr::mesh::quadedge::part_of_triangle(e1));

    umr::mesh::Edge e2 = umr::mesh::quadedge::connect(e1, e0);
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e0));
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e1));
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e2));
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e0.sym()));
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e1.sym()));
    BOOST_CHECK(umr::mesh::quadedge::part_of_triangle(e2.sym()));
}


BOOST_AUTO_TEST_CASE(test_law_of_cosines) {
    double tol = 1e-8;
    BOOST_CHECK_CLOSE(umr::math::loc_angle_deg(1, 1, 1), 60, tol);
    BOOST_CHECK_CLOSE(umr::math::loc_angle(1, 1, 1), std::numbers::pi / 3, tol);
}


BOOST_AUTO_TEST_SUITE_END()
