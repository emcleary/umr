#include <cmath>
#include "common.hpp"
#include "../src/triangle.hpp"

BOOST_AUTO_TEST_SUITE(test_point)


BOOST_AUTO_TEST_CASE(test_triangle) {
    // 30-60-90 triangle
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(0, sqrt(3));
    umr::mesh::Triangle triangle(p0, p1, p2);

    double tol = 1e-8;
    BOOST_CHECK_CLOSE(triangle.get_min_length(), 1, tol);
    BOOST_CHECK_CLOSE(triangle.get_max_length(), 2, tol);
    BOOST_CHECK_CLOSE(triangle.get_min_angle(), 30, tol);
    BOOST_CHECK_CLOSE(triangle.get_max_angle(), 90, tol);

    const umr::mesh::Point& circumcenter = triangle.get_circumcenter();
    BOOST_CHECK_CLOSE(circumcenter.x, 0.5, tol);
    BOOST_CHECK_CLOSE(circumcenter.y, sqrt(3) / 2, tol);

    BOOST_CHECK_CLOSE(triangle.get_radius(), 1, tol);
}


BOOST_AUTO_TEST_CASE(test_triangle_2) {
    // 45-45-90 triangle
    umr::mesh::Point p0 = umr::mesh::Point(1, 0);
    umr::mesh::Point p1 = umr::mesh::Point(0, 1);
    umr::mesh::Point p2 = umr::mesh::Point(-1, 0);
    umr::mesh::Triangle triangle(p0, p1, p2);

    double tol = 1e-8;
    BOOST_CHECK_CLOSE(triangle.get_min_length(), sqrt(2), tol);
    BOOST_CHECK_CLOSE(triangle.get_max_length(), 2, tol);
    BOOST_CHECK_CLOSE(triangle.get_min_angle(), 45, tol);
    BOOST_CHECK_CLOSE(triangle.get_max_angle(), 90, tol);

    const umr::mesh::Point& circumcenter = triangle.get_circumcenter();
    BOOST_CHECK_CLOSE(circumcenter.x, 0.0, tol);
    BOOST_CHECK_CLOSE(circumcenter.y, 0.0, tol);

    BOOST_CHECK_CLOSE(triangle.get_radius(), 1, tol);
}


BOOST_AUTO_TEST_SUITE_END()
