#include <cmath>
#include "common.hpp"
#include "../src/point.hpp"

BOOST_AUTO_TEST_SUITE(test_point)


BOOST_AUTO_TEST_CASE(test_distance) {
    umr::mesh::Point p0(0, 0);
    umr::mesh::Point p1(1, 0);
    umr::mesh::Point p2(0, 1);
    umr::mesh::Point p3(1, 1);

    double d00 = p0.distance_to(p0);
    double d01 = p0.distance_to(p1);
    double d02 = p0.distance_to(p2);
    double d03 = p0.distance_to(p3);

    BOOST_CHECK(d00 == 0);
    BOOST_CHECK(d01 == 1);
    BOOST_CHECK(d02 == 1);
    BOOST_CHECK(d03 == sqrt(2));

    BOOST_CHECK(d01 == p1.distance_to(p0));
    BOOST_CHECK(d02 == p2.distance_to(p0));
    BOOST_CHECK(d03 == p3.distance_to(p0));

    BOOST_CHECK(d00 == p0.distance_to(p0.x, p0.y));
    BOOST_CHECK(d01 == p0.distance_to(p1.x, p1.y));
    BOOST_CHECK(d02 == p0.distance_to(p2.x, p2.y));
    BOOST_CHECK(d03 == p0.distance_to(p3.x, p3.y));
}


BOOST_AUTO_TEST_SUITE_END()
