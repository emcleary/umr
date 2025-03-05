#include "common.hpp"
#include "../src/inequalities.hpp"


BOOST_AUTO_TEST_SUITE(test_inequalities)


BOOST_AUTO_TEST_CASE(test_close) {
    double x = 1;
    double rtol = 1e-6;
    double atol = 0;
    BOOST_CHECK(umr::inequalities::is_close(x, x, rtol, atol));
    BOOST_CHECK(umr::inequalities::is_close(x + rtol/2, x, rtol, atol));
    BOOST_CHECK(umr::inequalities::is_close(x - rtol/2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_close(x + rtol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_close(x - rtol*2, x, rtol, atol));

    rtol = 0;
    atol = 1e-6;
    BOOST_CHECK(umr::inequalities::is_close(x, x, rtol, atol));
    BOOST_CHECK(umr::inequalities::is_close(x + atol/2, x, rtol, atol));
    BOOST_CHECK(umr::inequalities::is_close(x - atol/2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_close(x + atol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_close(x - atol*2, x, rtol, atol));
}


BOOST_AUTO_TEST_CASE(test_lt) {
    double x = 1;
    double rtol = 1e-6;
    double atol = 0;
    BOOST_CHECK(!umr::inequalities::is_lt(x, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_lt(x, x + rtol*2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_lt(x, x + rtol/2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_lt(x - rtol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_lt(x - rtol/2, x, rtol, atol));

    rtol = 0;
    atol = 1e-6;
    BOOST_CHECK(!umr::inequalities::is_lt(x, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_lt(x, x + atol*2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_lt(x, x + atol/2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_lt(x - atol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_lt(x - atol/2, x, rtol, atol));
}


BOOST_AUTO_TEST_CASE(test_le) {
    double x = 1;
    double rtol = 1e-6;
    double atol = 0;
    BOOST_CHECK( umr::inequalities::is_le(x, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_le(x, x - rtol*2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_le(x, x - rtol/2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_le(x + rtol*2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_le(x + rtol/2, x, rtol, atol));

    rtol = 0;
    atol = 1e-6;
    BOOST_CHECK( umr::inequalities::is_le(x, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_le(x, x - atol*2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_le(x, x - atol/2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_le(x + atol*2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_le(x + atol/2, x, rtol, atol));
}


BOOST_AUTO_TEST_CASE(test_gt) {
    double x = 1;
    double rtol = 1e-6;
    double atol = 0;
    BOOST_CHECK(!umr::inequalities::is_gt(x, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_gt(x + rtol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_gt(x + rtol/2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_gt(x, x - rtol*2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_gt(x, x - rtol/2, rtol, atol));

    rtol = 0;
    atol = 1e-6;
    BOOST_CHECK(!umr::inequalities::is_gt(x, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_gt(x + atol*2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_gt(x + atol/2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_gt(x, x - atol*2, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_gt(x, x - atol/2, rtol, atol));
}


BOOST_AUTO_TEST_CASE(test_ge) {
    double x = 1;
    double rtol = 1e-6;
    double atol = 0;
    BOOST_CHECK( umr::inequalities::is_ge(x, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_ge(x - rtol*2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_ge(x - rtol/2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_ge(x, x + rtol*2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_ge(x, x + rtol/2, rtol, atol));

    rtol = 0;
    atol = 1e-6;
    BOOST_CHECK( umr::inequalities::is_ge(x, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_ge(x - atol*2, x, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_ge(x - atol/2, x, rtol, atol));
    BOOST_CHECK(!umr::inequalities::is_ge(x, x + atol*2, rtol, atol));
    BOOST_CHECK( umr::inequalities::is_ge(x, x + atol/2, rtol, atol));
}


BOOST_AUTO_TEST_SUITE_END()
