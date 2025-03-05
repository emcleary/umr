#include <cmath> // abs
#include <cassert>

#include "math.hpp"
#include "inequalities.hpp"


namespace umr {

namespace math {


double triangle_determinant(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2) {
    double det = p0.x * (p1.y - p2.y);
    det += p1.x * (p2.y - p0.y);
    det += p2.x * (p0.y - p1.y);
    return det;
};


std::pair<double, double> calculate_circumcenter(const mesh::Point& p0,
        const mesh::Point& p1,
        const mesh::Point& p2) {
    double x1 = p1.x - p0.x;
    double y1 = p1.y - p0.y;
    double x2 = p2.x - p0.x;
    double y2 = p2.y - p0.y;
    double r1 = x1 * (p1.x + p0.x) + y1 * (p1.y + p0.y);
    double r2 = x2 * (p2.x + p0.x) + y2 * (p2.y + p0.y);
    double det = x1 * y2 - y1 * x2;
    double xc = (y2 * r1 - y1 * r2) / det / 2;
    double yc = (x1 * r2 - x2 * r1) / det / 2;
    return {xc, yc};
}


bool is_ccw(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2) {
    double det = triangle_determinant(p0, p1, p2);
    return inequalities::is_gt(det, 0);
};


bool left_of(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1) {
    return is_ccw(p0, p1, q);
};


bool left_of(const mesh::Point& q, const mesh::Edge& e) {
    return left_of(q, e.org(), e.dest());
};


bool right_of(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1) {
    return is_ccw(p1, p0, q);
};


bool right_of(const mesh::Point& q, const mesh::Edge& e) {
    return right_of(q, e.org(), e.dest());
};


bool on_line(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1) {
    double det = triangle_determinant(p0, p1, q);
    return inequalities::is_close(det, 0);
};


bool on_line(const mesh::Point& q, const mesh::Edge& e) {
    return on_line(e.org(), e.dest(), q);
};


// is point q within the circumcircle of points p0, p1, p2?
bool in_circle(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2, const mesh::Point& q) {
    if (on_line(p0, p1, p2))
        return false;

    auto [xc, yc] = calculate_circumcenter(p0, p1, p2);

    auto calc_r2 = [&](const mesh::Point& p) {
        double dx = p.x - xc;
        double dy = p.y - yc;
        return dx * dx + dy * dy;
    };

    double r2_0 = calc_r2(p0);
    double r2_q = calc_r2(q);

    static const double rtol = 1e-12;
    static const double atol = 0; // assume r2_0 is ALWAYS significantly greater than rtol
    bool inside = inequalities::is_lt(r2_q, r2_0, rtol, atol);
    assert(inequalities::is_close(calc_r2(p1), r2_0) && "in_circle: r1 distance not close to r0");
    assert(inequalities::is_close(calc_r2(p2), r2_0) && "in_circle: r2 distance not close to r0");
    return inside;
}


bool in_diametral_circle(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& q) {
    double xmid = (p0.x + p1.x) / 2;
    double ymid = (p0.y + p1.y) / 2;
    double radius = p0.distance_to(xmid, ymid);
    double distance = q.distance_to(xmid, ymid);
    return inequalities::is_lt(distance, radius);
}


bool intersect(const mesh::Edge& e0, const mesh::Edge& e1) {
    auto det_sign = [&](const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2) {
        double d01_x = p1.x - p0.x;
        double d01_y = p1.y - p0.y;
        double d12_x = p2.x - p1.x;
        double d12_y = p2.y - p1.y;
        double det = d01_x * d12_y - d12_x * d01_y;
        if (inequalities::is_close(det, 0)) return 0;
        return std::signbit(det) ? -1 : 1;
    };

    const mesh::Point& a = e0.org();
    const mesh::Point& b = e0.dest();
    const mesh::Point& c = e1.org();
    const mesh::Point& d = e1.dest();

    double d_abc = det_sign(a, b, c);
    double d_abd = det_sign(a, b, d);
    double d_cda = det_sign(c, d, a);
    double d_cdb = det_sign(c, d, b);

    // Test if edges are segments of the same line
    // NB: edges will never overlap for this case, only intersect at points
    if (d_abc == 0 && d_abd == 0)
        return a == c || a == d || b == c || b == d;

    // Test if an endpoint of one edge is on another edge
    if (d_abc == 0 || d_abd == 0)
        return d_cda != d_cdb;
    else if (d_cda == 0 || d_cdb == 0)
        return d_abc != d_abd;

    // Test if each edge splits endpoints of the other edge
    bool cd_on_opp_sides_of_ab = d_abc != d_abd;
    bool ab_on_opp_sides_of_cd = d_cda != d_cdb;
    return cd_on_opp_sides_of_ab && ab_on_opp_sides_of_cd;
}


} // namespace math

} // namespace umr
