#ifndef MATH_H
#define  MATH_H

#include <memory> // hash
#include <ostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "triangle.hpp"
#include "quadedge.hpp"


namespace umr {

namespace math {


/**
 * @brief Calculate a triangle's circumcenter.
 *
 * Calculates the center of a circle enscribing three given points.
 *
 * @param p0 A point.
 * @param p1 A point.
 * @param p2 A point.
 * @return The (x, y) coordinates of the circumcircle.
 */
std::pair<double, double> calculate_circumcenter(const mesh::Point& p0,
        const mesh::Point& p1,
        const mesh::Point& p2);


/**
 * @brief Calculate the determinant of three points.
 *
 * Calculates the determinant of three points of a triangle.
 *
 * @param p0 A point.
 * @param p1 A point.
 * @param p2 A point.
 * @return The value of the determinant.
 */
double triangle_determinant(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2);


/**
 * @brief Check if three points are in counter-clockwise order.
 *
 * Checks if the three points are counter-clockwise by calculating the
 * determinant. A positive determinant corresponds to
 * counter-clockwise ordering.
 *
 * @param p0 A point.
 * @param p1 A point.
 * @param p2 A point.
 * @return True if counter-clockwise, false otherwise.
 */
bool is_ccw(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2);


/**
 * @brief Check if a point is to the left of a line.
 *
 * Check if the point q is to the left of the line joining points p0
 * and p1.  The line is in the direction from point p0 to p1.
 *
 * @param q The point being tested.
 * @param p0 A point on the line.
 * @param p1 A point on the line.
 * @return True if q is left of the line, false otherwise.
 */
bool left_of(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1);

/**
 * @brief Check if a point is to the left of an edge.
 *
 * @param q The point being tested.
 * @param e An edge.
 * @return True if the point is left of the edge, false otherwise.
 */
bool left_of(const mesh::Point& q, const mesh::Edge& e);

/**
 * @brief Check if a point is to the right of a line.
 *
 * Check if the point q is to the right of the line joining points p0
 * and p1.  The line is in the direction from point p0 to p1.
 *
 * @param q The point being tested.
 * @param p0 A point on the line.
 * @param p1 A point on the line.
 * @return True if q is right of the line, false otherwise.
 */
bool right_of(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1);

/**
 * @brief Check if a point is to the right of an edge.
 *
 * @param q The point being tested.
 * @param e An edge.
 * @return True if the point is right of the edge, false otherwise.
 */
bool right_of(const mesh::Point& q, const mesh::Edge& e);

/**
 * @brief Check if a point is to on a line.
 *
 * Check if the point q is on the line joining points p0 and p1.
 * The line is in the direction from point p0 to p1.
 *
 * @param q The point being tested.
 * @param p0 A point on the line.
 * @param p1 A point on the line.
 * @return True if q is on the line, false otherwise.
 */
bool on_line(const mesh::Point& q, const mesh::Point& p0, const mesh::Point& p1);


/**
 * @brief Check if a point is on an edge.
 *
 * @param q The point being tested.
 * @param e An edge.
 * @return True if the point is on the edge, false otherwise.
 */
bool on_line(const mesh::Point& q, const mesh::Edge& e);

/**
 * @brief Check if a point is within a circle.
 *
 * Check if a point q is within a circle that enscribes three points
 * p0, p1, and p2.
 *
 * @param p0 A point on the circle.
 * @param p1 A point on the circle.
 * @param p2 A point on the circle.
 * @param q The point being tested.
 * @return True if the point is within the circle, false otherwise.
 */
bool in_circle(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& p2, const mesh::Point& q);

/**
 * @brief Check if a point is within a diametral circle.
 *
 * Check if a point q is within a diametral circle that enscribes two points
 * p0, and p1. The diametral circle is the smallest circle going through
 * these points.
 *
 * @param p0 A point on the circle.
 * @param p1 A point on the circle.
 * @param q The point being tested.
 * @return True if the points is within the circle, false otherwise.
 */
bool in_diametral_circle(const mesh::Point& p0, const mesh::Point& p1, const mesh::Point& q);

/**
 * @brief Check if two edges intersect.
 *
 * Checks if two edges intersect. Returns true when edges intersect
 * or share an endpoint. Returns false when edges do not intersect,
 * even if not parallel. If one edge is contained within another edge
 * it will return false unless they share an endpoint.
 *
 * @param e0 An edge.
 * @param e1 An edge.
 * @return True if edges intersect.
 */
bool intersect(const mesh::Edge& e0, const mesh::Edge& e1);

/**
 * @brief Law of cosines to calculate angle.
 *
 * Given edge lengths a, b, and c, this calculates the angle (in
 * radians) opposite c.
 *
 * @param a Length of side A
 * @param b Length of side B
 * @param c Length of side C
 * @return Angle opposite side C in radians.
 */
inline double loc_angle(double a, double b, double c) {
    double cos_angle_rad = (a*a + b*b - c*c) / 2 / a / b;
    double angle_rad = std::acos(cos_angle_rad);
    assert(angle_rad > 0 && "Angle must be positive!\n");
    return angle_rad;
}

/**
 * @brief Law of cosines to calculate angle.
 *
 * Given edge lengths a, b, and c, this calculates the angle (in
 * degrees) opposite c.
 *
 * @param a Length of side A
 * @param b Length of side B
 * @param c Length of side C
 * @return Angle opposite side C in degrees.
 */
inline double loc_angle_deg(double a, double b, double c) {
    static const double to_deg = 180 / std::numbers::pi;
    return loc_angle(a, b, c) * to_deg;
}


} // namespace math

} // namespace umr

#endif // MATH_H
