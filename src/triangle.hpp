#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <ostream>

#include "point.hpp"


namespace umr {

namespace mesh {


/**
 * @class Triangle
 * @brief A triangle cell in the mesh.
 *
 * A triangle cell in the mesh. Points used are taken as references to
 * points in the mesh.
 */
class Triangle {
    friend std::ostream& operator<<(std::ostream& out, const Triangle& t);
public:

    /**
     * @brief Constructs a triangle with three points.
     *
     * @param p0 A point in the triangle.
     * @param p1 A point in the triangle.
     * @param p2 A point in the triangle.
     */
    Triangle(const Point& p0, const Point& p1, const Point& p2);

    ~Triangle() {}

    /** Get the radius of the triangle's circumcircle. */
    const double get_radius() const { return m_radius; }

    /** Get the circumcenter of the triangle's circumcircle. */
    const Point& get_circumcenter() const { return m_center; }

    /** Get the minimum angle of the triangle. */
    const double get_min_angle() const { return m_min_angle; }

    /** Get the maximum angle of the triangle. */
    const double get_max_angle() const { return m_max_angle; }

    /** Get the minimum side length of the triangle. */
    const double get_min_length() const { return m_min_length; }

    /** Get the maximum side length of the triangle. */
    const double get_max_length() const { return m_max_length; }

    /** Get the pointers to the triangle's points. */
    std::vector<const Point*> get_points() const { return {&m_p0, &m_p1, &m_p2}; }

private:
    const Point& m_p0;
    const Point& m_p1;
    const Point& m_p2;
    Point m_center;
    double m_radius;
    double m_min_angle;
    double m_max_angle;
    double m_min_length;
    double m_max_length;
};


/**
 * @brief Calculate a triangle's circumcenter.
 *
 * Calculates the center of a circle enscribing three given points.
 *
 * @param p0 A point.
 * @param p1 A point.
 * @param p2 A point.
 * @return A Point object.
 */
Point make_circumcenter(const Point& p0, const Point& p1, const Point& p2);


inline std::ostream& operator<<(std::ostream& out, const Triangle& t) {
    out << "Triangle: " << t.m_p0 << " " << t.m_p1 << " " << t.m_p2 << '\n';
    out << "   " << "Circumcenter: " << t.m_center << '\n';
    out << "   " << "Radius: " << +t.m_radius << '\n';
    out << "   " << "Angle min = " << +t.m_min_angle << ", Angle max = " << +t.m_max_angle;
    return out;
}


} // mesh

} // umr

#endif // TRIANGLE_H
