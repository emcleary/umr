#ifndef POINT_H
#define POINT_H

#include <memory>
#include <vector>
#include <ostream>


namespace umr {

namespace mesh {


/**
 * @class Point
 * @brief Point in space
 *
 * A 2-dimension point in space.
 */
struct Point {
    const double x;
    const double y;

    /**
     * @brief Creates a point at (0, 0).
     */
    Point();

    /**
     * @brief Creates a point at coordinates (x, y).
     *
     * @param x The x-coordinate.
     * @param y The y-coordinate.
     */
    Point(double x, double y);

    /**
     * @brief Creates a point at coordinates (x, y).
     *
     * @param p A pair of coordinates (x, y).
     */
    Point(std::pair<double, double> p);

    Point(const Point& p);
    Point(const Point&& p);

    Point& operator=(const Point& p) = delete;
    Point& operator=(const Point&& p) = delete;

    ~Point() {};

    /**
     * @brief Calculate the distance to a point.
     *
     * Calculates the Euclidian distance between the
     * object point and another point.
     *
     * @param p A point
     * @return The distance between points
     */
    double distance_to(const Point& p) const;

    /**
     * @brief Calculate the distance to a point.
     *
     * Calculates the Euclidian distance between the
     * object point and another point.
     *
     * @param x x-coordinate of a point
     * @param y y-coordinate of a point
     * @return distance between points
     */
    double distance_to(const double x, const double y) const;

    bool operator==(const Point& p) const;

    bool operator<(const Point& p) const;
};


std::ostream& operator<<(std::ostream& out, const Point& p);


struct PointHash {
    size_t operator()(const Point& p) const;
};


} // namespace mesh

} // namespace umr


template<>
struct std::hash<umr::mesh::Point> {
    size_t operator()(const umr::mesh::Point& p);
};


#endif // POINT_H
