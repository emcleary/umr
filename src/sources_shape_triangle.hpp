#ifndef SOURCES_SHAPE_TRIANGLE_H
#define SOURCES_SHAPE_TRIANGLE_H

#include "sources_shape_interface.hpp"

#include "point.hpp"


namespace umr {

namespace source {

/**
 * @class ShapeTriangle
 * @brief A class to handle triangles
 */
class ShapeTriangle : public IShape {
public:

    /**
     * Constructor for ShapeTriangle. Input points are assumed to
     * be in clockwise or counter-clockwise order, each pair of points
     * contructing a line.
     *
     * @param p0 A Point
     * @param p1 A Point
     * @param p2 A Point
     */
    ShapeTriangle(mesh::Point& p0, mesh::Point& p1, mesh::Point& p2);

    /** Destructor for ShapeTriangle */
    ~ShapeTriangle() {}

    /**
     * Set a number of subsegments for each segment. The input list
     * must be of size 3. Elements of the list correspond to lines
     * joining points in the order this object was constructed.
     */
    virtual void set_number_subsegments(std::initializer_list<size_t> nseg);

private:
    /** Confirms no input points are duplicates and that they are not all
        on a single line */
    void validate_construction(mesh::Point& p0, mesh::Point& p1, mesh::Point& p2);
};


} // namespace source

} // namespace umr

#endif // SOURCES_SHAPE_TRIANGLE_H
