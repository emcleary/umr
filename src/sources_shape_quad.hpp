#ifndef SOURCES_SHAPE_QUAD_H
#define SOURCES_SHAPE_QUAD_H

#include "sources_shape_interface.hpp"

#include "point.hpp"


namespace umr {

namespace source {

/**
 * @class ShapeQuadrilateral
 * @brief A class to handle quadrilateral shapes
 */
class ShapeQuadrilateral : public IShape {
public:

    /**
     * Constructor for ShapeQuadrilateral. Input points are assumed to
     * be in clockwise or counter-clockwise order, each pair of points
     * contructing a line.
     *
     * @param p0 A Point
     * @param p1 A Point
     * @param p2 A Point
     * @param p3 A Point
     */
    ShapeQuadrilateral(mesh::Point& p0, mesh::Point& p1,
            mesh::Point& p2, mesh::Point& p3);

    /** Destructor for ShapeQuadrilateral */
    ~ShapeQuadrilateral() {}

    /**
     * Set a number of subsegments for each segment. The input list
     * must be of size 4. Elements of the list correspond to lines
     * joining points in the order this object was constructed.
     */
    virtual void set_number_subsegments(std::initializer_list<size_t> nseg);

private:
    /** Confirms no input points are duplicates, and no lines cross */
    void validate_construction(mesh::Point& p0, mesh::Point& p1,
            mesh::Point& p2, mesh::Point& p3);
};


} // namespace source

} // namespace umr

#endif // SOURCES_SHAPE_QUAD_H
