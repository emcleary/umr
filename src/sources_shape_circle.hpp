#ifndef SOURCES_SHAPE_CIRCLE_H
#define SOURCES_SHAPE_CIRCLE_H

#include "sources_shape_interface.hpp"

#include "point.hpp"


namespace umr {

namespace source {


/**
 * @class ShapeQuadrilateral
 * @brief A class to handle circle
 */
class ShapeCircle : public IShape {
public:

    /**
     * Constructor for Circle
     *
     * @param pc Center Point
     * @param radius Radius of the circle
     */
    ShapeCircle(mesh::Point& pc, double radius);

    /** Destructor fo the ShapeCircle */
    ~ShapeCircle() {}

    /**
     * Set a number of subsegments for the circle's segment. The list
     * must have a size of 1.
     */
    virtual void set_number_subsegments(std::initializer_list<size_t> nseg);

};


} // namespace source

} // namespace umr

#endif // SOURCES_SHAPE_CIRCLE_H
