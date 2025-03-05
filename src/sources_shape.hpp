#ifndef SOURCES_SHAPE_H
#define SOURCES_SHAPE_H

#include "sources_shape_interface.hpp"

#include "point.hpp"


namespace umr {

namespace source {


/**
 * @class Shape
 * @brief A class to handle general shapes
 *
 * This class is designed to handle general shapes input by the user.
 * The class contains methods for adding Line and Arc parameterics, as
 * well any ParametricFunction. The initial number of subsegments of
 * each segment can be specified, which is useful for controlling
 * triangle density along those segments. For ParametricFunctions, it
 * is especially important to specify enough subsegments to resolve
 * all local features of the function. This facilitates any
 * optimization methods for splitting segments, calculating distances
 * between sources, etc.
 */
class Shape : public IShape {
public:

    /** Constructor for Shape */
    Shape() {}

    /** Destructor for Shape */
    ~Shape() {}

    /**
     * Constructs and adds a Line to the Shape
     *
     * @param p0 A Point
     * @param p1 A Point
     * @param nseg Number of subsegments of the line
     */
    void add_line(mesh::Point& p0, mesh::Point& p1, size_t nseg = 1);

    /**
     * Constructs and adds an Arc to the Shape
     *
     * @param pc The circle's center
     * @param radius The circle's radius
     * @param tmin Lower bound for the Arc
     * @param tmax Upper bound for the Arc
     * @param nseg Number of subsegments of the line
     */
    void add_arc(mesh::Point& pc, double radius, double tmin, double tmax, size_t nseg = 3);

    /**
     * Constructs and adds a Parametric to the Shape
     *
     * @param function The Parametric function
     * @param tmin Lower bound for the Arc
     * @param tmax Upper bound for the Arc
     * @param nseg Number of subsegments of the line
     */
    void add_parametric(parametric::ParametricFunction function, double tmin, double tmax, size_t nseg = 3);
};


} // namespace source

} // namespace umr

#endif // SOURCES_SHAPE_H
