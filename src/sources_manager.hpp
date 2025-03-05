#ifndef SOURCES_MANAGER_H
#define SOURCES_MANAGER_H

#include "sources_shape_interface.hpp"

#include <memory>

#include "point.hpp"
#include "parametrics.hpp"


namespace umr {

namespace source {


/**
 * @class SourceManager
 * @brief Manages a collection of IShape objects
 *
 * This class manages a collection of IShape object. An instance of
 * this class can be passed to the Builder and is used to generate a
 * set of input segments for building the mesh.
 *
 * The user is responsible to ensure all shapes and segments do not
 * intersect.
 */
class SourceManager {
public:

    /** Constructor for the SourceManager */
    SourceManager() {}

    /** Destructor for the SourceManager */
    ~SourceManager() {}

    /**
     * Adds an IShape to the SourceManager. The user is responsible
     * for ensuring only one shape acts as the exterior and all other
     * shapes added (if any) do not intersect with any other shape.
     */
    void add_shape(std::shared_ptr<IShape> shape);

    /**
     * Adds an IParametric segment to the SourceManager. It is the user's
     * responsibility to ensure each segment added does not intersect
     * with any other segment or shape.
     */
    void add_segment(std::shared_ptr<parametric::IParametric> segment);

    /**
     * Adds a point to the SourceManager which represents a source
     * point in the mesh. The user is responsible for ensure the point
     * does not overlap with any input point, segment or shape.
     */
    void add_point(double x, double y);

    /**
     * Adds a point to the SourceManager which represents a source
     * point in the mesh. The user is responsible for ensure the point
     * does not overlap with any input point, segment or shape.
     */
    void add_point(mesh::Point& p);

    /** Split all input segments and shapes as needed. */
    void finalize();

    /** Returns a copy of all points added. */
    std::vector<Real2> get_points() { return m_points; }

    /** Returns a list of all (split) segments (including those from shapes) */
    std::list<parametric::IParametric*> get_segments();

    /** Gets the number of shapes added to the ShapeManager. */
    unsigned int get_num_shapes() { return m_num_shapes; }

private:
    /** Ensures no precision errors between segments and points */
    void fix_point_segment_precision();

    /** Splits segments as needed, inserting subsegments into the list */
    void split_segments();

    std::list<std::shared_ptr<parametric::IParametric>> m_segments;
    std::vector<Real2> m_points;
    bool m_new_point = false;
    bool m_new_segment = false;
    unsigned int m_num_shapes = 0;
};


} // namespace source

} // namespace umr

#endif // SOURCES_MANAGER_H
