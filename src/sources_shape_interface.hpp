#ifndef SOURCES_SHAPE_INTERFACE_H
#define SOURCES_SHAPE_INTERFACE_H

#include <vector>
#include <memory>

#include "optimizers.hpp"
#include "parametrics.hpp"


namespace umr {

namespace source {


/**
 * @class IShape
 * @brief An interface to handle shape classes
 *
 * This class is designed to represent a set of IParametric objects as
 * a closed shape. It is capable of splitting segments uniformly into
 * subsegments.
 */
class IShape {
public:

    /** Constructor for IShape */
    IShape() {}

    /** Destructor for IShape */
    virtual ~IShape() {}

    /** Runs an optimizer to split IParametric objects */
    void refine();

    /** Set epsilon for the optimizer used when splitting */
    void set_optimizer_epsilon(double epsilon);

    /** Set the max number of iterations for the optimizer used when
        splitting */
    void set_optimizer_maxiter(size_t max_iter);

    /** Set a number of subsegments for each segment of the shape */
    virtual void set_number_subsegments(std::initializer_list<size_t> nseg);

    /**
     * Adds and IParametric pointer to the Shape
     *
     * @param ptr An IParametric object
     */
    void add_segment(std::shared_ptr<parametric::IParametric> ptr);

    /** Get a vector of subsegments representing the shape */
    std::vector<std::shared_ptr<parametric::IParametric>>& get_segments() { return m_segments; }

    /** Check if the shape is closed */
    bool is_closed() { return m_closed; }

protected:
    math::optimize::SplitSegmentUniformly m_optimizer;
    std::vector<std::shared_ptr<parametric::IParametric>> m_segments;

private:
    bool m_closed = false;
};


} // namespace source

} // namespace umr

#endif // SOURCES_SHAPE_INTERFACE_H
