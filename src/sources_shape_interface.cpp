#include "sources_shape_interface.hpp"

#include "inequalities.hpp"


namespace umr {

namespace source {


void IShape::set_number_subsegments(std::initializer_list<size_t> nseg) {
    assert(nseg.size() == m_segments.size());
    auto segment_iter = m_segments.begin();
    auto nseg_iter = nseg.begin();
    while (segment_iter != m_segments.end()) {
        (*segment_iter)->set_num_subsegments(*nseg_iter);
        ++segment_iter;
        ++nseg_iter;
    }
}

void IShape::refine() {
    for (std::shared_ptr<parametric::IParametric> ptr : m_segments) {
        std::vector<double> param = m_optimizer.run(ptr.get());
        ptr->set_parameters(param);
    }
}

void IShape::set_optimizer_epsilon(double epsilon) {
    m_optimizer.set_epsilon(epsilon);
}

void IShape::set_optimizer_maxiter(size_t max_iter) {
    m_optimizer.set_num_iter(max_iter);
}


/*
 * Here consectuve segments are joined together. It would be more
 * robust to, at a later place in the code, just loop over all
 * segments, compare all pairs of start and final points, adjusting
 * any points as needed. It would be more expensive, but it might be
 * easier for users since they wouldn't have to keep a close eye on
 * points themselves if multiple segments intersect.
 */
void IShape::add_segment(std::shared_ptr<parametric::IParametric> ptr) {
    if (is_closed()) {
        std::cout << "Shape is already closed. No more segments can be added.\n";
        return;
    }

    // Must ensure consecutive elements share the EXACT same point
    auto join_points = [](double x0, double y0, double x1, double y1,
            std::shared_ptr<parametric::IParametric> ptr) {
        if (x0 != x1 || y0 != y1) {
            if (inequalities::is_close(x0, x1) && inequalities::is_close(y0, y1)) {
                ptr->set_initial_point(x0, y0);
            } else {
                std::cout << "Shape::add_segment: consecutive segments "
                    "must end and start at the exact same point!\n";
                std::cout << std::format("Previous endpoint: x={}, y={}\n", x0, y0);
                std::cout << std::format("Current startpoint: x={}, y={}\n", x1, y1);
                exit(-1);
            }
        }
    };

    if (m_segments.size() > 0) {
        auto [x0, y0] = m_segments.back()->evaluate_tmax();
        auto [x1, y1] = ptr->evaluate_tmin();
        join_points(x0, y0, x1, y1, ptr);
    }

    m_segments.push_back(ptr);
    auto [x0, y0] = m_segments.back()->evaluate_tmax();
    auto [x1, y1] = m_segments.front()->evaluate_tmin();
    m_closed = inequalities::is_close(x0, x1) && inequalities::is_close(y0, y1);
    if (m_closed)
        join_points(x0, y0, x1, y1, m_segments.front());
}


} // namespace source

} // namespace umr
