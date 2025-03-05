#include "sources_manager.hpp"

#include <memory>

#include "inequalities.hpp"


namespace umr {

namespace source {


void SourceManager::add_shape(std::shared_ptr<IShape> shape) {
    if (!shape->is_closed()) {
        std::cout << "Shapes must be fully closed before adding them "
            "to the SourceManager.\n";
        exit(-1);
    }

    for (std::shared_ptr<parametric::IParametric> s : shape->get_segments())
        m_segments.push_back(s);

    m_new_segment = true;
    ++m_num_shapes;
}


void SourceManager::add_segment(std::shared_ptr<parametric::IParametric> segment) {
    m_segments.push_back(segment);
    m_new_segment = true;
}


void SourceManager::add_point(double x, double y) {
    m_points.push_back({x, y});
    m_new_point = true;
}


void SourceManager::add_point(mesh::Point& p) {
    m_points.push_back({p.x, p.y});
    m_new_point = true;
}


void SourceManager::finalize() {
    fix_point_segment_precision();
    split_segments();
    m_new_point = false;
    m_new_segment = false;
}


std::list<parametric::IParametric*> SourceManager::get_segments() {
    std::list<parametric::IParametric*> segments;
    for (auto segment : m_segments)
        segments.push_back(segment.get()->clone());
    return segments;
}


void SourceManager::fix_point_segment_precision() {
    if (m_new_point || m_new_segment) {
        for (auto [xp, yp] : m_points) {
            auto segment_iter = m_segments.begin();
            while (segment_iter != m_segments.end()) {
                auto [xmin, ymin] = (*segment_iter)->evaluate_tmin();
                if (inequalities::is_close(xp, xmin)
                        && inequalities::is_close(yp, ymin))
                    if (xp != xmin || yp != ymin)
                        (*segment_iter)->set_initial_point(xp, yp);
                auto [xmax, ymax] = (*segment_iter)->evaluate_tmax();
                if (inequalities::is_close(xp, xmax)
                        && inequalities::is_close(yp, ymax))
                    if (xp != xmax || yp != ymax)
                        (*segment_iter)->set_final_point(xp, yp);
                ++segment_iter;
            }
        }
    }
}


void SourceManager::split_segments() {
    if (m_new_segment) {
        auto segment_iter = m_segments.begin();
        while (segment_iter != m_segments.end()) {
            if ((*segment_iter)->get_num_subsegments() > 1) {
                (*segment_iter)->optimize_parameters();
                std::vector<std::shared_ptr<parametric::IParametric>>
                    subsegments = (*segment_iter)->get_subsegments();
                ++segment_iter;
                m_segments.erase(std::prev(segment_iter));
                m_segments.insert(segment_iter,
                        subsegments.begin(), subsegments.end());
            } else {
                ++segment_iter;
            }
        }
    }
}


} // namespace source

} // namespace umr
