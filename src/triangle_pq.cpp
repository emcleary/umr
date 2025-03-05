#include "triangle_pq.hpp"

#include <queue>

#include "inequalities.hpp"


namespace umr {

namespace mesh {

namespace algo {


void ITrianglePQ::deactivate(Triangle* t) {
    // Only deactivate triangles still in the queue!
    // Triangle can be removed already from get_next.
    if (m_active.contains(t))
        m_active[t] = false;
    else
        delete t;
}


size_t ITrianglePQ::size() {
    return m_active.size();
}


bool ITrianglePQ::contains(Triangle* t) {
    return m_active.contains(t);
}


bool ITrianglePQ::is_active(Triangle* t) {
    if (m_active.contains(t))
        return m_active[t] == true;
    return false;
}


bool ITrianglePQ::is_inactive(Triangle* t) {
    if (m_active.contains(t))
        return m_active[t] == false;
    return false;
}


bool ITrianglePQ::empty() {
    return top().first == nullptr;
}


bool TriCompMaxSize::operator() (TriangleEdge& a, TriangleEdge& b) {
    Triangle* ta = a.first;
    Triangle* tb = b.first;
    return ta->get_radius() < tb->get_radius();
}


bool TriCompMinSize::operator() (TriangleEdge& a, TriangleEdge& b) {
    Triangle* ta = a.first;
    Triangle* tb = b.first;
    return ta->get_radius() > tb->get_radius();
}


bool TriCompMaxAngleMin::operator() (TriangleEdge& a, TriangleEdge& b) {
    Triangle* ta = a.first;
    Triangle* tb = b.first;
    return ta->get_min_angle() < tb->get_min_angle();
}


bool TriCompMinAngleMin::operator() (TriangleEdge& a, TriangleEdge& b) {
    Triangle* ta = a.first;
    Triangle* tb = b.first;
    double aa = ta->get_min_angle();
    double ab = tb->get_min_angle();
    if (inequalities::is_close(aa, ab))
        return inequalities::is_lt(ta->get_radius(), tb->get_radius());
    return inequalities::is_gt(aa, ab);
}


} // namespace algo

} // namespace mesh

} // namespace umr
