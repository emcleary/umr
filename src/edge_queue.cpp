#include "edge_queue.hpp"


namespace umr {

namespace mesh {

namespace algo {


EdgeQueue::EdgeQueue(EdgeQueue& queue)
        : m_queue(queue.m_queue.begin(), queue.m_queue.end()) {}


EdgeQueue::EdgeQueue(EdgeQueue&& queue)
        : m_queue(queue.m_queue.begin(), queue.m_queue.end()) {}


void EdgeQueue::operator=(EdgeQueue& queue) {
    for (const Edge& e : queue.m_queue)
        insert(e);
}


void EdgeQueue::operator=(EdgeQueue&& queue) {
    for (const Edge& e : queue.m_queue)
        insert(e);
}


bool EdgeQueue::insert(const Edge& e) {
    if (m_queue.contains(e.sym()))
        return false;
    return m_queue.insert(e).second;
}


Edge EdgeQueue::top() {
    return *(std::prev(m_queue.end()));
}


void EdgeQueue::pop() {
    auto it = std::prev(m_queue.end());
    m_queue.erase(it);
}


bool EdgeQueue::empty() const {
    return m_queue.empty();
}


size_t EdgeQueue::size() const {
    return m_queue.size();
}


} // namespace algo

} // namespace mesh

} // namespace umr
