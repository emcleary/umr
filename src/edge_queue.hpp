#ifndef EDGE_QUEUE_H
#define EDGE_QUEUE_H

#include <set>
#include "quadedge.hpp"


namespace umr {

namespace mesh {

namespace algo {

/**
 * @class EdgeQueue
 * @brief A queue for Edge objects
 *
 * This is a queue for Edges, designed to skip inserting edges that
 * are symmetric to any edge already in the queue. The class stores
 * edges in a set, so no ordering is enforced.
 */
class EdgeQueue {
public:

    /** Constructor for EdgeQueue */
    EdgeQueue() {}

    /** Constructor for EdgeQueue */
    EdgeQueue(EdgeQueue& queue);

    /** Constructor for EdgeQueue */
    EdgeQueue(EdgeQueue&& queue);

    void operator=(EdgeQueue& queue);

    void operator=(EdgeQueue&& queue);

    /**
     * Inserts and Edge into the queue. The method returns false if
     * either the edge or its symmetric edge exists in the queue.
     *
     * @param e An Edge to be inserted
     * @return True if the edge is inserted, false otherwise
     */
    bool insert(const Edge& e);

    /**
     * Get the top element of the queue.
     *
     * @return The top Edge from the queue
     */
    Edge top();

    /**
     * Removes the top edge from the queue.
     */
    void pop();

    /**
     * Checks if the queue is empty.
     *
     * @return True if the queue is empty, false otherwise
     */
    bool empty() const;

    /**
     * Gets the size of the queue.
     *
     * @return The size of the queue
     */
    size_t size() const;

private:
    std::set<Edge> m_queue;
};


} // namespace algo

} // namespace mesh

} // namespace umr


#endif // EDGE_QUEUE_H
