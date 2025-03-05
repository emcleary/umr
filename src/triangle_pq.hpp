#ifndef TRIANGLE_PQ_H
#define TRIANGLE_PQ_H

#include <queue>
#include <unordered_map>
#include <cassert>

#include "quadedge.hpp"


namespace umr {

namespace mesh {

namespace algo {


/** Pairs a triangle with one of its edges */
using TriangleEdge = std::pair<Triangle*, Edge*>;


/**
 * @class ITrianglePQ
 * @brief A priority queue for triangles
 *
 * This is an interface for a priority queue for triangles.  It is
 * required for the refinement algorithms. Any triangle removed from
 * the mesh will get deactivated in the queue and deleted from memory
 * when popped.
 *
 * The interface allows for abstracting TrianglePQ based on its
 * priority type.
 */
class ITrianglePQ {
public:

    /** Constructor for ITrianglePQ */
    ITrianglePQ() {}

    /** Destructor for ITrianglePQ */
    virtual ~ITrianglePQ() {}

    /**
     * Pushes a triangle into the queue
     *
     * @param t A triangle
     * @param e Any Edge belonging to that triangle
     */
    virtual void push(Triangle* t, Edge* e) = 0;

    /**
     * Accesses the top active element of the queue,
     * deleting all inactive triangles before it.
     *
     * @param t A triangle
     * @param e Any Edge belonging to that triangle
     * @return Returns an active triangle paired with its edge.
     */
    virtual TriangleEdge top() = 0;

    /**
     * Removes the top active triangle from the queue,
     * deleting all inactive triangles before it.
     */
    virtual void pop() = 0;

    /**
     * Clears the queue, deleting all inactive triangles
     * in the process.
     */
    virtual void clear() = 0;

    /**
     * Gets the size of the priority queue, including both
     * active and inactive triangles.
     *
     * @return Size of the queue
     */
    size_t size();

    /**
     * Deactivates the given triangle in the queue.
     *
     * @param t A triangle
     */
    void deactivate(Triangle* t);

    /**
     * Checks if a given triangle is in the queue.
     *
     * @param t A triangle
     * @return True if the triangle is in the queue, false otherwise
     */
    bool contains(Triangle* t);

    /**
     * Checks if a triangle is active in the queue.
     *
     * @param t A triangle
     * @return True if the triangle is active in the queue, false otherwise
     */
    bool is_active(Triangle* t);

    /**
     * Checks if a triangle is inactive in the queue.
     *
     * @param t A triangle
     * @return True if the triangle is not active and in the queue,
     * false otherwise
     */
    bool is_inactive(Triangle* t);

    /**
     * Checks if the queue has any active triangles by running top.
     *
     * @return True if the queue has an active element, false otherwise
     */
    bool empty();

protected:
    std::unordered_map<Triangle*, bool> m_active;
};


/**
 * @class TrianglePQ
 * @brief A template class for triangle priority queue
 *
 * This class is a priority queue for triangles with a template
 * for comparison operators for the queue.
 *
 * @tparam Compare Comparison struct used in the priority queue
 */
template<class Compare>
class TrianglePQ : public ITrianglePQ {
public:

    /** Constructor for the TrianglePQ */
    TrianglePQ() {}

    /** Destructor for the TrianglePQ, deleting inactive triangles */
    virtual ~TrianglePQ() override {
        while (!m_pq.empty()) {
            auto te = m_pq.top();
            m_pq.pop();
            if (!m_active[te.first])
                delete te.first;
        }
    }

    /** Add a triangle to the priority queue */
    virtual void push(Triangle* t, Edge* e) {
        assert(!m_active.contains(t));
        m_active[t] = true;
        m_pq.push({t, e});
    }

    /** Returns the top triangle from the priority queue */
    virtual TriangleEdge top() {
        while (!m_pq.empty()) {
            TriangleEdge te = m_pq.top();
            if (m_active[te.first])
                return te;

            // triangle no longer active so delete!
            m_pq.pop();
            m_active.erase(te.first);
            delete te.first;
        }

        return {nullptr, nullptr};
    }

    /**
     * Removes the top active triangle in the queue, deleting
     * inactive triangles in the process
     */
    virtual void pop() {
        Triangle* t = top().first;
        if (t) {
            m_pq.pop();
            m_active.erase(t);
        }
    }

    /**
     * Clears the priority queue, deleting inactive triangles
     * in the process.
     */
    virtual void clear() {
        while (!m_pq.empty()) {
            auto te = m_pq.top();
            m_pq.pop();
            if (!m_active[te.first])
                delete te.first;
            m_active.erase(te.first);
        }
    }

private:
    std::priority_queue<TriangleEdge, std::vector<TriangleEdge>, Compare> m_pq;
};


/**
 * Comparison operators for the TrianglePQ
 */
struct TriCompMaxSize {
    bool operator() (TriangleEdge& a, TriangleEdge& b);
};


struct TriCompMinSize {
    bool operator() (TriangleEdge& a, TriangleEdge& b);
};


struct TriCompMaxAngleMin {
    bool operator() (TriangleEdge& a, TriangleEdge& b);
};


struct TriCompMinAngleMin {
    bool operator() (TriangleEdge& a, TriangleEdge& b);
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // TRIANGLE_PQ_H
