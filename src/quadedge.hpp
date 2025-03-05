#ifndef QUADEDGE_H
#define QUADEDGE_H

#include <memory>
#include <ostream>

#include "point.hpp"
#include "triangle.hpp"


namespace umr {

namespace mesh {

template<class A, class B> class QuadEdge;
class Edge;
class Dual;


/**
 * This namespaces is a set of tools for constructing, destructing,
 * and manipulating Edges.
 */
namespace quadedge {

///////////////////////////////////////////////////
// Constructors and destructor methods for Edges //
///////////////////////////////////////////////////

/**
 * The Edge constructor, creating Edges and their Duals joining two
 * points.
 *
 * @param p0 A Point
 * @param p1 A Point
 * @return An Edge created with p0 as its origin
 */
Edge& make_edge(Point& p0, Point& p1);

/**
 * An Edge constructor, creating Edges and their Duals by joining an
 * Edge destination with a Point. The Point must on the line of the
 * Edge beyond the destination, or anywhere to the left of the Edge.
 *
 * @param e0 The Edge
 * @param p The Point
 * @return An Edge created with e0's destination as its origin
 */
Edge& extend_edge(Edge& e0, Point& p);

/**
 * An Edge constructor, creating Edges and their Duals by joining two
 * Edges. The destination of Edge e0 is joined with the origin of Edge
 * e1.
 *
 * @param e0 An Edge
 * @param e1 An Edge
 * @return An Edge created with e0's destination as its origin
 */
Edge& connect(Edge& e0, Edge& e1);

/**
 * Construct a Triangle object belonging to the cell left of the given
 * edge. If the edge is not triangular, no Triangle object is created
 * and the method returns a nullptr.
 *
 * @param e The given Edge
 * @return The new Triangle or a nullptr
 */
Triangle* make_triangle(Edge& e);

/**
 * Delete Triangle cell objects corresponding to and edge, if
 * present.
 *
 * @param e The input Edge
 */
void delete_cell(Edge* e);

/**
 * Detaches an Edge from its neighboring edges. This function does not
 * free any memory used.
 *
 * @param e The Edge to detach
 */
void remove_edge(Edge* e);

/**
 * Delete an Edge, its symmetric Edge, and its Dual edges.
 *
 * @param e The Edge to delete
 */
void delete_edge(Edge* e);

/**
 * The Edge destructor, detaching an Edge from neighboring Edges and
 * freeing all memory belonging to the Edge: the Edges, Duals, and
 * Triangles.
 *
 * @param e The input Edge
 */
void delete_edge_all(Edge* e);


/////////////////////////////////////////////
// Manipulates Edges among their neighbors //
/////////////////////////////////////////////

/**
 * This method swaps onext pointers of the given Edges.  It is useful
 * for inserting and removed Edges from a Mesh, enabling Edge
 * rotations to be accurate. *
 *
 * @param e0 An Edge
 * @param e1 An Edge
 */
void splice(Edge& e0, Edge& e1);

/**
 * Swaps an Edge separating two Triangles.
 *
 * @param e An Edge
 */
void swap(Edge& e);


///////////////////
// Other methods //
///////////////////

/**
 * Checks if an Edge has a triangle-shaped cell on its left.
 *
 * @param e0 An Edge
 * @result True if the cell on the left is a triangle, false otherwise
 */
bool part_of_triangle(const Edge& e0);

/**
 * Checks if an Edge is attached to any other Edges.
 *
 * @param e An Edge
 * @result True if free (not attached), false otherwise
 */
bool is_free(const Edge& e);


} // namespace quadedge


/**
 * @class QuadEdge
 * @brief A class representing the QuadEdge datastructure
 *
 * This class represents the QuadEdge datastructure, handling math
 * operations for an Edge to move about its origin and destination,
 * and to move to Edges around its neighboring cells. It is adapted
 * from "Guibas and Stolfi, Primitives for the manipulation of general
 * subdivisions and the computation of Voronoi diagrams, ACM
 * Transactions on Graphics, Vol. 4, No. 2, 1985, 74-123."
 *
 * Template parameters must be (Edge, Dual) or (Dual Edge).
 *
 * @tparam A Either Edge or Dual
 * @tparam B Either Dual or Edge
 */
template<class A, class B>
class QuadEdge {
    friend Edge& quadedge::make_edge(const Point& p0, const Point& p1);

protected:
    /** Constructor for QuadEdge */
    QuadEdge() {}

public:

    /** \defgroup Topological Operators
     *
     * @brief Topological operators
     *
     * Topological operators used to move between Edges and Duals
     *
     * @{
     */

    /** @brief Rotate counter-clockwise */
    B& rot() const { return *m_right; }

    /** @brief Rotate clockwise */
    B& rot_inv() const { return this->left(); }

    /** @brief Get Edge/Dual to the right */
    B& right() const { return this->rot(); }

    /** @brief Get Edge/Dual to the left */
    B& left() const { return this->rot().rot().rot(); }

    /** @brief Get the symmetric Edge/Dual */
    A& sym() const { return this->rot().rot(); }

    /** @brief Rotate about the origin counter-clockwise */
    A& onext() const { return *m_onext; }

    /** @brief Rotate about the origin clockwise */
    A& oprev() const { return this->rot().onext().rot(); }

    /** @brief Rotate about the destination counter-clockwise */
    A& dnext() const { return this->sym().onext().sym(); }

    /** @brief Rotate about the destination clockwise */
    A& dprev() const { return this->rot_inv().onext().rot_inv(); }

    /** @brief Rotate around the left cell counter-clockwise */
    A& lnext() const { return this->rot_inv().onext().rot(); }

    /** @brief Rotate around the left cell clockwise */
    A& lprev() const { return this->onext().sym(); }

    /** @brief Rotate around the right cell counter-clockwise */
    A& rnext() const { return this->rot().onext().rot_inv(); }

    /** @brief Rotate around the right cell clockwise */
    A& rprev() const { return this->sym().onext(); }

    /** @} */

protected:

    A* m_onext;
    B* m_right;
};

/**
 * @class Edge
 * @brief An Edge belonging to the QuadEdge datastructure
 *
 * This class represents the Edges of the QuadEdge datastructure.  Its
 * constructors are private and only available through functions in
 * the quadedge namespace to ensure corresponding Edges and Duals are
 * all constructed together and properly linked.
 */
class Edge : public QuadEdge<Edge, Dual> {
    friend Edge& quadedge::make_edge(Point& p0, Point& p1);
    friend void quadedge::splice(Edge& e0, Edge& e1);
    friend void quadedge::swap(Edge& e);
    friend std::ostream& operator<<(std::ostream& out, const Edge& e);

public:
    /** Get origin Point */
    Point& org() const { return *m_org; }

    /** Get destination Point */
    Point& dest() const { return *m_dest; }

    bool operator==(const Edge& e) const;

    bool operator!=(const Edge& e) const;

    bool operator<(const Edge& e) const;

private:
    Edge() {}
    Edge(Point& p0, Point& p1) : m_org(&p0), m_dest(&p1) {}

    Point* m_org;
    Point* m_dest;
};

/**
 * @class Dual
 * @brief A Dual belonging to the QuadEdge datastructure
 *
 * This class represents the Duals of the QuadEdge datastructure.  Its
 * constructors are private and only available through functions in
 * the quadedge namespace to ensure corresponding Edges and Duals are
 * all constructed together and properly linked.
 *
 * Its m_data parameter represents a cell object. Currently only the
 * Triangle cell type is implemented.
 */
class Dual : public QuadEdge<Dual, Edge> {
    friend Edge& quadedge::make_edge(Point& p0, Point& p1);
    friend void quadedge::splice(Edge& e0, Edge& e1);
    friend void quadedge::delete_cell(Edge* e);
    friend void quadedge::delete_edge(Edge* e);
    friend void quadedge::remove_edge(Edge* e);
    friend std::ostream& operator<<(std::ostream& out, const Dual& d);
public:
    /** Get data of the cell */
    Triangle* get_data() { return m_data; }

    /** Sets data with a Triangle */
    void set_data(Triangle* cell);

    /** Clear data; no memory is freed */
    void clear_data() { m_data = nullptr; }

private:
    /** Constructor for Dual */
    Dual(Point& p0, Point& p1) : m_org(&p0), m_dest(&p1) {}

    /** Destructor for Dual */
    ~Dual();

    Triangle* m_data = nullptr;
    Point* m_org;
    Point* m_dest;
};


std::ostream& operator<<(std::ostream& out, const Edge& e);

std::ostream& operator<<(std::ostream& out, const Dual& d);


} // namespace mesh

} // namespace umr



template<>
struct std::hash<umr::mesh::Edge> {

    /** A hash calculation of an Edge */
    size_t operator()(const umr::mesh::Edge& e) {
        size_t h = 17;
        h += h * 31 + std::hash<double>()(e.org().x);
        h += h * 31 + std::hash<double>()(e.org().y);
        h += h * 31 + std::hash<double>()(e.dest().x);
        h += h * 31 + std::hash<double>()(e.dest().y);
        return h;
    }

};

#endif // QUADEDGE_H
