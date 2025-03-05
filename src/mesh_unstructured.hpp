#ifndef MESH_UNSTRUCTURED_H
#define MESH_UNSTRUCTURED_H

#include <unordered_set>

#include "point.hpp"
#include "quadedge.hpp"
#include "utilities.hpp"


namespace umr {

namespace mesh {


/**
 * @class MeshUnstructured
 * @brief A class for managing Point and Edges in a Mesh
 *
 * This class is a tool for building a mesh. It can be seen as a
 * database contains a set of points and edges, a set of contructors
 * and destructors for these objects, and a set of tools for joining
 * edges to build the mesh.
 *
 * The constructors and destructors for Points and Edges can be seen
 * as wrappers for those methods of the quadedge namespace. These
 * implementations are specifically designed to prevent any duplicates
 * and memory leaks.
 */
class MeshUnstructured {

public:
    /** Stores a unique set of Point pointers */
    using PointSet = std::unordered_set<Point*, pointer_hash<Point>, pointer_equal<Point>>;

    /** Stores a unique set of Edge pointers */
    using EdgeSet = std::unordered_set<Edge*, pointer_hash<Edge>, pointer_equal<Edge>>;

    /** Constructor for MeshUnstructured */
    MeshUnstructured() {}

    /** Destructor for MeshUnstructured */
    virtual ~MeshUnstructured();

    ///////////////////
    // Point methods //
    ///////////////////

    /**
     * Returns either a new Point object or a Point object already in
     * the mesh.
     *
     * @param x The x-coordinate
     * @param y The y-coordinate
     * @return A Point object
     */
    virtual Point& make_point(double x, double y);

    /**
     * Finds an Edge of the mesh with give Point, using the method
     * delaunay::locate.
     *
     * @param p A Point
     * @param e A initial guess for the Edge (optional)
     * @return True if an Edge has the Point, false otherwise
     */
    bool has_point(const Point* p, Edge* e = nullptr);

    /**
     * Inserts a Point into the point database.
     *
     * @return True if the point is new, false otherwise
     */
    bool insert_point(Point* p);

    /**
     * Erases a Point from the point database and frees its memory.
     */
    virtual void delete_point(Point* p);

    /**
     * Gets an unordered set of Point objects from the point database.
     *
     * return An unordered set of Points
     */
    const PointSet& get_points_set() const { return m_points_set; }

    //////////////////
    // Edge methods //
    //////////////////

    /**
     * Constructs and returns a new Edge or returns an existing Edge
     * from the database with both Points provided.
     *
     * @param p0 A Point
     * @param p1 A Point
     * @return An Edge object
     */
    Edge& make_edge(Point& p0, Point& p1);

    /**
     * Inserts a Edge into the edge database. The Edge is not joined
     * or detached from any other edges during insertion.
     *
     * @return True if the Edge is already in the database, false otherwise
     */
    bool insert_edge(Edge* e);

    /**
     * Checks if the Edge or its symmetric counterpart belongs in the
     * edge database.
     *
     * @param e An Edge
     * @return True if the Edge is in the database, false otherwise
     */
    bool has_edge(Edge* e);

    /**
     * An Edge constructor, creating Edges and their Duals by joining an
     * Edge destination with a Point. The Point must on the line of the
     * Edge beyond the destination, or anywhere to the left of the Edge.
     * The new Edge and the Point p are added to the mesh's database.
     *
     * @param e0 The Edge
     * @param p The Point
     * @return An Edge created with e0's destination as its origin
     */
    Edge& extend_edge(Edge& e0, Point& p);

    /**
     * An Edge constructor, creating Edges and their Duals by joining two
     * Edges. The destination of Edge e0 is joined with the origin of Edge
     * e1. The new Edge is added to the mesh's database.
     *
     * @param e0 An Edge
     * @param e1 An Edge
     * @return An Edge created with e0's destination as its origin
     */
    Edge& connect_edges(Edge& e0, Edge& e1);

    /**
     * Delete an Edge, its symmetric Edge, its Dual edges, and
     * neighboring cell Triangle objects, if any. The Edge gets
     * detached from all joined Edges and gets removed from the
     * database.
     *
     * @param e The Edge to delete
     */
    virtual void delete_edge(Edge* e);

    /**
     * Gets an unordered set of Edge objects from the edge database.
     *
     * return An unordered set of Edges
     */
    const EdgeSet& get_edges_set() const { return m_edges_set; }

    /**
     * Gets an Edge from the mesh. It should be assumed to be random.
     *
     * @return An Edge
     */
    Edge* get_initial_edge() { return *m_edges_set.begin(); }

    //////////////////////
    // Triangle methods //
    //////////////////////

    /** Constructs Triangle objects for each cell in the mesh. */
    void fill_triangles();

    /**
     * Constructs Triangle objects for each cell in the mesh. This is
     * intended for filling new cells after splitting a cavity. The
     * input argument designates a place to start filling the cells.
     *
     * @param edge An Edge
     */
    void fill_triangles(Edge& edge);

    /**
     * Frees memory of a Triangle object to the left of an Edge.
     * Pointers to the triangle of all 3 Edges in the cell get
     * cleared.
     *
     * @param e An Edge
     */
    virtual void delete_triangle(Edge* e);

    //////////////////
    // Uniform mesh //
    //////////////////

    /**
     * Adds 4 external points to the database. They are place at
     * distances of 3 times the width and height of the existing
     * points.
     *
     * @return The vector of external Points
     */
    const std::vector<Point*>& add_external_points();

    /**
     * Removes external Points and their Edges from the mesh.  Edges
     * get detached from other Edges and memory of Points and Edges
     * gets freed.
     */
    virtual void remove_external_edges();

    ///////////
    // Other //
    ///////////

    /**
     * Clears the Edge and Point databases in the mesh. No memory is
     * freed, and no edges are freed.
     */
    void clear();

protected:
    /**
     * Constructs a Triangle on the left side of the Edge.  A nullptr
     * gets returns if the cell is not a triangle or if the left side
     * of the Edge is in the exterior of the mesh.
     *
     * @param e An Edge
     * @return A Triangle pointer or a nullptr
     */
    virtual Triangle* make_triangle(Edge& e);

    PointSet m_points_set;
    EdgeSet m_edges_set;
    std::vector<Point*> m_external_points;
};


} // namespace mesh

} // namespace umr

#endif // MESH_UNSTRUCTURED_H
