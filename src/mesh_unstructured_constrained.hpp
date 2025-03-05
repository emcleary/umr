#ifndef MESH_UNSTRUCTURED_CONSTRAINED_H
#define MESH_UNSTRUCTURED_CONSTRAINED_H

#include <stack>
#include <unordered_set>

#include "mesh_unstructured.hpp"
#include "parametrics_interface.hpp"
#include "point.hpp"
#include "quadedge.hpp"
#include "triangle_pq.hpp"
#include "utilities.hpp"


namespace umr {

namespace mesh {


/**
 * @class MeshUnstructuredConstrained
 * @brief A class for managing Point and Edges in a Mesh with sources
 *
 * This class is a tool for building a mesh. It can be seen as a
 * database contains a set of points and edges, a set of contructors
 * and destructors for these objects, and a set of tools for joining
 * edges to build the mesh. In addition, this class contains source
 * Points and source Edges (segments), along with access to a Triangle
 * priority queue to facilitate with refinement algorithms and help
 * manage memory.
 *
 * The constructors and destructors for Points and Edges can be seen
 * as wrappers for those methods of the quadedge namespace. These
 * implementations are specifically designed to prevent any duplicates
 * and memory leaks.
 */
class MeshUnstructuredConstrained : public MeshUnstructured {
public:

    /** Stores a unique pair of segments and their parametrics */
    using SegmentSet = std::unordered_map<Edge*, parametric::IParametric*,
                                          pointer_hash<Edge>, pointer_equal<Edge>>;

    /**
     * Constructor for MeshUnstructuredConstrained. The TrianglePQ
     * type is set to algo::TriCompMinSize.
     */
    MeshUnstructuredConstrained()
            : m_triangle_pq(new algo::TrianglePQ<algo::TriCompMinSize>()) {}

    /**
     * Constructor for MeshUnstructuredConstrained.
     *
     * @param pq A ITrianglePQ pointer
     */
    MeshUnstructuredConstrained(algo::ITrianglePQ* const pq)
            : m_triangle_pq(pq) {}

    /** Destructor for MeshUnstructuredConstrained */
    virtual ~MeshUnstructuredConstrained();

    ///////////////////
    // Point methods //
    ///////////////////

    /**
     * Inserts a point to the mesh database. If the point exists in
     * the source database, than that point gets inserted. Otherwise a
     * new point is constructed.
     *
     * @param x The x-coordinate
     * @param y The y-coordinate
     * @return The Point
     */
    virtual Point& make_point(double x, double y) override;

    /**
     * Constructs and adds a new source Point to the database. If a
     * Point with the given coordinates already exists, it returns the
     * existing Point instead.
     *
     * @param x The x-coordinate
     * @param y The y-coordinate
     * @return The Point
     */
    Point* make_source_point(double x, double y);

    /**
     * Deletes a non-source point from the mesh and frees it from
     * memory. Source points get removed from the mesh and stored for
     * later use.
     *
     * @param p A Point
     */
    virtual void delete_point(Point* p) override;

    /**
     * Checks if a point is a source point
     *
     * @param p A Point
     * @return True if the Point is a source point, false otherwise
     */
    bool is_source_point(Point* p) const;

    /**
     * Check if any source points got deleted
     *
     * @return True if any source points have been deleted, false otherwise
     */
    bool any_deleted_source_points();

    /**
     * Deleted source points are stored in a stack. This method
     * returns top Point from the top of the stack, and pops the
     * element.
     *
     * @return A Point pointer or a nullptr if the stack is empty
     */
    Point* get_next_deleted_source_point();

    /** Gets an unordered set of source points */
    const PointSet& get_source_points() const { return m_source_points; }

    //////////////////
    // Edge methods //
    //////////////////

    /**
     * Constructs and adds a new source segment to the mesh
     * database. It does not get joined to any other edges.
     *
     * @param segment A parametric
     * @return A vector of new subsegments
     */
    std::vector<Edge*> make_source_edge(parametric::IParametric* segment);

    /**
     * Add an Edge and its IParametric segment to source segment
     * database
     *
     * @param e An Edge
     * @param segment A parametric
     */
    void add_source(Edge* e, parametric::IParametric* segment);

    /**
     * Removes the Edge form the mesh. If the Edge is a non-source
     * edge then it gets freed from memory. If it is a source segment
     * then it get stored.
     *
     * @param e An Edge
     */
    virtual void delete_edge(Edge* e) override;

    /**
     * Gets an instance of a source segment from the segment database
     * corresponding to the input Edge. Intended to convert any edge
     * copies (e.g. from a queue) into a corresponding edge object
     * from the database.
     *
     * @param edge An Edge
     * @param An Edge source object
     */
    Edge* get_source(Edge* edge);

    /**
     * Checks if an Edge is a source segment
     *
     * @param e An Edge
     * @return True if the Edge is a source, false otherwise
     */
    bool is_source(Edge* e) const;

    /**
     * Gets the IParametric object coupled with an Edge. Returns
     * nullptr if the Edge is not a source segment.
     *
     * @param edge An Edge
     * @return A IParametric pointer
     */
    parametric::IParametric* get_segment_parametric(Edge* edge);

    /**
     * Splits a source segment into n subsegments of equal length and
     * inserts them into the mesh. The existing segment and its
     * neighboring cells Triangles are freed from memory, creating
     * cavities. The method does not triangulate the cavities.
     *
     * @param edge A source segment to split
     * @param n The number of subsegments to split the segment into (optional)
     * @return A vector of new Edge segment pointers
     */
    std::vector<Edge*> split_edge(Edge& edge, int n = 0);

    /**
     * Check if any source segments got deleted
     *
     * @return True if any source segments have been deleted, false otherwise
     */
    bool any_deleted_source_edges();

    /**
     * Deleted source segments are stored in a stack. This method
     * returns top Edge from the top of the stack, and pops the
     * element.
     *
     * @return An Edge pointer or a nullptr if the stack is empty
     */
    Edge* get_next_deleted_source_edge();

    /** Gets a unordered map of source edges and their parametrics */
    const SegmentSet& get_source_edges() const { return m_source_edges; }

    //////////////////////
    // Triangle methods //
    //////////////////////

    /**
     * Constructs a Triangle on the left side of the Edge.  A nullptr
     * gets returns if the cell is not a triangle or if the left side
     * of the Edge is in the exterior of the mesh.
     *
     * If a triangle does get constructed, the triangle and the edge
     * get pushed to the triangle priority queue.
     *
     * @param e An Edge
     * @return A Triangle pointer or a nullptr
     */
    virtual Triangle* make_triangle(Edge& e);

    /**
     * Deactivates the triangle to the left of the given Edge in the
     * the triangle priority queue. All 3 edges of the cell get their
     * cell data cleared.
     *
     * @param e An Edge
     */
    virtual void delete_triangle(Edge* e);

    /** Gets the triangle priority queue */
    algo::ITrianglePQ& get_triangle_pq() { return *m_triangle_pq; }

    //////////////////
    // Uniform mesh //
    //////////////////

    /**
     * Removes all external Points and their Edges from the mesh up to
     * the source segments.  Edges get detached from other Edges and
     * memory of Points and Edges gets freed.
     */
    virtual void remove_external_edges();

private:
    /** Splits a source segment IParametric, it creates and return a
        set of new source segments corresponding to the parametric's
        subsegments. The parameterics parameters must be set before
        calling this method. */
    std::vector<Edge*> split_source(parametric::IParametric* segment);

    PointSet m_source_points;
    SegmentSet m_source_edges;

    std::stack<Point*> m_deleted_source_points;
    std::stack<Edge*> m_deleted_sources;

    algo::ITrianglePQ* const m_triangle_pq;
};


} // namespace mesh

} // namespace umr

#endif // MESH_UNSTRUCTURED_CONSTRAINED_H
