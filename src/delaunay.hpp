#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <vector>

#include "point.hpp"
#include "mesh.hpp"
#include "quadedge.hpp"


namespace umr {

namespace mesh {

namespace delaunay {


/**
 * Find and Edge closest to or containing a Point
 *
 * @param point A Point
 * @param edge An Edge serving as an initial guess
 * @return The Edge closest to or containing Point p
 */
Edge& locate(const Point& point, Edge& edge);

/**
 * Inserts and Edge using the Giftwrap algorithm, see "Steven
 * Fortune. A Sweepline Algorithm for Voronoi Diagrams. Algorithmica
 * 2(2):153–174, 1987."
 *
 * @param mesh A mesh
 * @param new_edge The Edge to be inserted
 * @param edge An Edge to start from in the mesh
 */
void insert_edge(MeshUnstructured& mesh, Edge& new_edge, Edge& edge);

/**
 * Make a cavity in the mesh per the Delaunay criteria
 *
 * @param mesh A mesh
 * @param point A Point to insert
 * @param edge An Edge nearby the point p
 * @return An Edge on the boundary of the cavity
 */
Edge& make_cavity(MeshUnstructured& mesh, const Point& point, Edge& edge);

/**
 * Splits a cavity using the Divide and Conquer algorithm
 *
 * @param mesh A mesh
 * @param edge An Edge on the boundary of the cavity
 * @param new_points Points within the cavity not currently in the mesh
 */
void split_cavity(MeshUnstructured& mesh, Edge& edge, std::vector<Point*> new_points);

/**
 * Inserts a point into the mesh, swapping edges as needed to meet the
 * Delaunay criteria. No edge deletion occurs (unless the point is
 * exactly on an existing edge).
 *
 * @param mesh A mesh
 * @param point A Point to be inserted
 * @param edge An Edge closest to the Point
 * @return An Edge
 */
Edge& insert_point(MeshUnstructured& mesh, Point& point, Edge& edge);

/**
 * Triangulates a mesh using the Bowyer-Watson algorithm.
 *
 * @param mesh A mesh
 */
void bowyer_watson(MeshUnstructured& mesh);

/**
 * Trianglulates a mesh using the Divide and Conquer algorithm.  The
 * implementation is done per the pseudocode in "Guibas and Stolfi,
 * Primitives for the manipulation of general subdivisions and the
 * computation of Voronoi diagrams, ACM Transactions on Graphics,
 * Vol. 4, No. 2, 1985, 74-123."
 *
 * @param mesh A mesh
 */
void divide_and_conquer(MeshUnstructured& mesh);


} // namespace delaunay

} // namespace mesh

} // namespace umr

#endif // DELAUNAY_H
