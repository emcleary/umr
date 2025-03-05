#ifndef REFINE_H
#define REFINE_H

#include <cmath>
#include <memory>

#include "edge_queue.hpp"
#include "delaunay.hpp"
#include "density.hpp"
#include "point.hpp"
#include "quadedge.hpp"
#include "sources.hpp"


namespace umr {

namespace mesh {

namespace refine {


struct ClusterMetrics {
    double length;
    double angle;
    Point apex;
    Edge* edge;
};


using PairFloatParam = std::pair<double, parametric::IParametric*>;


/**
 * Gets the upper and lower bounds of the domain of a mesh based on its source points.
 *
 * @param m A mesh
 * @return A pair of Point objects representing lower bounds and upper bounds
 */
std::pair<Point, Point> mesh_point_limits(const MeshUnstructuredConstrained& mesh);

/**
 * Prints quality measurements of the mesh: angles, lengths, number of triangles.
 *
 * @param m A mesh
 */
void mesh_quality(const MeshUnstructuredConstrained& mesh);

/**
 * Fills the mesh with a uniform grid
 *
 * @param mesh A mesh
 * @param pmin A Point representing the lower bounds of the mesh
 * @param pmax A Point representing the upper bounds of the mesh
 * @param hmin Edge length of the triangles
 * @param nextra Number of extra triangles to add beyond the bounds
 */
void uniform_grid(MeshUnstructured& mesh, const Point& pmin, const Point& pmax,
        double hmin, int nextra = 0);

/**
 * A function to make and split cavities using the Delaunay namespace
 * methods (Divide and Conquer). Does not check for source points or edges.
 *
 * @param mesh A mesh
 * @param edge An edge
 * @param hmin A lengthscale for the radius of a cavity to be made and split
 * @return An edge guaranteed to exist in the mesh
 */
Edge* refine_triangles(MeshUnstructuredConstrained& mesh, Edge& edge, double hmin);

/**
 * A function that inserts a point into the mesh. Neighboring
 * non-source points and edges within a distance of hmin will get deleted.
 * Remaining non-source edges will get swapped as needed to preserve Delaunay
 * triangularization. Triangularization uses the refine::make_and_split_cavity method.
 *
 * @param mesh A mesh
 * @param point The point to insert
 * @param edge An edge close to the point being inserted
 * @param hmin The length scale used to delete nearby points
 * @return An edge guaranteed to be in the mesh
 */
Edge* insert_point(MeshUnstructuredConstrained& mesh, const Point& point,
        Edge& edge, double hmin);

/**
 * A function to insert an edge into the mesh. This function uses the
 * method refine::refine_triangles.
 *
 * @param mesh A mesh
 * @param edge The edge to insert
 * @param hmin The length used to delete nearby points/edges
 * @return An edge guaranteed to be in the mesh
 */
Edge* insert_edge(MeshUnstructuredConstrained& mesh, Edge& edge, double hmin);

/**
 * Calculates the shortest lengthscale of a given input mesh using its
 * source points and segments.
 *
 * @param source_points A vector of 2d points
 * @param source_segments A list of parametric::IParametric objects
 * @return A lengthscale
 */
double calculate_minimum_length(const std::vector<Real2>& source_points,
        const std::list<parametric::IParametric*>& source_segments);

/**
 * Splits a set of source segments into subsegments of roughly uniform lengths.
 * Subsegments should fall in a range scaled by sqrt(3). Subsegments are inserted
 * into the supplied list.
 *
 * @param source_segments A list of parametric::IParametric objects
 * @return The smallest subsegment length
 */
double split_sources_uniformly(std::list<parametric::IParametric*>& source_segments,
        double target_hmin);

/**
 * Finds the edge closest to the provided target point or the first
 * source segment found on the way. If the target point is already in
 * the mesh, the edge returned will have that points as it org.
 *
 * @param mesh A mesh
 * @param point The Point being targeted
 * @param edge The starting Edge
 * @return An closest to or containing the target point
 */
Edge& locate_or_closest_edge(const MeshUnstructuredConstrained& mesh,
        const Point& point, Edge& edge);

/**
 * A convergence test comparing the given triangle to a given angle.
 *
 * @param mesh A mesh
 * @param triangle A Triangle to be tested
 * @param min_angle The convergence criterion
 * @result Returns true if the test passes
 */
bool convergence_test_angle(const Triangle& triangle, double min_angle);

/**
 * A convergence test comparing the given triangle to the density. It
 * is converged if the triangle's radius is less than or equal to the
 * density evaluated at the triangle's circumcenter.
 *
 * @param triangle A Triangle to be tested
 * @param density The convergence criterion
 * @result Returns true if the test passes
 */
bool convergence_test_density(const Triangle& triangle,
        const density::DensityManager& density);

/**
 * Fills a queue with segments encroached by the given point.
 *
 * @param mesh A mesh
 * @param edge A starting edge
 * @param pc The Point used for testing edges
 * @return An EdgeQueue filled with encroached edge segments
 */
algo::EdgeQueue fill_encroached_edge_queue(const MeshUnstructuredConstrained& mesh,
        Edge& edge, const Point& pc, double const triangle_radius);

/**
 * Check if a segment belongs to a cluster of input segments within a
 * given angle of each other.
 *
 * @param mesh A mesh
 * @param encroached An encroached Edge segment
 * @param min_angle_allowed Minimum angle to define a cluster, default to 60 degrees
 * @return Cluster data
 */
ClusterMetrics get_segment_cluster_metrics(MeshUnstructuredConstrained& mesh,
        Edge& encroached, double const min_angle_allowed = 60);

/**
 * Split a segment evenly into a given number of subsegments. Neighboring cells
 * are triangularized using refine::make_and_split_cavity.
 *
 * @param m A mesh
 * @param edge The segment to be split
 * @param num_subsegments The number of subsegments to split the segment into
 * @return A new subsegment, guaranteed to be in the mesh
 */
Edge* split_edge(MeshUnstructuredConstrained& mesh, Edge& edge, int num_segments = 2);

/**
 * Split an encroached segment per the Terminator algorithm. See "Shewchuck,
 * J. R. 2000. Mesh generation for domains with small
 * angles. Computational Geometry."
 *
 * Splitting only occurs if the segment length is less than
 * min_length_allowed in an effort to allow some density metric
 * control when the algorithm starts.
 *
 * @param mesh A mesh
 * @param encroached An encroached segment
 * @param min_length_allowed Lengthscale below which the Terminator algorithm starts
 * @param cluster Cluster segment data
 * @param num_segments Number of subsegments to split the segment into
 * if Terminator splitting does not occur.
 * @return Returns a vector newly split subsegments
 */
std::vector<Edge*> split_encroached_edge_terminator(
        MeshUnstructuredConstrained& mesh, Edge& encroached,
        double const min_length_allowed, const ClusterMetrics& cluster,
        int num_segments = 2);

/**
 * Split a segment into 2 subsegments at a given distance from a
 * segment cluster's apex.
 *
 * @param mesh A mesh
 * @param encroached An encroached segment
 * @param length The distance from the apex to split the segment
 * @param cluster Data of the segment cluster
 * @return Returns a vector of newly split subsegments
 */
std::vector<Edge*> segment_splitter_nonuniform(MeshUnstructuredConstrained& mesh,
        Edge& encroached, double length, const ClusterMetrics& cluster);

/**
 * Makes a cavity in the mesh near a given point (to be inserted). No
 * points get deleted. This is intended for use with Ruppert's
 * refinement algorithm when making cavities near encroached segment
 * and no existing encroaching point exist in the mesh.
 *
 * @param mesh A mesh
 * @param edge A starting Edge, typically an encroached subsegment
 * @param point A Point (to be inserted later)
 */
void make_cavity(MeshUnstructuredConstrained& mesh, Edge& edge, const Point& point);

/**
 * Splits a cavity starting from a given edge. It's implementation is
 * adapted from the "Rising Bubble" algorithm of "Guibas and Stolfi,
 * Primitives for the manipulation of general subdivisions and the
 * computation of Voronoi diagrams, ACM Transactions on Graphics,
 * Vol. 4, No. 2, 1985, 74-123." It's implementation preserved any input
 * source segments from getting swapped.
 *
 * @param mesh A mesh
 * @param basel The starting Edge
 */
void split_cavity(MeshUnstructuredConstrained& mesh, Edge& basel);

/**
 * Makes and splits a cavity starting from a given edge. It's
 * implementation is adapted from the "Rising Bubble" algorithm of
 * "Guibas and Stolfi, Primitives for the manipulation of general
 * subdivisions and the computation of Voronoi diagrams, ACM
 * Transactions on Graphics, Vol. 4, No. 2, 1985, 74-123." It's
 * implementation is designed to delete points too close to a point
 * recently inserted while splitting the cavity. It is designed to
 * preserve all source points and source segments.
 *
 * @param mesh A mesh
 * @param basel The starting Edge
 * @param min_length A distance from the Point within which other
 * points can get deleted.
 * @param point A Point at the center of the cavity made and split
 */
std::pair<Edge*, bool> make_and_split_cavity(MeshUnstructuredConstrained& mesh,
        Edge& basel, const double min_length, const Point& point);


} // namespace refine

} // namespace mesh

} // namespace umr

#endif // REFINE_H
