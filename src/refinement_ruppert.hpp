#ifndef REFINEMENT_RUPPERT_H
#define REFINEMENT_RUPPERT_H

#include "refinement_interface.hpp"

#include <list>
#include <vector>

#include "density.hpp"
#include "edge_queue.hpp"
#include "mesh.hpp"
#include "parametrics_interface.hpp"
#include "refine.hpp"


namespace umr {

namespace mesh {

namespace algo {


/**
 * @class RefinementRuppert
 * @brief Rupperts Refinement Algorithm
 *
 * This class generates a mesh using Ruppert's refinement
 * algorithm.  See "Ruppert, J. 1995. A Delaunay refinement algorithm
 * for quality 2-dimensional mesh generation. Journal of Algorithms."
 *
 * This algorithm has guaranteed convergence for minimum angles of up
 * to 20.7 degrees, but often converges for angles up to 30 degrees.
 *
 * Shewchuck's Terminator algorithm is implemented to ensure
 * convergence for geometries with small input angles. See "Shewchuck,
 * J. R. 2000. Mesh generation for domains with small
 * angles. Computational Geometry."
 */
class RefinementRuppert : public IRefinement {
public:
    /**
     * @brief Constructs RefinementChewUniform
     *
     * @param min_angle The targeted minimum angle for refinement
     * @param density A DensityManager for the algorithm
     * @param use_terminator Use the Terminator algorithm. On by default.
     */
    RefinementRuppert(double min_angle, density::DensityManager* density,
            bool use_terminator = true);

    /** Initializes the mesh */
    virtual void initialize(MeshUnstructuredConstrained& mesh,
            std::vector<Real2>& points,
            std::list<parametric::IParametric*>& segments) final;

    /** Run the refinement algorithm */
    virtual bool iteration(MeshUnstructuredConstrained& mesh) final;

    /** Execute any final steps of the refinement */
    virtual void finalize(MeshUnstructuredConstrained& mesh) final {}

private:

    /** Check if an edge is encroached by the opposing point in its triangle. */
    bool is_encroached(MeshUnstructuredConstrained& mesh, Edge& edge);

    /** Check if an edge is encroached by a point. */
    bool is_encroached(MeshUnstructuredConstrained& mesh,
            Edge& edge, Point& point);

    /** Splits a segment encroached by a proposed point. */
    bool split_encroached_segment(MeshUnstructuredConstrained& mesh,
            Edge& edge, double min_length);

    /** Fills encroached queue after splitting and before constructing triangles */
    void queue_encroached(MeshUnstructuredConstrained& mesh, Edge& edge);

    double m_hmin, m_min_angle;
    density::DensityManager* m_density;
    bool m_use_terminator;
    algo::EdgeQueue m_encroached_queue;
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // REFINEMENT_RUPPERT_H
