#ifndef REFINEMENT_CHEW_NONUNIFORM_H
#define REFINEMENT_CHEW_NONUNIFORM_H

#include "refinement_interface.hpp"

#include <cmath>
#include <memory>
#include <list>

#include "density.hpp"
#include "edge_queue.hpp"
#include "mesh.hpp"
#include "parametrics_interface.hpp"
#include "refine.hpp"
#include "refinement_none.hpp"


namespace umr {

namespace mesh {

namespace algo {


/**
 * @class RefinementChewNonuniform
 * @brief Chew's Nonuniform Refinement Algorithm
 *
 * This class generates a mesh using Chew's non-uniform refinement
 * algorithm.  See "Chew, L. P. 1993 Guaranteed-quality mesh
 * generation for curved surfaces. Proceedings of the Ninth Annual
 * Symposium for Computational Computing Geometry, 274-280."
 *
 * This algorithm has guaranteed convergence for minimum angles of up
 * to 26.5 degrees, but often converges for angles up to 30 degrees.
 *
 * Shewchuck's Terminator algorithm is implemented to ensure
 * convergence for geometries with small input angles. See "Shewchuck,
 * J. R. 2000. Mesh generation for domains with small
 * angles. Computational Geometry."
 */
class RefinementChewNonuniform : public IRefinement {
public:
    /**
     * @brief Constructs RefinementChewUniform
     *
     * @param min_angle The targeted minimum angle for refinement
     * @param density A DensityManager for the algorithm
     * @param use_terminator Use the Terminator algorithm. On by default.
     */
    RefinementChewNonuniform(double min_angle, density::DensityManager* density,
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

    /** Splits a segment encroached by a proposed point. */
    bool split_encroached_segment(MeshUnstructuredConstrained& mesh,
            Edge& edge, double min_length, int num_segments,
            bool terminator_only = false);

    double m_hmin, m_min_angle;
    density::DensityManager* m_density;
    bool m_use_terminator;
    algo::EdgeQueue m_encroached_queue;
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // REFINEMENT_CHEW_NONUNIFORM_H
