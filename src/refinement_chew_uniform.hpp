#ifndef REFINEMENT_CHEW_UNIFORM_h
#define REFINEMENT_CHEW_UNIFORM_h

#include "refinement_interface.hpp"

#include "density.hpp"
#include "mesh.hpp"
#include "parametrics.hpp"


namespace umr {

namespace mesh {

namespace algo {


/**
 * @class RefinementChewUniform
 * @brief Chew's Uniform Refinement Algorithm
 *
 * This class generates a mesh using Chew's uniform refinement
 * algorithm.  See "Chew, L. P. 1989 Constrained Delaunay
 * triangulations. Algorithmica 4, 97-108."
 *
 * The original algorithm allows for source segment lengths to be
 * between h and 2*h, but requires careful addition of extra source
 * points to fix all source segments of lengths between sqrt(3)*h and
 * 2*h. This is difficult to automate, and this implementation instead
 * shrinks h to be sufficiently small to work without these additional
 * points. The end result is a mesh finer than necessary.
 */
class RefinementChewUniform : public IRefinement {
public:

    /**
     * @brief Constructs RefinementChewUniform
     *
     * @param density A DensityManager for the algorithm
     */
    RefinementChewUniform(density::DensityManager* density) : m_density(density) {}

    /** Initializes the mesh */
    virtual void initialize(MeshUnstructuredConstrained& mesh,
            std::vector<Real2>& points,
            std::list<parametric::IParametric*>& segments) final;

    /** Run the refinement algorithm */
    virtual bool iteration(MeshUnstructuredConstrained& mesh) final;

    /** Execute any final steps of the refinement */
    virtual void finalize(MeshUnstructuredConstrained& mesh) final;

private:
    double m_hmin;
    Edge* m_edge;
    density::DensityManager* m_density;
    std::stack<Edge*> m_edge_stack;
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // REFINEMENT_CHEW_UNIFORM_h
