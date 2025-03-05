#ifndef REFINEMENT_NONE_H
#define REFINEMENT_NONE_H

#include "refinement_interface.hpp"

#include "mesh.hpp"


namespace umr {

namespace mesh {

namespace algo {


/**
 * @class Refinement
 * @brief No refinement
 *
 * This class triangulates a mesh with Delaunay triangulation.
 */
class Refinement : public IRefinement {
public:
    /** Constructor for Refinement */
    Refinement() {}

    /** Destructor for Refinement */
    virtual ~Refinement() {}

    /** Initializes the mesh */
    virtual void initialize(MeshUnstructuredConstrained& mesh,
            std::vector<Real2>& points,
            std::list<parametric::IParametric*>& segments) final;

    /** Run the refinement algorithm */
    virtual bool iteration(MeshUnstructuredConstrained& mesh) final;

    /** Execute any final steps of the refinement */
    virtual void finalize(MeshUnstructuredConstrained& mesh) final {}
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // REFINEMENT_NONE_H
