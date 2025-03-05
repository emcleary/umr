#ifndef REFINEMENT_INTERFACE_H
#define REFINEMENT_INTERFACE_H

#include <list>

#include "mesh.hpp"
#include "parametrics_interface.hpp"


namespace umr {

namespace mesh {

namespace algo {


/**
 * @class IRefinement
 * @brief Interface for the Refinement classes
 *
 * This class serves as an interface for the Refinement classes
 * responsible for refining meshes. Its purpose is to make abstract
 * the steps of the algorithms: initialize, iterate, and finalize.
 */
class IRefinement {
public:
    /** Constructor for IRefinement */
    IRefinement() {}

    /** Destructor for IRefinement */
    virtual ~IRefinement() {}

    /** Initializes the mesh */
    virtual void initialize(MeshUnstructuredConstrained& mesh, std::vector<Real2>& points, std::list<parametric::IParametric*>& segments) = 0;

    /** Run the refinement algorithm iteratively */
    virtual bool iteration(MeshUnstructuredConstrained& mesh) = 0;

    /** Execute any final steps of the refinement */
    virtual void finalize(MeshUnstructuredConstrained& mesh) = 0;
};


} // namespace algo

} // namespace mesh

} // namespace umr

#endif // REFINEMENT_INTERFACE_H
