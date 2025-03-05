#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "parametrics.hpp"
#include "refinement.hpp"


namespace umr {


/**
 * @class MeshGenerator
 * @brief Builds, refines, and dumps a mesh
 *
 * This object gets built by the Builder and is responsible
 * for generating the mesh, dumping mesh data for plotting,
 * and printing mesh quality.
 */
class MeshGenerator {
public:

    /**
     * @brief MeshGenerator constructor
     *
     * @param mesh Mesh pointer
     * @param point_sources Input source points
     * @param segment_sources Input source segments
     * @param refinement Refinement object
     * @param data_frequency Dump mesh data every N iterations
     */
    MeshGenerator(mesh::MeshUnstructuredConstrained* const mesh,
            std::vector<Real2> point_sources,
            std::list<parametric::IParametric*> segment_sources,
            mesh::algo::IRefinement* const refinement,
            int data_frequency);

    /** Descructor for MeshGenerator */
    ~MeshGenerator();

    /**
     * @brief Triangulates the mesh
     *
     * Triangulates the mesh, either running for n iterations or until
     * refinement converges, whichever comes first.
     *
     * @param n Max number of iterations
     * @return A pointer to the generated mesh
     */
    const mesh::MeshUnstructuredConstrained& triangulate(int n = 0);

    /** Initializes the mesh */
    void initialize();

    /**
     * @brief Triangulates the mesh
     *
     * Triangulates the mesh, either running for n iterations or until
     * refinement converges, whichever comes first.
     *
     * @param n Max number of iterations
     */
    int run(int n = 0);

    /** Finalizes the mesh */
    const mesh::MeshUnstructuredConstrained& finalize();

    /** Checks if refinement is converged */
    bool is_converged();

    /** Returns the number of iterations */
    int get_num_iterations();

    /**
     * Dumps mesh data to a file. The type is determined by its
     * extension" "vtk" or "json".
     */
    void dump(std::string&& filename);

    /** Prints measurements of mesh quality */
    void quality();

private:
    mesh::MeshUnstructuredConstrained* const m_mesh;
    std::vector<Real2> m_point_sources;
    std::list<parametric::IParametric*> m_segment_sources;
    mesh::algo::IRefinement* const m_refinement;
    int m_data_frequency;
    bool m_initialized = false;
    bool m_finished = false;
    bool m_finalized = false;
    int m_iterations = 0;
};


} // namespace umr

#endif // MESH_GENERATOR_H
