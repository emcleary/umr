#include "mesh_generator.hpp"

#include "io.hpp"


namespace umr {


MeshGenerator::MeshGenerator(mesh::MeshUnstructuredConstrained* const mesh,
        std::vector<Real2> point_sources,
        std::list<parametric::IParametric*> segment_sources,
        mesh::algo::IRefinement* const refinement,
        int data_frequency)
        : m_mesh(mesh),
          m_point_sources(std::move(point_sources)),
          m_segment_sources(std::move(segment_sources)),
          m_refinement(refinement),
          m_data_frequency(data_frequency) {
}


MeshGenerator::~MeshGenerator() {
    if (m_mesh)
        delete m_mesh;

    if (m_refinement)
        delete m_refinement;

    for (auto segment : m_segment_sources)
        delete segment;
}


const mesh::MeshUnstructuredConstrained& MeshGenerator::triangulate(int n) {
    initialize();
    run(n);
    return finalize();
}


void MeshGenerator::initialize() {
    if (m_initialized)
        return;

    m_refinement->initialize(*m_mesh, m_point_sources, m_segment_sources);
    m_initialized = true;
}


int MeshGenerator::run(int n) {
    if (m_finished)
        return 0;

    if (m_data_frequency && m_iterations == 0)
        mesh::io::dump(*m_mesh, std::format("iteration_{}.vtk", m_iterations));

    do {
        --n;
        ++m_iterations;
        std::cout << "Iteration " << m_iterations << '\n';
        m_finished = m_refinement->iteration(*m_mesh);
        if (m_finished) {
            std::cout << "Finished refinement at iteration " << m_iterations << '\n';
            break;
        }
        if (m_data_frequency && m_iterations % m_data_frequency == 0)
            mesh::io::dump(*m_mesh, std::format("iteration_{}.vtk", m_iterations));
    } while (n != 0);

    return m_iterations;
}


const mesh::MeshUnstructuredConstrained& MeshGenerator::finalize() {
    if (m_finished && !m_finalized) {
        m_refinement->finalize(*m_mesh);
        m_finalized = true;
    }
    return *m_mesh;
}


bool MeshGenerator::is_converged() {
    return m_finished;
}


int MeshGenerator::get_num_iterations() {
    return m_iterations;
}


void MeshGenerator::dump(std::string&& filename) {
    if (is_vtk_file(filename) || is_json_file(filename))
        mesh::io::dump(*m_mesh, filename, false);
    else
        mesh::io::dump(*m_mesh, filename + ".vtk", false);
}


void MeshGenerator::quality() {
    mesh::refine::mesh_quality(*m_mesh);
}


} // namespace umr
