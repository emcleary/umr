#include "refinement_chew_uniform.hpp"
#include "refine.hpp"
#include "point.hpp"


namespace umr {

namespace mesh {

namespace algo {


void RefinementChewUniform::initialize(MeshUnstructuredConstrained& mesh,
        std::vector<Real2>& points,
        std::list<parametric::IParametric*>& segments) {

    m_hmin = refine::calculate_minimum_length(points, segments);
    if (m_density)
        m_hmin = std::min(m_hmin, m_density->get_minimum_density());
    m_hmin = refine::split_sources_uniformly(segments, m_hmin);

    // Makes a database of source points
    for (auto [x, y] : points)
        mesh.make_source_point(x, y);

    // Makes a database of source edges (aka segments)
    for (auto segment : segments) {
        std::vector<Edge*> ev = mesh.make_source_edge(segment);
        for (Edge* e : ev)
            m_edge_stack.push(e);
    }

    // Generate the initial mesh
    auto [pmin, pmax] = refine::mesh_point_limits(mesh);
    refine::uniform_grid(mesh, pmin, pmax, m_hmin, 3);

    // Add the source points to the mesh
    m_edge = mesh.get_initial_edge();
    assert(m_edge != nullptr
            && "RefineChewUniform::initialize: no initial edge found!");
    for (Real2 r : points) {
        Point* p = &mesh.make_point(r.first, r.second);
        m_edge = refine::insert_point(mesh, *p, *m_edge, m_hmin);
        m_edge = refine::refine_triangles(mesh, *m_edge, m_hmin);
    }
}


bool RefinementChewUniform::iteration(MeshUnstructuredConstrained& mesh) {

    Edge* ed;
    if (mesh.any_deleted_source_edges()) {
        ed = mesh.get_next_deleted_source_edge();
    } else if (!m_edge_stack.empty()) {
        ed = m_edge_stack.top();
        m_edge_stack.pop();
    } else {
        std::cout <<
            "RefinementChewUniform::iteration: all edges already inserted.\n";
        return true;
    }

    m_edge = refine::insert_point(mesh, ed->org(), *m_edge, m_hmin);
    m_edge = refine::refine_triangles(mesh, *m_edge, m_hmin);

    m_edge = refine::insert_point(mesh, ed->dest(), *m_edge, m_hmin);
    m_edge = refine::refine_triangles(mesh, *m_edge, m_hmin);

    assert(mesh.has_point(&ed->org()) && "Somehow org got deleted!");
    assert(mesh.has_point(&ed->dest()) && "Somehow dest got deleted!");

    // Insert ed if necessary (it might already be present from point
    // insertion)
    if (quadedge::is_free(*ed)) {
        delaunay::insert_edge(mesh, *ed, *m_edge);
        m_edge = refine::refine_triangles(mesh, *ed, m_hmin);
    }

    assert(mesh.has_point(&ed->org()) && "Somehow org got deleted!");
    assert(mesh.has_point(&ed->dest()) && "Somehow dest got deleted!");
    assert(!mesh.get_next_deleted_source_point()
            && "Source point somehow got deleted!");
    assert(mesh.has_edge(ed) && "Edge not inserted properly!");

    m_edge = ed;
    return m_edge_stack.empty() && !mesh.any_deleted_source_edges();
}


void RefinementChewUniform::finalize(MeshUnstructuredConstrained& mesh) {
    assert(!mesh.any_deleted_source_edges()
            && "RefinementChewUniform::finalize: algorith unfinished!");
    mesh.remove_external_edges();
}


} // namespace algo

} // namespace mesh

} // namespace umr
