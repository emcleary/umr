#include "mesh_unstructured_constrained.hpp"

#include <iostream>
#include <queue>
#include <unordered_set>

#include "delaunay.hpp"
#include "math.hpp"
#include "mesh_unstructured.hpp"
#include "triangle_pq.hpp"


namespace umr {

namespace mesh {


MeshUnstructuredConstrained::~MeshUnstructuredConstrained() {
    for (auto [edge, parametric] : m_source_edges)
        delete parametric;

    while (!m_deleted_sources.empty()) {
        Edge* e = m_deleted_sources.top();
        m_deleted_sources.pop();
        m_edges_set.erase(e);
        m_edges_set.erase(&e->sym());
        mesh::quadedge::delete_edge_all(e);
    }

    while (!m_deleted_source_points.empty()) {
        Point* p = m_deleted_source_points.top();
        m_deleted_source_points.pop();
        m_points_set.erase(p);
        delete p;
    }

    delete m_triangle_pq;
}


///////////////////
// Point methods //
///////////////////

Point& MeshUnstructuredConstrained::make_point(double x, double y) {
    Point* p = new Point(x, y);
    auto it = m_source_points.find(p);
    if (it != m_source_points.end()) {
        delete p;
        m_points_set.insert(*it);
        return **it;
    }
    auto ret = m_points_set.insert(p);
    if (!ret.second)
        delete p;
    return **ret.first;
}


Point* MeshUnstructuredConstrained::make_source_point(double x, double y) {
    Point* p = new Point(x, y);
    auto ret = m_source_points.insert(p);
    if (!ret.second)
        delete p;
    else {
        auto it = m_points_set.find(p);
        if (it != m_points_set.end())
            assert(*ret.first == *it
                    && "Point and source point are different objects!");
    }
    return *ret.first;
}


void MeshUnstructuredConstrained::delete_point(Point* p) {
    assert(p);
    if (m_source_points.contains(p)) {
        m_deleted_source_points.push(p);
        m_points_set.erase(p);
    } else {
        MeshUnstructured::delete_point(p);
    }

}


bool MeshUnstructuredConstrained::is_source_point(Point* p) const {
    assert(p);
    auto it = m_source_points.find(p);
    if (it == m_source_points.end())
        return false;
    return p == *it;
}


bool MeshUnstructuredConstrained::any_deleted_source_points() {
    return !m_deleted_source_points.empty();
}


Point* MeshUnstructuredConstrained::get_next_deleted_source_point() {
    if (m_deleted_source_points.empty())
        return nullptr;
    Point* p = m_deleted_source_points.top();
    m_deleted_source_points.pop();
    return p;
}


//////////////////
// Edge methods //
//////////////////


std::vector<Edge*> MeshUnstructuredConstrained::make_source_edge(
        parametric::IParametric* segment) {

    assert(segment);
    std::vector<double>& params = segment->get_parameters();
    assert(params.size() >= 2);
    Real2 r0 = segment->evaluate(params[0]);
    Point* p0 = make_source_point(r0.first, r0.second);
    make_point(r0.first, r0.second);
    std::vector<Edge*> edges;
    for (size_t i = 1; i < params.size(); ++i) {
        auto [x, y] = segment->evaluate(params[i]);
        Point* p1 = make_source_point(x, y);
        make_point(x, y);
        Edge* e = &make_edge(*p0, *p1);
        edges.push_back(e);
        parametric::IParametric* subseg = segment->clone();
        subseg->set_parameter_bounds(params[i-1], params[i]);
        add_source(e, subseg);
        p0 = p1;
    }
    return edges;
}


void MeshUnstructuredConstrained::add_source(Edge* e, parametric::IParametric* segment) {
    assert(e);
    if (m_source_edges.contains(e)) {
        std::cout << "Sources already contains edge\n";
        exit(-1);
    }
    m_source_edges[e] = segment;
}


void MeshUnstructuredConstrained::delete_edge(Edge* e) {
    assert(e);
    if (e->left().get_data())
        delete_triangle(e);
    if (e->right().get_data())
        delete_triangle(&e->sym());

    m_edges_set.erase(e);
    m_edges_set.erase(&e->sym());
    mesh::quadedge::remove_edge(e);
    if (is_source(e)) {
        m_deleted_sources.push(e);
    } else {
        mesh::quadedge::delete_edge(e);
    }
}


Edge* MeshUnstructuredConstrained::get_source(Edge* edge) {
    assert(edge);
    auto it = m_source_edges.find(edge);
    if (it != m_source_edges.end())
        return it->first;
    auto itsym = m_source_edges.find(&edge->sym());
    if (itsym != m_source_edges.end())
        return itsym->first;
    return nullptr;
}


bool MeshUnstructuredConstrained::is_source(Edge* e) const {
    assert(e);
    Edge* source = nullptr;
    auto itsym = m_source_edges.find(&e->sym());
    if (itsym == m_source_edges.end()) {
        auto it = m_source_edges.find(e);
        if (it != m_source_edges.end())
            source = it->first;
    } else {
        source = &(itsym->first->sym());
    }
    return e == source;
}


parametric::IParametric* MeshUnstructuredConstrained::get_segment_parametric(Edge* edge) {
    assert(edge);
    auto it = m_source_edges.find(edge);
    if (it != m_source_edges.end())
        return it->second;
    return nullptr;
}


std::vector<Edge*> MeshUnstructuredConstrained::split_edge(Edge& edge, int n) {
    assert(n != 1 && "Must split a segment into at least 2 subsegments!");

    // ensure the edge has no cells
    if (edge.left().get_data())
        delete_triangle(&edge);

    if (edge.right().get_data())
        delete_triangle(&edge.sym());

    // ensure the split edge is a source
    auto it = m_source_edges.find(&edge);
    bool flip = false;
    if (it == m_source_edges.end()) {
        it = m_source_edges.find(&edge.sym());
        flip = true;
    }
    assert(it != m_source_edges.end());

    auto segment = (*it).second;
    if (n != 0) {
        segment->set_num_subsegments(n);
        segment->optimize_parameters();
    }
    std::vector<Edge*> new_edges = split_source(segment);

    if (flip) {
        for (size_t i = 0; i < new_edges.size(); ++i)
            new_edges[i] = &new_edges[i]->sym();
        std::reverse(new_edges.begin(), new_edges.end());
    }

    // link edges together
    for (size_t i = 0; i < new_edges.size()-1; ++i) {
        assert(&new_edges[i+1]->onext() == new_edges[i+1]);
        mesh::quadedge::splice(new_edges[i]->sym(), *new_edges[i+1]);
    }

    // add new edges to the mesh
    assert(edge.org() == new_edges.front()->org());
    assert(edge.dest() == new_edges.back()->dest());
    mesh::quadedge::splice(*new_edges.front(), edge);
    mesh::quadedge::splice(new_edges.back()->sym(), edge.sym().oprev());

    // remove the edge from sources -- ensures it doesn't get stored upon deletion
    delete it->second;
    m_source_edges.erase(it);

    // remove the old edge
    delete_edge(&edge);

    return new_edges;
}


bool MeshUnstructuredConstrained::any_deleted_source_edges() {
    return !m_deleted_sources.empty();
}


Edge* MeshUnstructuredConstrained::get_next_deleted_source_edge() {
    if (m_deleted_sources.empty())
        return nullptr;
    Edge* e = m_deleted_sources.top();
    m_deleted_sources.pop();
    return e;
}


std::vector<Edge*> MeshUnstructuredConstrained::split_source(
        parametric::IParametric* segment) {
    assert(segment);
    std::vector<Edge*> new_edges;
    std::vector<double>& params = segment->get_parameters();
    assert(params.size() >= 2);

    Real2 r0 = segment->evaluate(params[0]);
    Point* p0 = make_source_point(r0.first, r0.second);
    insert_point(p0);
    for (size_t i = 1; i < params.size(); ++i) {
        auto [x, y] = segment->evaluate(params[i]);
        Point* p1 = make_source_point(x, y);
        insert_point(p1);

        Edge* e = &make_edge(*p0, *p1);
        new_edges.push_back(e);
        parametric::IParametric* subseg =
            segment->split(params[i-1], params[i]);
        add_source(e, subseg);
        p0 = p1;
    }

    return new_edges;
}


//////////////////////
// Triangle methods //
//////////////////////


mesh::Triangle* MeshUnstructuredConstrained::make_triangle(Edge& e) {
    mesh::Triangle* t = mesh::quadedge::make_triangle(e);
    if (t)
        m_triangle_pq->push(t, &e);
    return t;
}


void MeshUnstructuredConstrained::delete_triangle(Edge* e) {
    assert(mesh::quadedge::part_of_triangle(*e));
    m_triangle_pq->deactivate(e->left().get_data());
    e->left().clear_data();
    e->lnext().left().clear_data();
    e->lnext().lnext().left().clear_data();
}


//////////////////
// Uniform mesh //
//////////////////

void MeshUnstructuredConstrained::remove_external_edges() {
    // Find edge on the boundary of the mesh
    Edge* edge = get_initial_edge();
    for (Edge* e : m_edges_set) {
        if (e->org().x < edge->org().x)
            edge = e;
    }

    while (true) {
        bool is_exterior = !mesh::quadedge::part_of_triangle(edge->sym());
        is_exterior |= !math::is_ccw(edge->sym().org(), edge->sym().lnext().org(),
                edge->sym().lnext().lnext().org());
        if (is_exterior)
            break;
        edge = &edge->onext();
    }

    // Confirm 1 side is triangle and the other is not (interior) triangle
    bool is_interior = mesh::quadedge::part_of_triangle(*edge);
    is_interior &= math::is_ccw(edge->org(), edge->lnext().org(),
            edge->lnext().lnext().org());
    bool is_exterior = !mesh::quadedge::part_of_triangle(edge->sym());
    is_exterior |= !math::is_ccw(edge->sym().org(), edge->sym().lnext().org(),
            edge->sym().lnext().lnext().org());
    assert((is_interior == is_exterior)
            && "MeshUnstructuredConstrained::remove_external_edges: "
            "Starting edge not on the boundary!");
    // need edge to be interior
    if (is_interior) // sym is interior
        edge = &edge->sym();

    // must queue 1 exterior edge per region
    std::unordered_set<Edge*> in_queue;
    std::queue<Edge*> queue;
    bool on_source = true; // ensure at least 1 edge gets queued if all edges are non-sources
    Edge* efin = edge;
    do {
        bool curr_source = is_source(edge);
        if (!curr_source && on_source) {
            queue.push(&edge->sym());
            in_queue.insert(&edge->sym());
        }
        on_source = curr_source;
        edge = &edge->lnext();
    } while (edge != efin);

    if (queue.empty()) {
        std::cout << "MeshUnstructuredConstrained::remove_external_edges: "
            "No exterior edges found!\n";
        return;
    }

    // remove edges
    auto add_to_queue = [&](Edge* e) {
        if (!is_source(e)) {
            if (!(in_queue.contains(e) || in_queue.contains(&e->sym()))) {
                queue.push(e);
                in_queue.insert(e);
            }
        }
    };

    while (!queue.empty()) {
        Edge* e = queue.front();
        queue.pop();
        //in_queue.erase(e);
        add_to_queue(&e->lnext().sym());
        add_to_queue(&e->lnext().lnext().sym());
        Point* po = e == &e->onext() ? &e->org() : nullptr;
        Point* pd = &e->sym() == &e->sym().onext() ? &e->dest() : nullptr;
        delete_edge(e);
        if (po) delete_point(po);
        if (pd) delete_point(pd);
    }

    m_external_points.clear();
}


} // namespace mesh

} // namespace umr
