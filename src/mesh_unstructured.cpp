#include "mesh_unstructured.hpp"

#include <iostream>
#include <queue>
#include <unordered_set>

#include "delaunay.hpp"
#include "math.hpp"


namespace umr {

namespace mesh {


MeshUnstructured::~MeshUnstructured() {
    for (auto e : m_edges_set)
        mesh::quadedge::delete_edge_all(e);
    for (auto p : m_points_set) {
        delete p;
    }
}


///////////////////
// Point methods //
///////////////////


Point& MeshUnstructured::make_point(double x, double y) {
    Point* p = new Point(x, y);
    auto ret = m_points_set.insert(p);
    if (!ret.second)
        delete p;
    return **ret.first;
}


bool MeshUnstructured::has_point(const Point* p, Edge* e) {
    assert(p);
    if (!e)
        e = get_initial_edge();
    Edge& ef = delaunay::locate(*p, *e);
    return ef.org() == *p;
}


bool MeshUnstructured::insert_point(Point* p) {
    assert(p);
    auto ret = m_points_set.insert(p);
    return ret.second;
}


void MeshUnstructured::delete_point(Point* p) {
    assert(p);
    m_points_set.erase(p);
    delete p;
}


//////////////////
// Edge methods //
//////////////////


Edge& MeshUnstructured::make_edge(Point& p0, Point& p1) {
    Edge* e = &mesh::quadedge::make_edge(p0, p1);
    auto iter = m_edges_set.find(&e->sym());
    if (iter == m_edges_set.end()) {
        auto ret = m_edges_set.insert(e);
        assert(has_edge(e) && "make_edge: Missing edge after inserted");
        if (!ret.second) // element already exists
            mesh::quadedge::delete_edge_all(e);
        return **ret.first;
    }
    assert(m_edges_set.find(e) == m_edges_set.end()
            && "Only one direction of each edge should be tracked!");
    mesh::quadedge::delete_edge_all(e);
    return (**iter).sym();
}


bool MeshUnstructured::insert_edge(Edge* e) {
    assert(e);
    bool has_e = m_edges_set.contains(e);
    bool has_esym = m_edges_set.contains(&e->sym());
    if (has_e || has_esym)
        return false;
    m_edges_set.insert(e);
    return true;
}


bool MeshUnstructured::has_edge(Edge* e) {
    assert(e);
    if (m_edges_set.contains(e))
        return true;
    return m_edges_set.contains(&e->sym());
}


Edge& MeshUnstructured::extend_edge(Edge& e0, Point& p) {
    assert(!math::right_of(p, e0) && "extend_edge must have point "
            "to the left of or on the line containing the edge");
    Edge& e1 = this->make_edge(e0.dest(), p);

    if (math::on_line(p, e0))
        mesh::quadedge::splice(e0.sym(), e1);
    else
        mesh::quadedge::splice(e0.sym().oprev(), e1);
    return e1;
}


Edge& MeshUnstructured::connect_edges(Edge& e0, Edge& e1) {
    Edge& e2 = make_edge(e0.dest(), e1.org());
    mesh::quadedge::splice(e2, e0.lnext());
    mesh::quadedge::splice(e2.sym(), e1);
    return e2;
}


void MeshUnstructured::delete_edge(Edge* e) {
    assert(e);
    m_edges_set.erase(e);
    m_edges_set.erase(&e->sym());
    mesh::quadedge::delete_edge_all(e);
}


//////////////////////
// Triangle methods //
//////////////////////


void MeshUnstructured::fill_triangles() {
    Edge* e = get_initial_edge();
    fill_triangles(*e);
}


void MeshUnstructured::fill_triangles(Edge& edge) {
    std::queue<Edge*> queue;
    queue.push(&edge);
    queue.push(&edge.sym());
    while (!queue.empty()) {
        Edge* e = queue.front();
        queue.pop();

        mesh::Triangle* t0 = e->left().get_data();
        if (t0) continue;

        mesh::Triangle* t = make_triangle(*e);
        if (!t) continue; // nullptr if not ccw triangle

        e->left().set_data(t); // automatically set to neighboring duals

        queue.push(&e->lnext().sym());
        queue.push(&e->lnext().lnext().sym());
    }
}


mesh::Triangle* MeshUnstructured::make_triangle(Edge& e) {
    return mesh::quadedge::make_triangle(e);
}


void MeshUnstructured::delete_triangle(Edge* e) {
    assert(e);
    assert(mesh::quadedge::part_of_triangle(*e));
    delete e->left().get_data();
    e->left().clear_data();
    e->lnext().left().clear_data();
    e->lnext().lnext().left().clear_data();
}


//////////////////
// Uniform mesh //
//////////////////


const std::vector<Point*>& MeshUnstructured::add_external_points() {
    if (!m_external_points.empty()) {
        std::cout << "External points already exist.\n";
        return m_external_points;
    }

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    for (auto p : m_points_set) {
        xmin = std::min(xmin, p->x);
        xmax = std::max(xmax, p->x);
        ymin = std::min(ymin, p->y);
        ymax = std::max(ymax, p->y);
    }

    double dx = 3.0 * (xmax - xmin);
    double dy = 3.0 * (ymax - ymin);
    xmin -= dx;
    xmax += dx;
    ymin -= dy;
    ymax += dy;

    Point* p0 = &make_point(xmin, ymin);
    Point* p1 = &make_point(xmax, ymin);
    Point* p2 = &make_point(xmax, ymax);
    Point* p3 = &make_point(xmin, ymax);
    m_external_points.push_back(p0);
    m_external_points.push_back(p1);
    m_external_points.push_back(p2);
    m_external_points.push_back(p3);
    return m_external_points;
}


void MeshUnstructured::remove_external_edges() {
    if (m_external_points.size() == 0) {
        std::cout << "No external points exist. No edges can be removed.\n";
        return;
    }

    assert(m_external_points.size() == 4);
    Point* p0 = m_external_points[0];
    Point* p1 = m_external_points[1];
    Point* p2 = m_external_points[2];
    Point* p3 = m_external_points[3];

    Edge* e = get_initial_edge();
    e = &delaunay::locate(*p0, *e);
    assert(e->org() == *p0);
    while (e->dest() != *p1 && e->dest() != *p2)
        e = &e->onext();

    bool is_ccw = math::is_ccw(e->org(), e->dest(), *p3);
    if (is_ccw)
        e = &e->oprev();

    while (e != &e->oprev())
        delete_edge(&e->oprev());
    e = &e->lnext();
    while (e != &e->oprev())
        delete_edge(&e->oprev());
    e = &e->lnext();
    while (e != &e->oprev())
        delete_edge(&e->oprev());
    e = &e->lnext();
    while (e != &e->oprev())
        delete_edge(&e->oprev());
    delete_edge(e);

    for (Point* p : m_external_points)
        delete_point(p);
    m_external_points.clear();
}


///////////
// Other //
///////////

void MeshUnstructured::clear() {
    m_edges_set.clear();
    m_points_set.clear();
}


} // namespace mesh

} // namespace umr
