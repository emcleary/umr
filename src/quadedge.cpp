#include "quadedge.hpp"

#include <cassert>

#include "math.hpp"
#include "point.hpp"


namespace umr {

namespace mesh {

namespace quadedge {


Edge& make_edge(Point& p0, Point& p1) {
    Edge* e0 = new Edge(p0, p1);
    Edge* e1 = new Edge(p1, p0);
    Dual* d0 = new Dual(p0, p1);
    Dual* d1 = new Dual(p1, p0);
    e0->m_right = d0;
    d0->m_right = e1;
    e1->m_right = d1;
    d1->m_right = e0;
    e0->m_onext = e0;
    e1->m_onext = e1;
    d0->m_onext = d1;
    d1->m_onext = d0;
    return *e0;
};


Edge& extend_edge(Edge& e0, Point& p) {
    assert(!math::right_of(p, e0) && "extend_edge must have point "
            "to the left of or on the line containing the edge");
    Edge& e1 = make_edge(e0.dest(), p);

    if (math::on_line(p, e0))
        splice(e0.sym(), e1);
    else
        splice(e0.sym().oprev(), e1);
    return e1;
};


Edge& connect(Edge& e0, Edge& e1) {
    Edge& e2 = make_edge(e0.dest(), e1.org());
    splice(e2, e0.lnext());
    splice(e2.sym(), e1);
    return e2;
};


Triangle* make_triangle(Edge& e) {
    if (!part_of_triangle(e))
        return nullptr;
    Point& p0 = e.org();
    Point& p1 = e.dest();
    Point& p2 = e.lnext().dest();
    if (!math::is_ccw(p0, p1, p2))
        return nullptr;
    return new Triangle(p0, p1, p2);
}


void delete_cell(Edge* e) {
    // Free memory -- cell data
    auto delete_cell = [](Dual* d) {
        if (d->m_data) {
            size_t count = 0;
            Triangle* t = d->m_data;
            while (d->m_data) {
                ++count;
                d->m_data = nullptr;
                d = &d->onext();
            }
            assert(count == 3 && "Only setup for triangle!");
            delete t;
        }
    };
    delete_cell(&e->left());
    delete_cell(&e->right());
}


void remove_edge(Edge* e) {
    // Detach edge from neighbors
    splice(*e, e->oprev());
    splice(e->sym(), e->sym().oprev());
}


void delete_edge(Edge* e) {
    // Free memory -- edges and duals
    delete &e->left();
    delete &e->sym();
    delete &e->right();
    delete e;
};


void delete_edge_all(Edge* e) {
    delete_cell(e);
    remove_edge(e);
    delete_edge(e);
};


void splice(Edge& e0, Edge& e1) {
    Dual& d0 = e0.onext().rot();
    Dual& d1 = e1.onext().rot();

    Edge& a0 = e0.onext();
    Edge& a1 = e1.onext();
    Dual& b0 = d0.onext();
    Dual& b1 = d1.onext();

    e0.m_onext = &a1;
    e1.m_onext = &a0;
    d0.m_onext = &b1;
    d1.m_onext = &b0;
};


void swap(Edge& e) {
    assert(part_of_triangle(e) && "swap: edge must be part of a triangle");
    assert(part_of_triangle(e.sym()) && "swap: edge.sym must be part of a triangle");
    Edge& a = e.oprev();
    Edge& b = e.sym().oprev();
    splice(e, a);
    splice(e.sym(), b);
    splice(e, a.lnext());
    splice(e.sym(), b.lnext());
    e.m_org = a.m_dest;
    e.m_dest = b.m_dest;
}


bool part_of_triangle(const Edge& e0) {
    Edge& e1 = e0.lnext();
    Edge& e2 = e0.lnext().lnext();
    if (e0.org() == e1.org())
        return false;
    if (e1.org() == e2.org())
        return false;
    if (e2.org() == e0.org())
        return false;
    if (e0.dest() != e1.org())
        return false;
    if (e1.dest() != e2.org())
        return false;
    if (e2.dest() != e0.org())
        return false;
    return true;
}


bool is_free(const Edge& e) {
    bool free_edge = e == e.onext();
    if (free_edge)
        assert(e.sym() == e.sym().onext());
    else
        assert(e.sym() != e.sym().onext());
    return free_edge;
}


} // namespace quadedge


bool Edge::operator==(const Edge& e) const {
    return this->m_org->x == e.m_org->x
        && this->m_org->y == e.m_org->y
        && this->m_dest->x == e.m_dest->x
        && this->m_dest->y == e.m_dest->y;
}


bool Edge::operator!=(const Edge& e) const {
    return !operator==(e);
}


bool Edge::operator<(const Edge& e) const {
    if (this->org() == e.org())
        return this->dest() < e.dest();
    return this->org() < e.org();
}


Dual::~Dual() {
    if (m_data) {
        delete m_data;
        m_data = nullptr;
    }
}


void Dual::set_data(Triangle* cell) {
    Dual* d = this;
    size_t count = 0;
    while (d->m_data == nullptr) {
        d->m_data = cell;
        d = &d->onext();
        count++;
    }
    assert(count == 3);
}


std::ostream& operator<<(std::ostream& out, const Edge& e) {
    const Point* const p0 = e.m_org;
    const Point* const p1 = e.m_dest;
    out << "Edge: " << *p0 << " -> " << *p1;
    return out;
}


std::ostream& operator<<(std::ostream& out, const Dual& d) {
    const Point* const p0 = d.m_org;
    const Point* const p1 = d.m_dest;
    out << "Dual: " << *p0 << " -> " << *p1;
    return out;
}

} // namespace mesh

} // namespace umr
