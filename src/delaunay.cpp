#include "delaunay.hpp"

#include <cassert>
#include <iostream>

#include "math.hpp"


namespace umr {

namespace mesh {

namespace delaunay {


/////////////////////////////////
// Local function declarations //
/////////////////////////////////

std::pair<Edge*, Edge*> _divide_and_conquer_iteration(MeshUnstructured& mesh,
        std::vector<Point*>& points,
        size_t imin, size_t imax);

void _rising_bubble(MeshUnstructured& mesh, Edge* basel);

void _giftwrap(MeshUnstructured& mesh, Edge& edge);

/////////////////////////////////
/////////////////////////////////


Edge& locate(const Point& point, Edge& edge) {
    size_t iteration = 0;
    Edge* e = &edge;
    while (true) {
        if (iteration++ >= 10000) {
            std::cout << "locate: Edge not found after 10000 iterations\n";
            exit(-1);
        }

        if (e->org() == point || e->dest() == point) {
            break;
        } else if (math::right_of(point, *e)) {
            e = &e->sym();
        } else if (math::left_of(point, e->onext())) {
            e = &e->onext();
        } else if (math::left_of(point, e->dprev())) {
            e = &e->dprev();
        } else {
            if (math::on_line(point, e->onext()))
                e = &e->onext();
            else if (math::on_line(point, e->dprev()))
                e = &e->dprev();
            break;
        }
    }

    // Make sure the point is either on a triangle corner
    // or within the triangle
    Point& pa = e->org();
    Point& pb = e->dest();
    Point& pc = e->lnext().dest();
    if (pc == point)
        e = &e->lprev();
    else if (pb == point)
        e = &e->lnext();
    else if (pa != point)
        assert(math::in_circle(pa, pb, pc, point) \
                && "locate: target point " \
                "is not within the found edge's triangle");

    return *e;
}


// Effectively a "make_cavity" for inserting an edge; when done
// calling giftwrap as a "split_cavity".
void insert_edge(MeshUnstructured& mesh, Edge& new_edge, Edge& edge) {
    bool already_joined_to_something = new_edge.onext() != new_edge;
    already_joined_to_something |= new_edge.sym().onext() != new_edge.sym();
    assert(!already_joined_to_something
            && "insert_edge: edge to be inserted must not be joined to anything");

    // Find edges with correct points as org
    Edge* ei = &locate(new_edge.org(), edge);
    Edge* ef = &locate(new_edge.dest(), *ei);
    while (ei->org() != new_edge.org())
        ei = &ei->lnext();
    while (ef->org() != new_edge.dest())
        ef = &ef->lnext();

    // Rotate edges as needed for splicing
    if (ei != &new_edge) {
        while (!math::right_of(new_edge.dest(), *ei))
            ei = &ei->oprev();
        while (math::right_of(new_edge.dest(), *ei))
            ei = &ei->oprev();
    }

    if (ei != &new_edge) {
        while (!math::right_of(new_edge.org(), *ef))
            ef = &ef->oprev();
        while (math::right_of(new_edge.org(), *ef))
            ef = &ef->oprev();
    }

    if (*ei == new_edge) {
        // Edge already exists; remove it and replace it with the new_edge
        // Important in case new_edge pointer is stored elsewhere (e.g. sources)
        assert(*ef == new_edge.sym());
        ei = &ei->oprev();
        ef = &ef->oprev();
        mesh.delete_edge(&ei->onext());
        mesh.insert_edge(&new_edge);
        quadedge::splice(*ei, new_edge);
        quadedge::splice(*ef, new_edge.sym());
    } else {
        // Edge does not exist; must add it and refine with giftwrap as necessary
        assert(*ei != new_edge);
        assert(*ef != new_edge.sym());
        Edge* ej = &ei->lnext();
        while (ej->dest() != new_edge.dest()) {
            if (math::intersect(*ej, new_edge)) {
                ej = &ej->oprev();
                mesh.delete_edge(&ej->onext());
            } else {
                ej = &ej->lnext();
            }
        }

        mesh.insert_edge(&new_edge);
        quadedge::splice(*ei, new_edge);
        quadedge::splice(*ef, new_edge.sym());

        _giftwrap(mesh, new_edge);
        _giftwrap(mesh, new_edge.sym());
    }
}


Edge& make_cavity(MeshUnstructured& mesh, const Point& point, Edge& edge) {
    Edge* er = &edge;
    Edge* e = &edge;

    // Delete the triangle cell now, if any, in case no edges
    // get removed.
    if (edge.left().get_data())
        mesh.delete_triangle(&edge);

    while (true) {
        // Ensure edge is not on the domain boundary
        if (e->sym() != e->lnext()) { // ensure not a dangling edge
            if (!quadedge::part_of_triangle(e->sym())) {
                e = &e->lnext();
                if (e == er)
                    break;
                continue;
            }
        }

        bool valid_triangle = quadedge::part_of_triangle(e->sym());
        if (valid_triangle) {
            bool contains = math::in_circle(e->sym().org(), e->sym().dest(),
                    e->sym().lnext().dest(), point);
            if (contains) {
                er = &e->oprev();
                mesh.delete_edge(e);
                e = er;
                continue;
            }
        }

        e = &e->lnext();
        if (e == er)
            break;
    }

    return *er;
}


// Messy implementation using Divide and Conquer to make a new mesh,
// then insert the new mesh into the cavity.
void split_cavity(MeshUnstructured& mesh, Edge& edge, std::vector<Point*> new_points) {

    // Construct temporary mesh to triangulate the cavity
    MeshUnstructuredConstrained m;

    // Add cavity boundary points to the mesh
    Edge* efin = &edge;
    Edge* e = &edge;
    do {
        m.insert_point(&e->org());
        e = &e->lnext();
    } while (e != efin);

    // Add external points to the mesh
    m.add_external_points();

    // Add new points to be placed within the cavity
    for (auto p : new_points)
        m.insert_point(p);

    // Triangularize the cavity
    divide_and_conquer(m);

    // Not all cavity edges are guaranteed to exist.
    // Add them in as needed.
    assert(e == efin);
    do {
        Point& p0 = m.make_point(e->org().x, e->org().y);
        Point& p1 = m.make_point(e->dest().x, e->dest().y);
        Edge* en;
        if (!m.has_edge(e)) {
            // make edge externally (NOT using the mesh, then insert)
            en = &quadedge::make_edge(p0, p1);
            Edge* einit = m.get_initial_edge();
            insert_edge(m, *en, *einit);
        } else { // more like "get" than "make"
            en = &m.make_edge(p0, p1);
        }
        // SHOULD never need to split, currently "cheating" with nullptr
        m.add_source(en, nullptr);
        e = &e->lnext();
    } while (e != efin);

    // Locate an edge on the cavity for merging
    Edge* em = m.get_initial_edge();
    em = &locate(e->org(), *em);
    while (em->org() != e->org())
        em = &em->lnext();
    while (*em != *e)
        em = &em->onext();

    // Delete external points
    m.remove_external_edges();

    // Merge meshes -- bounds of m with cavity of mesh
    while (true) {
        Edge* enext = &e->lnext();
        Edge* emnext = &em->sym().onext();
        // Join e and em together in a loop, THEN delete em
        // Should prevent issues with cavity at external boundary
        quadedge::splice(*e, em->sym().lnext());
        if (e->lnext() == em->sym()) {
            m.delete_edge(em);
            break;
        }
        m.delete_edge(em);
        e = enext;
        em = emnext;
    }

    // Insert the remaining edges and points into the mesh
    for (Edge* e : m.get_edges_set())
        mesh.insert_edge(e);

    for (Point* p : m.get_points_set())
         mesh.insert_point(p);

    // All are now added; clear them from the mesh so
    // no object memory gets freed
    m.clear();

}


Edge& insert_point(MeshUnstructured& mesh, Point& point, Edge& edge) {
    Edge* e = &edge;
    if (e->org() == point)
        return *e;
    if (e->dest() == point)
        return e->sym();

    if (math::on_line(point, *e)) {
        Edge* t = &e->oprev();
        mesh.delete_edge(e);
        e = t;
    }
    Edge* base = &mesh.make_edge(e->org(), point);
    quadedge::splice(*base, *e);
    Point& first = e->org();
    while (true) {
        base = &mesh.connect_edges(*e, base->sym());
        e = &base->oprev();
        if (e->dest() == first)
            break;
    }

    assert(base->dest() == point && "base dest is not p");
    Edge* eret = &base->sym();

    e = &base->oprev();
    while (true) {
        Edge* t = &e->oprev();
        if (math::right_of(point, *e)
                && math::in_circle(e->org(), t->dest(), e->dest(), point)) {
            if (eret == e || eret == &e->sym())
                eret = &eret->onext();
            quadedge::swap(*e);
            e = t;
        } else if (e->org() == first) {
            break;
        } else {
            e = &e->onext().lprev();
        }
    }

    assert(eret->org() == point && "eret org is not point");
    return *eret;
}


void bowyer_watson(MeshUnstructured& mesh) {
    auto& points = mesh.get_points_set();
    const std::vector<Point*>& external = mesh.add_external_points();
    Edge& e0 = mesh.make_edge(*external[0], *external[1]);
    Edge& e1 = mesh.extend_edge(e0, *external[3]);
    mesh.connect_edges(e1, e0);
    Edge& e2 = mesh.extend_edge(e1.sym(), *external[2]);
    mesh.connect_edges(e2, e1.sym());

    Edge* e = &e0;
    for (Point* p : points) {
        e = &locate(*p, *e);
        e = &make_cavity(mesh, *p, *e);
        e = &insert_point(mesh, *p, *e);
    }

    mesh.fill_triangles();
}


void divide_and_conquer(MeshUnstructured& mesh) {
    const double tol = 1e-7;

    auto& pts = mesh.get_points_set();
    std::vector<Point*> points(pts.begin(), pts.end());
    std::sort(points.begin(), points.end(), [&](Point* const p0, Point* const p1) {
        if (std::abs(p0->x - p1->x) <= tol)
            return p0->y < p1->y - tol;
        return p0->x < p1->x - tol;
    });

    size_t n_pts = points.size();
    _divide_and_conquer_iteration(mesh, points, 0, n_pts);
}


////////////////////////////////
// Local function definitions //
////////////////////////////////

void _rising_bubble(MeshUnstructured& mesh, Edge* basel) {
    while (true) {
        Edge* lcand = &basel->sym().onext();
        bool valid_lcand = math::right_of(lcand->dest(), *basel);
        if (valid_lcand) {
            while (math::in_circle(basel->dest(), basel->org(),
                            lcand->dest(), lcand->onext().dest())) {
                Edge& t = lcand->onext();
                mesh.delete_edge(lcand);
                lcand = &t;
            }
        }

        Edge* rcand = &basel->oprev();
        bool valid_rcand = math::right_of(rcand->dest(), *basel);
        if (valid_rcand) {
            while (math::in_circle(basel->dest(), basel->org(),
                            rcand->dest(), rcand->oprev().dest())) {
                Edge& t = rcand->oprev();
                mesh.delete_edge(rcand);
                rcand = &t;
            }
        }

        valid_lcand = math::right_of(lcand->dest(), *basel);
        valid_rcand = math::right_of(rcand->dest(), *basel);
        if (!valid_lcand && !valid_rcand) {
            break;
        } else if (valid_lcand && !valid_rcand) {
            basel = &mesh.connect_edges(basel->sym(), lcand->sym());
        } else if (!valid_lcand && valid_rcand) {
            basel = &mesh.connect_edges(*rcand, basel->sym());
        } else {
            bool left_violates_delaunay = math::in_circle(
                    lcand->dest(), lcand->org(), rcand->org(), rcand->dest());
            if (left_violates_delaunay) {
                basel = &mesh.connect_edges(*rcand, basel->sym());
            } else {
                basel = &mesh.connect_edges(basel->sym(), lcand->sym());
            }
        }
    }
}


std::pair<Edge*, Edge*> _divide_and_conquer_iteration(MeshUnstructured& mesh,
        std::vector<Point*>& points,
        size_t imin, size_t imax) {
    size_t num = imax - imin;
    assert(num > 1 && "_divide_and_conquer_iteration must have >= 2 points");

    if (num == 2) {
        Edge& e = mesh.make_edge(*points[imin], *points[imin+1]);
        return std::make_pair(&e, &e.sym());
    }

    if (num == 3) {
        Edge& e0 = mesh.make_edge(*points[imin], *points[imin+1]);
        Edge& e1 = mesh.make_edge(*points[imin+1], *points[imin+2]);
        quadedge::splice(e0.sym(), e1);
        if (math::left_of(e1.dest(), e0)) {
            mesh.connect_edges(e1, e0);
            return std::make_pair(&e0, &e1.sym());
        }
        if (math::right_of(e1.dest(), e0)) {
            Edge& e2 = mesh.connect_edges(e0.sym(), e1.sym());
            return std::make_pair(&e2, &e2.sym());
        }
        return std::make_pair(&e0, &e1.sym());
    }

    size_t imid = (imin + imax) >> 1;
    std::pair<Edge*, Edge*> el =
        _divide_and_conquer_iteration(mesh, points, imin, imid);
    std::pair<Edge*, Edge*> er =
        _divide_and_conquer_iteration(mesh, points, imid, imax);
    Edge* ldo = el.first;
    Edge* ldi = el.second;
    Edge* rdi = er.first;
    Edge* rdo = er.second;
    while (true) {
        if (math::left_of(rdi->org(), *ldi))
            ldi = &ldi->lnext();
        else if (math::right_of(ldi->org(), *rdi))
            rdi = &rdi->rprev();
        else
            break;
    }

    Edge* basel = &mesh.connect_edges(rdi->sym(), *ldi);
    if (ldi->org() == ldo->org())
        ldo = &basel->sym();
    if (rdi->org() == rdo->org())
        rdo = basel;

    _rising_bubble(mesh, basel);

    return std::make_pair(ldo, rdo);
}


void _giftwrap(MeshUnstructured& mesh, Edge& edge) {
    if (quadedge::part_of_triangle(edge))
        return;

    Edge* e = &edge.lnext();
    Edge* ref = nullptr;
    while (e->dest() != edge.org()) {
        if (!ref || math::in_circle(edge.org(), edge.dest(), ref->dest(), e->dest()))
            ref = e;
        e = &e->lnext();
    }

    if (ref) {
        Edge& left_bound = edge.lprev();
        Edge& right_bound = edge.lnext();
        Edge* before = ref;
        Edge* after = &ref->lnext();
        if (*before != edge.lnext()) {
            assert(math::is_ccw(before->org(), before->dest(), edge.lnext().org()));
            mesh.connect_edges(*before, edge.lnext());
        }
        if (*after != edge.lprev()) {
            assert(math::is_ccw(edge.lprev().org(), edge.lprev().dest(), after->org()));
            mesh.connect_edges(edge.lprev(), *after);
        }
        assert(quadedge::part_of_triangle(edge));
        assert(quadedge::part_of_triangle(edge.lnext()));
        assert(quadedge::part_of_triangle(edge.lnext().lnext()));
        if (edge.lnext() != right_bound)
            _giftwrap(mesh, edge.lnext().sym());
        if (edge.lprev() != left_bound)
            _giftwrap(mesh, edge.lnext().lnext().sym());
    }
}


////////////////////////////////
////////////////////////////////


} // namespace delaunay

} // namespace mesh

} // namespace umr
