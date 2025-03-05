#include "refinement_ruppert.hpp"

#include <numeric>

#include "inequalities.hpp"
#include "math.hpp"

namespace umr {

namespace mesh {

namespace algo {


static const double MIN_ANGLE_CONVERGE = 20.7;


RefinementRuppert::RefinementRuppert(double min_angle,
        density::DensityManager* density, bool use_terminator)
        : m_min_angle(min_angle), m_density(density),
          m_use_terminator(use_terminator) {

    if (m_min_angle > MIN_ANGLE_CONVERGE) {
        std::cout << "WARNING: The Ruppert refinement algorithm only guarantees "
            "convergence with a minimum angle of " << MIN_ANGLE_CONVERGE << ". "
            "No convergence is guaranteed with the input minimum angle "
            "of " << min_angle << ". Proceeding anyway.\n";
    }
}


void RefinementRuppert::initialize(MeshUnstructuredConstrained& mesh,
        std::vector<Real2>& points, std::list<parametric::IParametric*>& segments) {

    m_hmin = refine::calculate_minimum_length(points, segments);

    // Makes a database of source points
    for (auto [x, y] : points) {
        mesh.make_source_point(x, y);
        mesh.make_point(x, y);
    }

    // Makes a database of source edges (aka segments)
    std::vector<Edge*> segment_edges;
    for (auto segment : segments) {
        std::vector<Edge*> ev = mesh.make_source_edge(segment);
        segment_edges.insert(segment_edges.end(), ev.begin(), ev.end());
    }

    // Triangulates the initial mesh, using points and segments
    // from the source databases whenever possible
    mesh.add_external_points();
    delaunay::divide_and_conquer(mesh);

    // Any segment not in the mesh must be inserted now
    Edge* ei = mesh.get_initial_edge();
    for (Edge* e : segment_edges) {
        if (quadedge::is_free(*e))
            delaunay::insert_edge(mesh, *e, *ei);
        ei = e;
    }

    // Divide and Conquer can insert and then delete a segment!
    // Empty the deleted edge tracker now
    while(mesh.get_next_deleted_source_edge())
        continue;

    // Finalize the initial mesh
    mesh.remove_external_edges();
    mesh.fill_triangles();

    // Fill the initial encroachment queue
    for (auto se : mesh.get_source_edges()) {
        Edge* e = se.first;
        Triangle* t = e->left().get_data();
        if (t && is_encroached(mesh, *e))
            m_encroached_queue.insert(*e);

        e = &e->sym();
        t = e->left().get_data();
        if (t && is_encroached(mesh, *e))
            m_encroached_queue.insert(*e);
    }
}


bool RefinementRuppert::iteration(MeshUnstructuredConstrained& mesh) {

    while (!m_encroached_queue.empty()) {
        Edge e = m_encroached_queue.top();
        Edge* esrc = mesh.get_source(&e);
        m_encroached_queue.pop();

        Triangle* triangle_on_left = esrc->left().get_data();
        Triangle* triangle_on_right = esrc->right().get_data();

        double min_length_left = triangle_on_left
            ? triangle_on_left->get_min_length()
            : std::numeric_limits<double>::max();
        double min_length_right = triangle_on_right
            ? triangle_on_right->get_min_length()
            : std::numeric_limits<double>::max();

        double min_angle_left = triangle_on_left
            ? triangle_on_left->get_min_angle()
            : std::numeric_limits<double>::max();
        double min_angle_right = triangle_on_right
            ? triangle_on_right->get_min_angle()
            : std::numeric_limits<double>::max();

        double min_length = std::min(min_length_left, min_length_right);
        double min_angle = std::min(min_angle_left, min_angle_right);

        bool angle_small_enough = inequalities::is_lt(min_angle, m_min_angle);
        bool not_on_boundary = triangle_on_left && triangle_on_right;
        if (angle_small_enough || not_on_boundary) {
            bool success = split_encroached_segment(mesh, *esrc, min_length);
            if (success)
                return false;
        }
    }

    algo::ITrianglePQ& triangle_pq = mesh.get_triangle_pq();
    if (triangle_pq.empty())
        return true;

    Triangle* triangle = nullptr;
    Edge* triangle_edge = nullptr;
    while (!triangle_pq.empty()) {
        Triangle* t = triangle_pq.top().first;
        triangle_edge = triangle_pq.top().second;

        bool passes_size_test =
            refine::convergence_test_density(*t, *m_density);
        if (!passes_size_test) {
            triangle = t;
            break;
        }

        bool passes_shape_test =
            refine::convergence_test_angle(*t, m_min_angle);
        if (!passes_shape_test) {
            triangle = t;
            break;
        }

        triangle_pq.pop();
    }

    // Convergence tests succeeded for all triangles in the queue
    if (triangle == nullptr)
        return true;

    Point proposed_point = triangle->get_circumcenter();
    double prev_length = triangle->get_min_length();
    double triangle_radius = triangle->get_radius();
    double min_length = std::min(prev_length,
            m_density->evaluate(proposed_point.x, proposed_point.y));

    // get edge nearest proposed point, if possible
    Edge* ei = &refine::locate_or_closest_edge(mesh, proposed_point, *triangle_edge);

    // case where proposed point it out of bounds
    if (mesh.is_source(ei)) {
        bool contains_circumcenter = math::left_of(proposed_point, *ei);
        contains_circumcenter &= math::left_of(proposed_point, ei->lnext());
        contains_circumcenter &= math::left_of(proposed_point, ei->lnext().lnext());
        if (!contains_circumcenter)
            m_encroached_queue.insert(*ei);
    }

    // case where proposed point encroaches segments
    if (m_encroached_queue.size() == 0)
        m_encroached_queue = refine::fill_encroached_edge_queue(mesh,
                *triangle_edge, proposed_point, triangle_radius);

    bool success = false;
    if (m_encroached_queue.empty()) {
        min_length = std::min(min_length, triangle_radius);
        refine::insert_point(mesh, proposed_point, *ei, min_length);
        success = true;
    } else {
        // split ONE encroached segment; the rest are split
        // on subsequent iterations
        while (!success && !m_encroached_queue.empty()) {
            Edge e = m_encroached_queue.top();
            Edge* esrc = mesh.get_source(&e);
            m_encroached_queue.pop();
            // TODO: should min_length be recalibrated with the esrc midpoint?
            success = split_encroached_segment(mesh, *esrc, min_length);
        }
    }

    // Delete the triangle manually if no splitting occured
    if (!success)
        triangle_pq.pop();

    return triangle_pq.empty();
}


bool RefinementRuppert::is_encroached(MeshUnstructuredConstrained& mesh, Edge& edge) {

    Point& p = edge.lnext().dest();
    return is_encroached(mesh, edge, p);
}


bool RefinementRuppert::is_encroached(MeshUnstructuredConstrained& mesh,
        Edge& edge, Point& point) {

    if (!mesh.is_source_point(&point))
        return false;
    return math::in_diametral_circle(edge.org(), edge.dest(), point);
}


bool RefinementRuppert::split_encroached_segment(
        MeshUnstructuredConstrained& mesh, Edge& edge, double min_length) {

    bool split_left = edge.left().get_data();
    bool split_right = edge.right().get_data();
    assert(split_left || split_right);

    double xmid = (edge.org().x + edge.dest().x) / 2;
    double ymid = (edge.org().y + edge.dest().y) / 2;
    Point proposed_point(xmid, ymid);
    refine::ClusterMetrics cluster =
        refine::get_segment_cluster_metrics(mesh, edge, m_min_angle);

    std::vector<Edge*> new_edges;
    bool use_terminator = inequalities::is_le(cluster.angle, 60);
    if (m_use_terminator && use_terminator) {
        new_edges = refine::split_encroached_edge_terminator(mesh, edge, min_length, cluster);
        if (new_edges.size() == 0)
            return false;
    } else {
        new_edges = mesh.split_edge(edge, 2);
    }

    Point& midpoint = new_edges[0]->dest();

    if (split_left) {
        refine::make_cavity(mesh, *new_edges[1], midpoint);
        refine::split_cavity(mesh, new_edges[0]->sym());
    }

    if (split_right) {
        refine::make_cavity(mesh, new_edges[0]->sym(), midpoint);
        refine::split_cavity(mesh, *new_edges[1]);
    }

    queue_encroached(mesh, *new_edges[0]);
    mesh.fill_triangles(*new_edges[0]);

    return true;
}


// fill any newly encroached segments splitting but BEFORE constructing triangles
void RefinementRuppert::queue_encroached(MeshUnstructuredConstrained& mesh, Edge& edge) {
    bool none_left = edge.left().get_data() == nullptr;
    bool none_right = edge.right().get_data() == nullptr;
    assert(none_left || none_right);
    bool triangle_left = quadedge::part_of_triangle(edge);
    bool triangle_right = quadedge::part_of_triangle(edge.sym());
    assert(triangle_left || triangle_right);

    std::queue<Edge*> queue;
    std::unordered_set<Edge*> visited;
    if (triangle_left && none_left) {
        queue.push(&edge);
        visited.insert(&edge);
    }
    if (triangle_right && none_right) {
        queue.push(&edge.sym());
        visited.insert(&edge.sym());
    }

    while (!queue.empty()) {
        Edge* e0 = queue.front();
        queue.pop();

        Edge* ei = e0;
        Edge* efin = ei;
        do {
            if (mesh.is_source(ei))
                if (is_encroached(mesh, *ei))
                    m_encroached_queue.insert(*ei);
            ei = &ei->lnext();
        } while (ei != efin);

        Edge* e1 = &e0->lnext().sym();
        Edge* e2 = &e0->lnext().lnext().sym();
        if (e1->left().get_data() == nullptr && quadedge::part_of_triangle(*e1)) {
            if (!visited.contains(e1)) {
                queue.push(e1);
                visited.insert(e1);
                visited.insert(&e1->sym());
            }
        }
        if (e2->left().get_data() == nullptr && quadedge::part_of_triangle(*e2)) {
            if (!visited.contains(e2)) {
                queue.push(e2);
                visited.insert(e2);
                visited.insert(&e2->sym());
            }
        }
    }
}


} // namespace algo

} // namespace mesh

} // namespace umr
