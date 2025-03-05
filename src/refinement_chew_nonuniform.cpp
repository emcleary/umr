#include "refinement_chew_nonuniform.hpp"

#include "inequalities.hpp"
#include "math.hpp"


namespace umr {

namespace mesh {

namespace algo {


static const double MIN_ANGLE_CONVERGE = 26.5;


RefinementChewNonuniform::RefinementChewNonuniform(double min_angle,
        density::DensityManager* density, bool use_terminator)
        : m_min_angle(min_angle), m_density(density),
          m_use_terminator(use_terminator) {

    if (m_min_angle > MIN_ANGLE_CONVERGE) {
        std::cout << "WARNING: The Chew-Nonuniform refinement algorithm "
            "only guarantees convergence with a minimum angle "
            "of " << MIN_ANGLE_CONVERGE << " deg. No convergence is guaranteed "
            "with the input minimum angle of " << min_angle << " deg. "
            "Proceeding anyway.\n";
    }
}


void RefinementChewNonuniform::initialize(MeshUnstructuredConstrained& mesh,
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
}


bool RefinementChewNonuniform::iteration(MeshUnstructuredConstrained& mesh) {
    algo::ITrianglePQ& triangle_pq = mesh.get_triangle_pq();

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

    // find edge nearest triangle center, starting from an edge of the triangle
    // if circumcenter is on a line, this function will return the edge on the line!
    Edge* ei = &refine::locate_or_closest_edge(mesh, triangle->get_circumcenter(), *triangle_edge);
    bool on_source = mesh.is_source(ei);
    bool contains_circumcenter = math::left_of(triangle->get_circumcenter(), *ei);
    contains_circumcenter &= math::left_of(triangle->get_circumcenter(), ei->lnext());
    contains_circumcenter &= math::left_of(triangle->get_circumcenter(), ei->lnext().lnext());

    // propose edge midpoint or triangle circumcenter?
    Point proposed_point = on_source && !contains_circumcenter
        ? Point((ei->org().x + ei->dest().x) / 2, (ei->org().y + ei->dest().y) / 2)
        : Point(triangle->get_circumcenter());

    double prev_length = triangle->get_min_length();
    double triangle_radius = triangle->get_radius();

    // Terminator algorithm
    algo::EdgeQueue encroached_queue = refine::fill_encroached_edge_queue(
            mesh, *triangle_edge, proposed_point, triangle_radius);
    bool small_angle_cluster_found = false;
    bool success = false;
    if (m_use_terminator) {
        while (!success && !encroached_queue.empty()) {
            Edge e = encroached_queue.top();
            Edge* encroached_edge = mesh.get_source(&e);
            encroached_queue.pop();
            double xmid = (encroached_edge->org().x + encroached_edge->dest().x) / 2;
            double ymid = (encroached_edge->org().y + encroached_edge->dest().y) / 2;
            double min_length = std::min(m_hmin, m_density->evaluate(xmid, ymid));
            min_length = std::min(prev_length, min_length);
            refine::ClusterMetrics cluster = refine::get_segment_cluster_metrics(
                    mesh, *encroached_edge, m_min_angle);
            if (inequalities::is_le(cluster.angle, m_min_angle)) {
                success |= split_encroached_segment(
                        mesh, *encroached_edge, min_length, 2, true);
                small_angle_cluster_found = true;
            }
        }
    }

    if (!small_angle_cluster_found) {
        // Chew's nonuniform algorithm
        bool on_line = math::on_line(triangle->get_circumcenter(), *ei);
        bool on_line_but_not_source = on_line && !on_source;
        if (contains_circumcenter || on_line_but_not_source) {
            // add the triangle's circumcenter to the mesh
            const Point& pc = triangle->get_circumcenter();
            double min_length = std::min(m_hmin, m_density->evaluate(pc.x, pc.y));
            // including the prev_length here can prevent going back
            // and forth in a cycle; keeping the radius of the
            // inserted point smaller than its surroundings, so only
            // points coincidentally too close will be deleted
            min_length = std::min(prev_length, min_length);
            min_length = std::min(triangle_radius, min_length);
            refine::insert_point(mesh, pc, *ei, min_length);
        } else {
            // split segment
            assert(mesh.is_source(ei));
            double xmid = (ei->org().x + ei->dest().x) / 2;
            double ymid = (ei->org().y + ei->dest().y) / 2;
            double min_length = std::min(m_hmin, m_density->evaluate(xmid, ymid));
            double edge_length = ei->org().distance_to(ei->dest());

            bool split_to_three = inequalities::is_lt(edge_length, 4 * min_length)
                && inequalities::is_gt(edge_length, 2 * sqrt(3) * min_length);
            int num_segments = split_to_three ? 3 : 2;

            min_length = std::min(prev_length, min_length);
            split_encroached_segment(mesh, *ei, min_length, num_segments);
        }
    } else if (!success) {
        // no splitting occured so the triangle must be deleted manually
        triangle_pq.pop();
    }

    return triangle_pq.empty();
}


bool RefinementChewNonuniform::split_encroached_segment(
        MeshUnstructuredConstrained& mesh, Edge& edge,
        double min_length, int num_segments, bool terminator_only) {

    bool split_left = edge.left().get_data();
    bool split_right = edge.right().get_data();
    assert(split_left || split_right);

    double xmid = (edge.org().x + edge.dest().x) / 2;
    double ymid = (edge.org().y + edge.dest().y) / 2;
    Point proposed_point(xmid, ymid);
    refine::ClusterMetrics cluster =
        refine::get_segment_cluster_metrics(mesh, edge, m_min_angle);

    std::vector<Edge*> new_edges;
    bool use_terminator = terminator_only || inequalities::is_le(cluster.angle, 60);
    if (m_use_terminator && use_terminator) {
        new_edges = refine::split_encroached_edge_terminator(
                mesh, edge, min_length, cluster, num_segments);
        if (new_edges.size() == 0)
            return false;
    } else {
        new_edges = mesh.split_edge(edge, num_segments);
    }

    double min_new_segment_length = new_edges[0]->org().distance_to(new_edges[0]->dest());
    for (auto it = new_edges.begin()+1; it != new_edges.end(); ++it) {
        double segment_length = (*it)->org().distance_to((*it)->dest());
        min_new_segment_length = std::min(min_new_segment_length, segment_length);
    }

    if (split_left)
        refine::make_and_split_cavity(mesh, new_edges[1]->sym(),
                min_new_segment_length, new_edges[1]->org());
    if (split_right)
        refine::make_and_split_cavity(mesh, *new_edges[1],
                min_new_segment_length, new_edges[1]->org());

    mesh.fill_triangles(*new_edges[0]);

    return true;
}


} // namespace algo

} // namespace mesh

} // namespace umr
