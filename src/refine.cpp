#include "refine.hpp"

#include <cassert>

#include <memory>
#include <unordered_set>

#include "math.hpp"
#include "delaunay.hpp"
#include "integrators.hpp"
#include "inequalities.hpp"


namespace umr {

namespace mesh {

namespace refine {


using PairFloatParam = std::pair<double, parametric::IParametric*>;


std::pair<Point, Point> mesh_point_limits(const MeshUnstructuredConstrained& mesh) {

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();

    for (auto p : mesh.get_source_points()) {
        xmin = std::min(xmin, p->x);
        xmax = std::max(xmax, p->x);
        ymin = std::min(ymin, p->y);
        ymax = std::max(ymax, p->y);
    }

    Point pmin(xmin, ymin);
    Point pmax(xmax, ymax);
    return {pmin, pmax};
}


// used for all algorithms
void mesh_quality(const MeshUnstructuredConstrained& mesh) {
    double max_radius = 0;
    double min_radius = std::numeric_limits<double>::max();
    double max_angle = 0;
    double min_angle = 180;
    double min_angle_geometry = 180;
    double max_length = 0;
    double min_length = std::numeric_limits<double>::max();
    std::unordered_set<Triangle*> triangles;

    auto update_quality = [&](Edge& e) {
        Triangle* t = e.left().get_data();
        assert(t);

        triangles.insert(t);

        Point& pa = e.org();
        Point& pb = e.dest();
        Point& pc = e.lnext().dest();
        double ab = pa.distance_to(pb);
        double bc = pb.distance_to(pc);
        double ca = pc.distance_to(pa);
        double angle = math::loc_angle_deg(ab, bc, ca);
        if (mesh.is_source(&e) && mesh.is_source(&e.lnext())) {
            min_angle_geometry = std::min(min_angle_geometry, angle);
        } else {
            max_length = std::max(max_length, t->get_max_length());
            min_length = std::min(min_length, t->get_min_length());
            max_radius = std::max(max_radius, t->get_radius());
            min_radius = std::min(min_radius, t->get_radius());
            max_angle = std::max(max_angle, angle);
            min_angle = std::min(min_angle, angle);
        }
    };

    for (Edge* e : mesh.get_edges_set()) {
        update_quality(*e);
        update_quality(e->sym());
    }

    std::cout << "Mesh quality\n";
    if (min_angle < 180) {
        // Only print this data if there are interior triangles.
        // Corner only triangles can corrupt this data
        // i.e. input min angle of the geometry should be
        // kept separate from the min angle of the refined triangles.
        std::cout << "   min radius = " << min_radius << '\n';
        std::cout << "   max radius = " << max_radius << '\n';;
        std::cout << "   min length = " << min_length << '\n';
        std::cout << "   max length = " << max_length << '\n';
        std::cout << "   min angle = " << min_angle << " deg\n";
        std::cout << "   max angle = " << max_angle << " deg\n";
    }
    std::cout << "   min angle geometry = " << min_angle_geometry << " deg\n";
    std::cout << "   number of triangles = " << triangles.size() << '\n';
    std::cout << '\n';
}


// used by chew uniform and tests 4-6
void uniform_grid(MeshUnstructured& mesh, const Point& pmin, const Point& pmax, double hmin, int nextra) {
    double length = hmin * std::sqrt(3);
    double xmin = pmin.x - nextra * length;
    double xmax = pmax.x + nextra * length;
    double ymin = pmin.y - nextra * length;
    double ymax = pmax.y + nextra * length;
    double dx = length;
    double dy = length * std::sqrt(3) / 2;
    const size_t n_pts_x = 1 + std::ceil((xmax - xmin) / dx);
    const size_t n_pts_y = 1 + std::ceil((ymax - ymin) / dy);
    assert(n_pts_x >= 2);
    assert(n_pts_y >= 2);

    std::cout << "UNIFORM GRID\n";
    std::cout << "   x in [" << xmin << ", " << xmax << "]\n";
    std::cout << "   y in [" << ymin << ", " << ymax << "]\n";
    std::cout << '\n';

    auto add_line = [&](double x, double y, const size_t npts) {
        Point* p0 = &mesh.make_point(x, y);
        x += dx;
        Point* p1 = &mesh.make_point(x, y);
        Edge* e = &mesh.make_edge(*p0, *p1);
        for (size_t i = 0; i < npts - 2; i++) {
            x += dx;
            Point* p = &mesh.make_point(x, y);
            e = &mesh.extend_edge(*e, *p);
        }
        return &e->sym();
    };

    auto join_lines = [&](Edge* e0, Edge* e1) {
        while (true) {
            mesh.connect_edges(*e0, *e1);
            if (e0->dest().x < e1->org().x) {
                e0 = &e1->onext();
                e1 = &e0->onext().sym();
                if (e1->org().y < e1->dest().y)
                    break;
            } else {
                e1 = &e1->onext();
                e0 = &e1->lnext();
                if (e0->org().x < e0->dest().x)
                    break;
            }
        }
    };

    double y = ymin;
    Edge* e0 = add_line(xmin, y, n_pts_x);

    for (size_t i = 0; i < n_pts_y - 1; i += 2) {
        y += dy;
        Edge* e1 = add_line(xmin + dx/2, y, n_pts_x-1);
        y += dy;
        Edge* e2 = add_line(xmin, y, n_pts_x);

        Edge* eb = &mesh.connect_edges(e0->sym(), *e1);
        join_lines(eb, &e0->sym());

        eb = &mesh.connect_edges(e1->sym(), *e2);
        join_lines(e2, eb);

        e0 = e2;
    }

    mesh.fill_triangles();
}


// used by chew uniform only
//
// TODO: this function, or at least where it gets used, can likely be
// replaced by new make/split cavity functions
Edge* refine_triangles(MeshUnstructuredConstrained& mesh, Edge& edge, double hmin) {
    mesh.fill_triangles(edge);
    mesh.fill_triangles(edge.sym());
    algo::ITrianglePQ& triangle_pq = mesh.get_triangle_pq();
    Edge* ei = &edge;
    while (!triangle_pq.empty()) {
        Triangle* t = triangle_pq.top().first;
        if (inequalities::is_gt(t->get_radius(), hmin)) {
            Point* p = new Point(t->get_circumcenter());
            ei = &delaunay::locate(*p, *ei);
            ei = &delaunay::make_cavity(mesh, *p, *ei);
            std::vector<Point*> new_points;
            new_points.push_back(p);
            delaunay::split_cavity(mesh, *ei, new_points);
            mesh.fill_triangles(*ei);
            assert(!triangle_pq.is_active(t));
        } else {
            triangle_pq.pop();
        }
    }
    return ei;
}


// used by chew uniform and chew nonuniform
Edge* insert_point(MeshUnstructuredConstrained& mesh, const Point& point, Edge& edge, double hmin) {

    Edge* ei = &edge;
    bool has_point = mesh.has_point(&point, ei);
    if (!has_point) {

        ei = &delaunay::locate(point, *ei);

        Point& p = mesh.make_point(point.x, point.y);
        if (math::on_line(p, *ei)) {
            Edge* etmp = ei;
            ei = &ei->oprev();
            mesh.delete_edge(etmp);
        } else {
            assert(ei->left().get_data());
            mesh.delete_triangle(ei);
            // Ensure point is contained within current triangle
            assert(math::left_of(p, *ei));
            assert(math::left_of(p, ei->lnext()));
            assert(math::left_of(p, ei->lnext().lnext()));
        }

        // Insert the new point into the mesh
        Edge* eprev = &mesh.make_edge(ei->org(), p);
        quadedge::splice(*ei, *eprev);
        Edge* efin = &eprev->lprev();
        while (ei != efin) {
            eprev = &mesh.connect_edges(*ei, eprev->sym());
            ei = &eprev->oprev();
        }
        mesh.fill_triangles(*ei);

        ei = &ei->lprev();

        // Triangularize cells around the inserted point.
        // Must keep track of each edge used while looping
        // around the point since each edge can get swapped.
        std::unordered_set<Point, PointHash> visited;
        while (true) {
            bool unvisited_edge = visited.insert(ei->dest()).second;
            if (!unvisited_edge)
                break;
            make_and_split_cavity(mesh, ei->sym(), hmin, p);
            mesh.fill_triangles(*ei);
            ei = &ei->onext();
        }
    }

    return ei;
}


// used in regression test 06
// TODO: adapt test and/or delaunay::insert_edge to work with test and remove this function
Edge* insert_edge(MeshUnstructuredConstrained& mesh, Edge& edge, double hmin) {
    Edge* ei = mesh.get_initial_edge();

    bool has_edge = mesh.has_edge(&edge);
    if (!has_edge) {
        ei = insert_point(mesh, edge.org(), *ei, hmin);
        ei = refine_triangles(mesh, *ei, hmin);

        ei = insert_point(mesh, edge.dest(), *ei, hmin);
        ei = refine_triangles(mesh, *ei, hmin);

        assert(mesh.has_point(&edge.org()));
        assert(mesh.has_point(&edge.dest()));
    }

    // Do this regardless of whether the mesh already has the edge or not.
    // This will ensure the new edge (a source) will overwrite whatever else
    // is there (presumably NOT a source).
    delaunay::insert_edge(mesh, edge, *ei);
    // DO NOT refine triangles here, right after inserting the edge.
    // Do it elsewhere in a refinement algorithm.
    // Otherwise longer edge insertion will always get shortened.

    return ei;
}


// used by chew uniform, chew nonuniform, and ruppert
double calculate_minimum_length(const std::vector<Real2>& source_points,
        const std::list<parametric::IParametric*>& source_segments) {

    // Find the best hmin for the given sources
    double hmin = std::numeric_limits<double>::max();
    std::cout << "Initial hmin = " << hmin << '\n';

    math::optimize::ShortestDistanceParametricPoint dist_to_point;
    for (const Real2& p0 : source_points) {
        for (const Real2& p1 : source_points) {
            if (p0 == p1) continue;
            double dx = p0.first - p1.first;
            double dy = p0.second - p1.second;
            double d = sqrt(dx * dx + dy * dy);
            hmin = std::min(hmin, d);
        }
    }
    std::cout << "hmin after points to points = " << hmin << '\n';

    for (const Real2& p : source_points) {
        dist_to_point.set_point(p);
        for (parametric::IParametric* segment : source_segments) {
            dist_to_point.set_parametric(segment);
            double t = dist_to_point.run(&dist_to_point)[0];
            auto [x, y] = segment->evaluate(t);
            double dx = p.first - x;
            double dy = p.second - y;
            double d2 = dx * dx + dy * dy;
            if (inequalities::is_close(d2, 0)) continue;
            hmin = std::min(hmin, sqrt(d2));
        }
    }
    std::cout << "hmin after points to segments = " << hmin << '\n';

    math::optimize::ShortestDistanceParametricParametric dist_to_segments;
    for (parametric::IParametric* s0 : source_segments) {
        dist_to_segments.set_parametric_0(s0);
        auto [xi, yi] = s0->evaluate_tmin();
        auto [xf, yf] = s0->evaluate_tmax();
        double dxi = xf - xi;
        double dyi = yf - yi;
        double length = sqrt(dxi * dxi + dyi * dyi);
        hmin = std::min(hmin, length);
        for (parametric::IParametric* s1 : source_segments) {
            if (s0 == s1) continue;
            if (s0->evaluate_tmin() == s1->evaluate_tmin()) continue;
            if (s0->evaluate_tmin() == s1->evaluate_tmax()) continue;
            if (s0->evaluate_tmax() == s1->evaluate_tmin()) continue;
            if (s0->evaluate_tmax() == s1->evaluate_tmax()) continue;
            dist_to_segments.set_parametric_1(s1);
            std::vector<double> t = dist_to_segments.run(&dist_to_segments);
            auto [x0, y0] = s0->evaluate(t[0]);
            auto [x1, y1] = s1->evaluate(t[1]);
            double dx = x1 - x0;
            double dy = y1 - y0;
            double d = sqrt(dx * dx + dy * dy);
            hmin = std::min(hmin, d);
        }
    }
    std::cout << "hmin after segments = " << hmin << '\n';
    return hmin;
}


// used by chew uniform only
double split_sources_uniformly(std::list<parametric::IParametric*>& source_segments,
        double target_hmin) {

    math::integrate::SegmentLength integrator;
    static const double rtol = 1e-5;
    static const double atol = 0;
    auto comp = [](PairFloatParam& a, PairFloatParam& b) { return a.first < b.first; };
    std::priority_queue<PairFloatParam, std::vector<PairFloatParam>, decltype(comp)> pq(comp);
    double hmin = target_hmin;
    double hinit = target_hmin;

    for (auto segment : source_segments) {
        // get continuous length of parametric function
        double length = integrator.run(segment);
        // estimate number of splits
        double ns = length / hinit;
        double nsr = std::round(ns);
        bool is_close = inequalities::is_close(ns, nsr, rtol, atol);
        int ni = is_close ? nsr : std::ceil(ns);
        size_t n = std::max(ni, 1); // ensure never zero!
        // split segment with optimizer
        if (n > 1) {
            segment->set_num_subsegments(n);
            segment->optimize_parameters();
        }
        // check min/max segment lengths and confirm to ensure
        // subsegments are close in length
        const std::vector<double>& parameters = segment->get_parameters();
        double lmin = std::numeric_limits<double>::max();
        double lmax = std::numeric_limits<double>::min();
        for (size_t i = 0; i < parameters.size()-1; i++) {
            auto [xmin, ymin] = segment->evaluate(parameters[i]);
            auto [xmax, ymax] = segment->evaluate(parameters[i+1]);
            double dx = xmax - xmin;
            double dy = ymax - ymin;
            double len = sqrt(dx*dx + dy*dy);
            lmin = std::min(lmin, len);
            lmax = std::max(lmax, len);
        }
        bool bounds_close = inequalities::is_close(lmax, lmin, rtol, atol);
        if (!bounds_close) {
            std::cout << "WARNING: refine::split_sources_uniformly: A segment "
                "could not be split uniformly. A likely cause is using a nonlinear "
                "function to define a segment without initially splitting it "
                "fine enough to resolve all its features.\n"
                "Proceeding anyway.";
        }
        // queue segments based on (edge) length of subsegments
        pq.push({lmin, segment});
        // store queue min metrics
        hmin = std::min(hmin, lmin);
    }
    std::cout << "hmin after source (first pass) = " << hmin << '\n';

    static const double rtol_length = 1e-8;
    while (true) {
        auto [length, segment] = pq.top();
        double hmax = length;
        double hmax_allowed = hmin * sqrt(3);
        bool is_close = inequalities::is_close(hmax, hmax_allowed, rtol_length, atol);
        if (!is_close && hmax < hmax_allowed) break;
        pq.pop();
        while (true) {
            // refine further
            size_t nsub = segment->get_num_subsegments();
            segment->set_num_subsegments(nsub + 1);
            segment->optimize_parameters();

            // confirm subsegment lengths are close
            const std::vector<double>& parameters = segment->get_parameters();
            double lmin = std::numeric_limits<double>::max();
            double lmax = std::numeric_limits<double>::min();
            for (size_t i = 0; i < parameters.size()-1; i++) {
                auto [xmin, ymin] = segment->evaluate(parameters[i]);
                auto [xmax, ymax] = segment->evaluate(parameters[i+1]);
                double dx = xmax - xmin;
                double dy = ymax - ymin;
                double len = sqrt(dx*dx + dy*dy);
                lmin = std::min(lmin, len);
                lmax = std::max(lmax, len);
            }
            bool bounds_close = inequalities::is_close(lmax, lmin, rtol, atol);
            if (!bounds_close) {
                std::cerr << "ERROR: refine::split_sources_uniformly: segment was not split uniformly\n";
                std::cerr << "Try initially splitting the segment.\n";
                exit(-1);

            }

            is_close = inequalities::is_close(lmin, hmax_allowed, rtol_length, atol);
            if (!is_close && lmin < hmax_allowed) {
                hmin = std::min(hmin, lmin);
                pq.push({lmin, segment});
                break;
            }
        }

    }
    std::cout << "hmin after source (second pass) = " << hmin << '\n';

    // split segments
    auto segment_iter = source_segments.begin();
    while (segment_iter != source_segments.end()) {
        if ((*segment_iter)->get_num_subsegments() > 1) {
            (*segment_iter)->optimize_parameters();
            std::vector<std::shared_ptr<parametric::IParametric>> subsegments
                = (*segment_iter)->get_subsegments();
            parametric::IParametric* x = *segment_iter;
            ++segment_iter;
            delete x;
            source_segments.erase(std::prev(segment_iter));
            for (auto& subseg : subsegments) {
                if (subseg->get_num_subsegments() > 1)
                    subseg->optimize_parameters();
                source_segments.insert(segment_iter, subseg.get()->clone());
            }
        } else {
            ++segment_iter;
        }
    }

    return hmin;
}


// used by both ruppert and chew nonuniform
Edge& locate_or_closest_edge(const MeshUnstructuredConstrained& mesh, const Point& point, Edge& edge) {
    // nothing to do if the point is contained within the triangle
    if (math::left_of(point, edge) && math::left_of(point, edge.lnext())
            && math::left_of(point, edge.lnext().lnext()))
        return edge;

    // don't start on source
    Edge* e = &edge;
    if (mesh.is_source(e)) {
        e = &e->lnext();
        if (mesh.is_source(e)) {
            e = &e->lnext();
            if (mesh.is_source(e)) {
                // All 3 sides of the triangle are sources.
                // Point MUST be outside the triangle since
                // containment was already checked.
                // Edge closest to the point MUST be largest edge.
                Edge* e1 = &e->lnext();
                Edge* e2 = &e1->lnext();
                double l0 = e->org().distance_to(e->dest());
                double l1 = e1->org().distance_to(e1->dest());
                double l2 = e2->org().distance_to(e2->dest());
                if (l0 > l1 && l0 > l2)
                    return *e;
                if (l1 > l0 && l1 > l2)
                    return *e1;
                return *e2;
            }
        }
    }

    int iteration = 0;
    while (true) {
        if (iteration++ >= 10000) {
            std::cout << "delaunay::locate: Edge not found after 10000 iterations\n";
            exit(-1);
        }
        // TODO: can these point checks be included in the "else" break?
        //       must test left/right_of to see if endpoints get included
        if (e->org() == point || e->dest() == point) {
            break;
        } else if (mesh.is_source(e)) {
            if (!math::on_line(point, *e)) {
                if (math::left_of(point, *e))
                    e = &e->sym();
            }
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

    if (!quadedge::part_of_triangle(*e)) {
        Point& pa = e->org();
        Point& pb = e->dest();
        if (point == pa) {
            e = &e->onext();
        } else if (point == pb) {
            e = &e->sym();
        } else {
            assert(!math::right_of(point, *e)); // left of OR on line
        }
    } else {
        // If possible, make sure returned edge.org is the target point
        if (point == e->lnext().dest())
            e = &e->lnext().lnext();
        else if (point == e->dest())
            e = &e->lnext();
    }
    return *e;
}


// used by chew nonuniform and ruppert
bool convergence_test_angle(const Triangle& triangle, double min_angle) {
    return inequalities::is_gt(triangle.get_min_angle(), min_angle);
}


// used by chew nonuniform and ruppert
bool convergence_test_density(const Triangle& triangle, const density::DensityManager& density) {
    // Compares triangle radius to density.  This makes sense because
    // the triangle radius will be equal to the shortest distance from
    // the proposed point to a (remaining) point on the triangle.
    // If the density is greater than what would otherwise be the a new edge
    // length, than do not bother with the triangle.
    const Point& pc = triangle.get_circumcenter();
    return inequalities::is_le(triangle.get_radius(), density.evaluate(pc.x, pc.y));
}


algo::EdgeQueue fill_encroached_edge_queue(const MeshUnstructuredConstrained& mesh,
        Edge& edge, const Point& pc, double const triangle_radius) {

    auto is_encroached = [&](Edge* e) {
        if (mesh.is_source(e)) {
            double radius = e->org().distance_to(e->dest()) / 2;
            double x = (e->org().x + e->dest().x) / 2;
            double y = (e->org().y + e->dest().y) / 2;
            Point midpoint(x, y);
            double distance = pc.distance_to(midpoint);
            return inequalities::is_lt(distance, radius);
        }
        return false;
    };

    auto is_within_circumcircle = [&](Edge* e) {
        double distance = pc.distance_to(e->dest());
        return inequalities::is_le(distance, triangle_radius);
    };

    algo::EdgeQueue encroached_queue;

    if (is_encroached(&edge))
        encroached_queue.insert(edge);
    if (is_encroached(&edge.lnext()))
        encroached_queue.insert(edge.lnext());
    if (is_encroached(&edge.lnext().lnext()))
        encroached_queue.insert(edge.lnext().lnext());

    /*
      A major limitation of this is that segments blocked by other
      segments, e.g. a narrow cone, could still be added to the
      encroached edge queue.  On "quick" fix is to queue them anyway,
      and afterward filter out segments visible by testing with a
      locate function. Any segment that reaches the point are visible;
      any segment stuck at another segment is not visible. A more
      rigorous fix would be to check for intersections between the
      point, the segment midpoint, and every other segment in
      the encroached queue closer to the point.
    */
    std::set<Point*> visited;
    std::queue<Edge*> queue;
    queue.push(&edge);
    queue.push(&edge.lnext());
    queue.push(&edge.lnext().lnext());
    while (!queue.empty()) {
        Edge* e = queue.front();
        queue.pop();
        auto it = visited.insert(&e->dest());
        if (!it.second)
            continue;
        if (!is_within_circumcircle(e))
            continue;

        Edge* ei = &e->sym();
        if (!mesh.is_source_point(&ei->org())) {
            Edge* ef = ei;
            // Not a source point ==> no edges will be encroached
            do {
                if (!visited.contains(&ei->dest()))
                    queue.push(ei);
                ei = &ei->onext();
            } while (ei != ef);
        } else {
            // NB: source point DOES NOT mean source edges will exist.
            // It does mean that we can find up to 2 candidate segments.
            Edge* ef = &ei->oprev();
            while (!mesh.is_source(ei) && ei != ef) {
                if (!visited.contains(&ei->dest()))
                    queue.push(ei);
                ei = &ei->onext();
            }
            if (is_encroached(ei))
                encroached_queue.insert(*ei);
            // only try the other direction if CCW did not go all the way around
            if (ei != ef) {
                ei = &e->sym();
                ef = &ei->onext();
                if (mesh.is_source(ei))
                    ei = &ei->oprev();
                while (!mesh.is_source(ei) && ei != ef) {
                    if (!visited.contains(&ei->dest()))
                        queue.push(ei);
                    ei = &ei->oprev();
                }
                if (is_encroached(ei))
                    encroached_queue.insert(*ei);
            }
        }
    }

    return encroached_queue;
}


ClusterMetrics get_segment_cluster_metrics(MeshUnstructuredConstrained& mesh,
        Edge& encroached, double const min_angle_allowed) {

    assert(mesh.is_source(&encroached));

    // apex opposes side c
    auto get_angle = [](double a, double b, double c) {
        double cos_angle_rad = (a*a + b*b - c*c) / 2 / a / b;
        if (inequalities::is_close(cos_angle_rad, 1))
            cos_angle_rad = 1.0;
        if (inequalities::is_close(cos_angle_rad, -1))
            cos_angle_rad = -1.0;
        double angle_rad = std::acos(cos_angle_rad);
        assert(angle_rad > 0 && "Angle must be positive!\n");
        return angle_rad * 180 / std::numbers::pi;
    };

    // TODO: add checks to ensure angle only interior angles are measured
    auto get_segment_cluster = [&](Edge* e) -> std::pair<double, double> {
        std::deque<Edge*> cluster;
        Edge* prev = nullptr;
        Edge* curr = nullptr;
        double min_angle = 180;
        double min_length = std::numeric_limits<double>::max();

        // counterclockwise
        curr = &e->onext();
        prev = e;
        // stop after full loop or when outside the domain
        while (curr != e && curr->right().get_data()) {
            // get angle
            if (mesh.is_source(curr)) {
                double a = curr->org().distance_to(curr->dest());
                double b = prev->org().distance_to(prev->dest());
                double c = curr->dest().distance_to(prev->dest());
                double angle = get_angle(a, b, c);
                if (inequalities::is_gt(angle, 60))
                    break;
                min_angle = std::min(min_angle, angle);
                min_length = std::min(min_length, a);
                min_length = std::min(min_length, b);
                prev = curr;
            }
            curr = &curr->onext();
        }

        // clockwise
        curr = &e->oprev();
        prev = e;
        // stop after full loop or when outside the domain
        while (curr != e && curr->left().get_data()) {
            // get angle
            if (mesh.is_source(curr)) {
                double a = curr->org().distance_to(curr->dest());
                double b = prev->org().distance_to(prev->dest());
                double c = curr->dest().distance_to(prev->dest());
                double angle = get_angle(a, b, c);
                if (inequalities::is_gt(angle, 60))
                    break;
                min_angle = std::min(min_angle, angle);
                min_length = std::min(min_length, a);
                min_length = std::min(min_length, b);
                prev = curr;
            }
            curr = &curr->oprev();
        }

        return {min_length, min_angle};
    };

    auto [min_length, min_angle] = get_segment_cluster(&encroached);
    auto [min_length_sym, min_angle_sym] = get_segment_cluster(&encroached.sym());
    if (min_angle < min_angle_sym)
        return ClusterMetrics(min_length, min_angle, encroached.org(), &encroached);
    return ClusterMetrics(min_length_sym, min_angle_sym, encroached.dest(), &encroached.sym());
}


// used by chew nonuniform
Edge* split_edge(MeshUnstructuredConstrained& mesh, Edge& edge, int num_segments) {
    bool split_left = edge.left().get_data();
    bool split_right = edge.right().get_data();
    assert(split_left || split_right);

    std::vector<Edge*> new_edges = mesh.split_edge(edge, num_segments);

    // Segments should all be the same length. Get the min length
    // on the off chance they aren't (e.g. splitting some nonlinear multivalued function).
    double min_new_segment_length = new_edges[0]->org().distance_to(new_edges[0]->dest());
    for (auto it = new_edges.begin()+1; it != new_edges.end(); ++it)
        min_new_segment_length = std::min(min_new_segment_length, (*it)->org().distance_to((*it)->dest()));

    if (split_left)
        make_and_split_cavity(mesh, new_edges[1]->sym(), min_new_segment_length, new_edges[1]->org());
    if (split_right)
        make_and_split_cavity(mesh, *new_edges[1], min_new_segment_length, new_edges[1]->org());

    mesh.fill_triangles(*new_edges[0]);

    return new_edges[0];
}


std::vector<Edge*> split_encroached_edge_terminator(
        MeshUnstructuredConstrained& mesh,
        Edge& encroached, double const min_length_allowed,
        const ClusterMetrics& cluster, int num_segments) {

    assert(mesh.is_source(&encroached));

    auto get_isosceles_distance = [](double length, double angle) {
        double x = 2 * (1 - cos(angle * std::numbers::pi / 180));
        return length * sqrt(x);
    };

    auto is_power_of_2 = [](double length) -> std::pair<double, double> {
        assert(length > 0);
        int exponent;
        double significand = std::frexp(length, &exponent);
        bool already_power_of_2;
        double new_length = length;
        if (inequalities::is_close(significand, 0.5, 1e-4)) {
            already_power_of_2 = true;
        } else if (inequalities::is_close(significand, 1.0, 1e-4)) {
            already_power_of_2 = true;
        } else {
            already_power_of_2 = false;
            new_length = std::pow(2, exponent - 1);
        }
        return {new_length, already_power_of_2};
    };

    double encroached_edge_length = encroached.org().distance_to(encroached.dest());
    auto [encroached_power_of_2, encroached_already_power_of_2] = is_power_of_2(encroached_edge_length);
    auto [min_length_power_of_2, min_length_already_power_of_2] = is_power_of_2(cluster.length);

    double _min_length = get_isosceles_distance(cluster.length, cluster.angle);
    bool allow_terminator = inequalities::is_le(_min_length, 2 * min_length_allowed);
    if (!allow_terminator) {
        return mesh.split_edge(encroached, num_segments);
    }

    if (!min_length_already_power_of_2) {
        double length = is_power_of_2(cluster.length / 2).first;
        return segment_splitter_nonuniform(mesh, encroached, length, cluster);
    }

    if (!encroached_already_power_of_2) {
        if (inequalities::is_le(min_length_power_of_2, encroached_edge_length / 2)) {
            return segment_splitter_nonuniform(mesh, encroached, min_length_power_of_2, cluster);
        } else {
            return segment_splitter_nonuniform(mesh, encroached, min_length_power_of_2 / 2, cluster);
        }
    }

    if (inequalities::is_lt(1.5 * min_length_power_of_2, encroached_power_of_2))
        return segment_splitter_nonuniform(mesh, encroached, min_length_power_of_2, cluster);

    // Check min_length of isosceles triangle after splitting encroached length.
    // Comparing this to the min length allowed (i.e. density, previous triangle length)
    // makes the final decision whether to split or terminate.
    // ALTERNATIVELY could use half the encroached side length.
    // TODO: test the alternative and decide which is more robust
    double min_length = get_isosceles_distance(encroached_power_of_2 / 2, cluster.angle);
    if (inequalities::is_gt(min_length, min_length_allowed)) {
        assert(inequalities::is_close(encroached_power_of_2, encroached_edge_length));
        // if above always holds, can change this to mesh.split_edge
        // ALWAYS stick to 2 splits though, never 3 (an option for chew-nonuniform)
        // 2 splits is part of the convergence algorithm
        return segment_splitter_nonuniform(mesh, encroached, encroached_power_of_2 / 2, cluster);
    }

    return {};
}


std::vector<Edge*> segment_splitter_nonuniform(MeshUnstructuredConstrained& mesh,
        Edge& encroached, double length, const ClusterMetrics& cluster) {

    parametric::IParametric* parametric = mesh.get_segment_parametric(&encroached);
    if (parametric == nullptr)
        parametric = mesh.get_segment_parametric(&encroached.sym());
    assert(parametric);

    // Set parametric where to split at distance "length" from
    // the cluster's apex
    auto [x0, y0] = parametric->evaluate_tmin();
    if (cluster.apex.x == x0 && cluster.apex.y == y0) {
        parametric->point_from_start(length);
    } else {
        parametric->evaluate_tmax();
        parametric->point_from_end(length);
    }

    bool split_right = encroached.left().get_data() == nullptr;
    Edge& edge_to_split = split_right ? encroached.sym() : encroached;
    std::vector<Edge*> new_edges = mesh.split_edge(edge_to_split);

    if (split_right) {
        assert(new_edges.size() == 2);
        Edge* e = new_edges[0];
        new_edges[0] = &new_edges[1]->sym();
        new_edges[1] = &e->sym();
    }

    return new_edges;
}


/*
  This implementation deletes edges belonging to triangles that encircle the point.
  Source edges are skipped and will not get deleted.

  It has no sense on encroachement of points on the give edge, so no points will be deleted.
*/
void make_cavity(MeshUnstructuredConstrained& mesh, Edge& edge, const Point& point) {
    Edge* ei = &edge;
    Edge* efin = ei;
    do {
        // delete edge if right triangle encircles midpoint
        // since triangles are already (constrained) delaunay,
        // we only need to check if they encircle the new point
        if (!mesh.is_source(ei)) {
            Edge* er = &ei->sym();
            if (math::in_circle(er->org(), er->dest(), er->lnext().dest(), point)) {
                Edge* etmp = &ei->oprev();
                mesh.delete_edge(ei);
                ei = etmp;
            } else {
                ei = &ei->lnext();
            }
        } else {
            ei = &ei->lnext();
        }
    } while (ei != efin);
}


/*
  This is the "Rising Bubble" algorithm (see Guidbas et al.), adapted
  to work with constrained triangulations by stopping edge deletion
  at any source edge.

  No source points get deleted when splitting.
*/
void split_cavity(MeshUnstructuredConstrained& mesh, Edge& initial_edge) {
    Edge* basel = &initial_edge;

    while (true) {
        Edge* lcand = &basel->sym().onext();
        bool valid_lcand = math::right_of(lcand->dest(), *basel);
        if (valid_lcand) {
            while (math::in_circle(basel->dest(), basel->org(),
                            lcand->dest(), lcand->onext().dest())) {
                if (mesh.is_source(lcand))
                    break;
                mesh.delete_edge(lcand);
                lcand = &basel->sym().onext();
            }
        }

        Edge* rcand = &basel->oprev();
        bool valid_rcand = math::right_of(rcand->dest(), *basel);
        if (valid_rcand) {
            while (math::in_circle(basel->dest(), basel->org(),
                            rcand->dest(), rcand->oprev().dest())) {
                if (mesh.is_source(rcand))
                    break;
                mesh.delete_edge(rcand);
                rcand = &basel->oprev();
            }
        }

        if (lcand->dest() == rcand->dest())
            break;

        valid_lcand = math::right_of(lcand->dest(), *basel);
        valid_rcand = math::right_of(rcand->dest(), *basel);
        if (!valid_lcand && !valid_rcand) {
            // e.g. multiple segments on straight line cannot be joined
            basel = &basel->lprev();
        } else if (valid_lcand && !valid_rcand) {
            basel = &mesh.connect_edges(basel->sym(), lcand->sym());
        } else if (!valid_lcand && valid_rcand) {
            basel = &mesh.connect_edges(*rcand, basel->sym());
        } else {
            bool left_violates_delaunay = math::in_circle(lcand->dest(), lcand->org(),
                    rcand->org(), rcand->dest());
            if (left_violates_delaunay) {
                basel = &mesh.connect_edges(*rcand, basel->sym());
            } else {
                basel = &mesh.connect_edges(basel->sym(), lcand->sym());
            }
        }
    }
}


/*
  This is an implementation of the Rising Bubble algorithm (see Guibas et al.)
  that simultaneously makes and splits a cavity.

  Points within a distance of min_length from the provided point will get deleted.
  Care is taken to ensure no source points get deleted.

  Edges get deleted when their neighboring triangle encircles a point before or after
  in the mesh. Deleting edges stops at source edges to ensure they do not get deleted.

  The provided edge basel will not get deleted. Care must be take when selecting which
  edge to provide. It can take multiple runs for Delaunay criteria to be achieved.
  See insert_point.

  The returned basel will be different from the input basel if any edges get deleted.
*/
std::pair<Edge*, bool> make_and_split_cavity(MeshUnstructuredConstrained& mesh, Edge& basel_init,
        const double min_length, const Point& point) {

    Edge* basel = &basel_init;

    auto too_close = [&](Point& p) {
        if (mesh.is_source_point(&p))
            return false;
        if (point == p)
            return false;
        if (basel->org() == p)
            return false;
        if (basel->dest() == p)
            return false;
        if (basel_init.org() == p or basel_init.dest() == p)
            return false;
        double distance = point.distance_to(p);
        return inequalities::is_lt(distance, min_length);
    };

    bool any_points_deleted = false;
    auto remove_point = [&](Edge* e) {
        if (!too_close(e->dest()))
            return;
        if (mesh.is_source_point(&e->dest()))
            return;
        e = &e->sym();
        while (e != &e->onext())
            mesh.delete_edge(&e->onext());
        Point* p = &e->org();
        mesh.delete_edge(e);
        mesh.delete_point(p);
        any_points_deleted = true;
    };

    auto get_lcand = [&]() {
        Edge* lcand = &basel->sym().onext();
        while (too_close(lcand->dest())) {
            remove_point(lcand);
            lcand = &basel->sym().onext();
        }
        return lcand;
    };

    auto get_rcand = [&]() {
        Edge* rcand = &basel->oprev();
        while (too_close(rcand->dest())) {
            remove_point(rcand);
            rcand = &basel->oprev();
        }
        return rcand;
    };

    while (true) {
        Edge* lcand = get_lcand();
        bool valid_lcand = math::right_of(lcand->dest(), *basel);
        if (valid_lcand) {
            while (math::in_circle(basel->dest(), basel->org(), lcand->dest(), lcand->onext().dest())) {
                if (mesh.is_source(lcand) || lcand == &basel_init)
                    break;
                mesh.delete_edge(lcand);
                lcand = get_lcand();
            }
        }

        Edge* rcand = get_rcand();
        bool valid_rcand = math::right_of(rcand->dest(), *basel);
        if (valid_rcand) {
            while (math::in_circle(basel->dest(), basel->org(), rcand->dest(), rcand->oprev().dest())) {
                if (mesh.is_source(rcand) || rcand == &basel_init)
                    break;
                mesh.delete_edge(rcand);
                rcand = get_rcand();
            }
        }

        if (lcand->dest() == rcand->dest())
            break;

        valid_lcand = math::right_of(lcand->dest(), *basel);
        valid_rcand = math::right_of(rcand->dest(), *basel);
        if (!valid_lcand && !valid_rcand) {
            // This case occurs when angles between the basel and lcand/rcand
            // are both >= 180 degrees. Since this algorithm deletes edges and points,
            // just shifting basel is not ideal since the initial basel is
            // intended not to get removed.

            // Check if the next edges beyond lcand and rcand are < 180 degrees
            // with respect to those edges. If yes, new edges can be added.
            bool valid_left = math::right_of(lcand->sym().lprev().org(), *lcand);
            bool valid_right = math::left_of(rcand->lnext().dest(), *rcand);

            if (valid_left) {
                Edge* ei = &lcand->rprev();
                bool encircles_next = math::in_circle(ei->org(), ei->dest(), ei->lnext().dest(), lcand->org());
                if (encircles_next) {
                    bool next_source = mesh.is_source(ei);
                    if (!next_source)
                        mesh.delete_edge(ei);
                }
                mesh.connect_edges(lcand->sym(), lcand->sym().lprev());
            }

            if (valid_right) {
                Edge* ei = &rcand->lnext().sym();
                bool encircles_next = math::in_circle(ei->org(), ei->dest(),
                        ei->lnext().dest(), rcand->org());
                if (encircles_next) {
                    bool next_source = mesh.is_source(ei);
                    if (!next_source)
                        mesh.delete_edge(ei);
                }
                mesh.connect_edges(rcand->lnext(), *rcand);
            }

            // The only possible way neither are valid is if
            // all points tested are on a line. The only way
            // to proceed is by shifting basel.
            if (!valid_left && !valid_right)
                basel = &basel->rprev();

        } else if (valid_lcand && !valid_rcand) {
            basel = &mesh.connect_edges(basel->sym(), lcand->sym());
        } else if (!valid_lcand && valid_rcand) {
            basel = &mesh.connect_edges(*rcand, basel->sym());
        } else {
            bool left_violates_delaunay = math::in_circle(lcand->dest(), lcand->org(),
                    rcand->org(), rcand->dest());
            if (left_violates_delaunay) {
                basel = &mesh.connect_edges(*rcand, basel->sym());
            } else {
                basel = &mesh.connect_edges(basel->sym(), lcand->sym());
            }
        }

        // The triangle to the left of the new basel is not guaranteed
        // to satify the Delaunay criteria! Must flip it as necessary!
        // NB: basel can never be basel_init; candidates can be but they
        // have no risk of getting deleted here
        while (!mesh.is_source(basel)) {
            lcand = get_lcand();
            rcand = get_rcand();
            bool left_violates_delaunay = math::in_circle(basel->org(), basel->dest(),
                    basel->lnext().dest(), lcand->dest());
            bool right_violates_delaunay = math::in_circle(basel->org(), basel->dest(),
                    basel->lnext().dest(), rcand->dest());
            if (left_violates_delaunay) {
                mesh.delete_edge(basel);
                basel = &mesh.connect_edges(lcand->sym().lnext(), lcand->sym());
            } else if (right_violates_delaunay) {
                mesh.delete_edge(basel);
                basel = &mesh.connect_edges(*rcand, rcand->lprev());
            } else {
                break;
            }
        }
    }

    return {basel, any_points_deleted};
}


} // namespace refine

} // namespace mesh

} // namespace umr
