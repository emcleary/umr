#include "optimizers.hpp"

#include "nonlinear_least_squares.hpp"


namespace umr {

namespace math {

namespace optimize {


void segment_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr) {

    parametric::IParametric* parametric = static_cast<parametric::IParametric*>(ptr);

    double avg = 0;
    auto [xprev, yprev] = parametric->evaluate_tmin();
    for (alglib::ae_int_t i = 0; i < t.length(); i++) {
        auto [x, y] = parametric->evaluate(t[i]);
        double dx = x - xprev;
        double dy = y - yprev;
        fi[i] = sqrt(dx * dx + dy * dy);
        xprev = x;
        yprev = y;
        avg += fi[i];
    }
    auto [x, y] = parametric->evaluate_tmax();
    double dx = x - xprev;
    double dy = y - yprev;
    fi[fi.length()-1] = sqrt(dx * dx + dy * dy);
    avg += fi[fi.length()-1];

    avg /= fi.length();
    for (alglib::ae_int_t i = 0; i < fi.length(); i++)
        fi[i] -= avg;
}


void shortest_distance_parametric_parametric_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr) {
    ShortestDistanceParametricParametric* models =
        static_cast<ShortestDistanceParametricParametric*>(ptr);
    const parametric::IParametric* const param0 = models->get_parametric_0();
    const parametric::IParametric* const param1 = models->get_parametric_1();
    auto [x0, y0] = param0->evaluate(t[0]);
    auto [x1, y1] = param1->evaluate(t[1]);
    double dx = x1 - x0;
    double dy = y1 - y0;
    fi[0] = dx;
    fi[1] = dy;
}


void shortest_distance_parametric_point_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr) {
    ShortestDistanceParametricPoint* models =
        static_cast<ShortestDistanceParametricPoint*>(ptr);
    const parametric::IParametric* const param = models->get_parametric();
    const Real2& point = models->get_point();
    auto [x0, y0] = param->evaluate(t[0]);
    double dx = point.first - x0;
    double dy = point.second - y0;
    fi[0] = dx;
    fi[1] = dy;
}


void split_segment_at_length_residual(const alglib::real_1d_array &t,
        alglib::real_1d_array &fi, void *ptr) {
    SplitSegmentAtLength* models = static_cast<SplitSegmentAtLength*>(ptr);
    const parametric::IParametric* const param = models->get_parametric();
    const Real2& point = models->get_point();
    const double length = models->get_length();
    auto [x0, y0] = param->evaluate(t[0]);
    double dx = point.first - x0;
    double dy = point.second - y0;
    double d2 = dx * dx + dy * dy;
    fi[0] = length * length - d2;
}


void SplitSegmentUniformly::initialize(void* ptr) {
    parametric::IParametric* parametric = static_cast<parametric::IParametric*>(ptr);
    assert(m_num_resid == m_num_params + 1
            && "SplitSegmentUniformly::initialize: params and resid sizes not set correctly!");

    const double tmin = parametric->get_tmin();
    const double tmax = parametric->get_tmax();
    for (int i = 0; i < m_num_params; i++) {
        m_params[i] = tmin + (tmax - tmin) * (i+1) / (m_num_params + 1);
        m_scale[i] = 1;
        m_lower_bound[i] = tmin;
        m_upper_bound[i] = tmax;
    }
}

std::vector<double> SplitSegmentUniformly::finalize(void* ptr) {
    parametric::IParametric* parametric = static_cast<parametric::IParametric*>(ptr);
    std::vector<double> sol(m_num_params + 2);
    sol[0] = parametric->get_tmin();
    for (int i = 0; i < m_num_params; i++)
        sol[i+1] = m_params[i];
    sol[m_num_params+1] = parametric->get_tmax();
    return sol;
}


void SplitSegmentUniformly::set_num_subsegments(size_t n_subsegments) {
    assert(n_subsegments > 0 && "SplitSegmentUniformly::set_num_subsegments: "
            "must split into at least 1 subsegment!");
    set_num_params(n_subsegments-1);
    set_num_resid(n_subsegments);
}



void ShortestDistanceParametricParametric::initialize(void* ptr) {
    ShortestDistanceParametricParametric* models =
        static_cast<ShortestDistanceParametricParametric*>(ptr);
    const parametric::IParametric* const param0 = models->get_parametric_0();
    const parametric::IParametric* const param1 = models->get_parametric_1();
    const double t0_min = param0->get_tmin();
    const double t0_max = param0->get_tmax();
    const double t1_min = param1->get_tmin();
    const double t1_max = param1->get_tmax();
    m_params[0] = (t0_min + t0_max) / 2;
    m_params[1] = (t1_min + t1_max) / 2;
    m_lower_bound[0] = t0_min;
    m_upper_bound[0] = t0_max;
    m_lower_bound[1] = t1_min;
    m_upper_bound[1] = t1_max;
    m_scale[0] = 1;
    m_scale[1] = 1;
}


std::vector<double> ShortestDistanceParametricParametric::finalize(void* ptr) {
    std::vector<double> sol(2);
    sol[0] = m_params[0];
    sol[1] = m_params[1];
    return sol;
}


void ShortestDistanceParametricPoint::initialize(void* ptr) {
    ShortestDistanceParametricPoint* models = static_cast<ShortestDistanceParametricPoint*>(ptr);
    const parametric::IParametric* const param = models->get_parametric();
    const double t0_min = param->get_tmin();
    const double t0_max = param->get_tmax();
    m_params[0] = (t0_min + t0_max) / 2;
    m_lower_bound[0] = t0_min;
    m_upper_bound[0] = t0_max;
    m_scale[0] = 1;
}


std::vector<double> ShortestDistanceParametricPoint::finalize(void* ptr) {
    std::vector<double> sol(1);
    sol[0] = m_params[0];
    return sol;
}


void SplitSegmentAtLength::initialize(void* ptr) {
    SplitSegmentAtLength* models = static_cast<SplitSegmentAtLength*>(ptr);
    const parametric::IParametric* const param = models->get_parametric();
    m_params[0] = models->get_initial_parameter();
    m_lower_bound[0] = param->get_tmin();
    m_upper_bound[0] = param->get_tmax();
    m_scale[0] = 1;
}


std::vector<double> SplitSegmentAtLength::finalize(void* ptr) {
    std::vector<double> sol(1);
    sol[0] = m_params[0];
    return sol;
}


} // namespace optimize

} // namespace math

} // namespace umr
