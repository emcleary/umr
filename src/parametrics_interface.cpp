#include <cassert>

#include "parametrics_interface.hpp"
#include "optimizers.hpp"
#include "inequalities.hpp"


namespace umr {

namespace parametric {


IParametric::IParametric(IParametric& ip) : m_x0_b(ip.m_x0_b), m_y0_b(ip.m_y0_b),
					    m_x1_b(ip.m_x1_b), m_y1_b(ip.m_y1_b),
					    m_num_subsegments(ip.m_num_subsegments) {
    m_parameters.insert(m_parameters.begin(), ip.m_parameters.begin(), ip.m_parameters.end());
}


IParametric::IParametric(IParametric& ip, const double tmin, const double tmax)
        : m_x0_b(ip.m_x0_b), m_y0_b(ip.m_y0_b),
          m_x1_b(ip.m_x1_b), m_y1_b(ip.m_y1_b) {
    assert(inequalities::is_ge(tmin, ip.get_tmin())
            && "IParametric::IParametric: tmin must be >= parametric's tmin");
    assert(inequalities::is_le(tmax, ip.get_tmax())
            && "IParametric::IParametric: tmax must be <= parametric's tmax");
    m_parameters.push_back(tmin);
    m_parameters.push_back(tmax);
}


IParametric::IParametric() {
    m_parameters.push_back(0);
    m_parameters.push_back(1);
}


IParametric::IParametric(const double tmin, const double tmax) {
    m_parameters.push_back(tmin);
    m_parameters.push_back(tmax);
}


Real2 IParametric::evaluate(const double t) const {
    static const double tol = 1e-15;
    assert(validate_parameter_in_bounds(t));
    if (inequalities::is_le(t, get_tmin(), tol, 0)) return evaluate_tmin();
    if (inequalities::is_ge(t, get_tmax(), tol, 0)) return evaluate_tmax();
    return evaluate_model(t);
}


Real2 IParametric::differentiate_model(double t) const {
    assert(validate_parameter_in_bounds(t));
    static const double dt = 1e-8;
    auto [xp, yp] = evaluate(t + dt);
    auto [xm, ym] = evaluate(t - dt);
    double dxdt = (xp - xm) / dt / 2;
    double dydt = (yp - ym) / dt / 2;
    return {dxdt, dydt};
}


void IParametric::set_num_subsegments(size_t n) {
    double tmax = m_parameters[m_num_subsegments];
    m_num_subsegments = n;
    m_parameters.resize(m_num_subsegments + 1);
    m_parameters[m_num_subsegments] = tmax;
}


std::vector<std::shared_ptr<IParametric>> IParametric::get_subsegments() {
    std::vector<std::shared_ptr<IParametric>> subsegments;
    for (size_t i = 0; i < m_num_subsegments; ++i) {
        std::shared_ptr<IParametric> subsegment;
        subsegment.reset(split(m_parameters[i], m_parameters[i+1]));
        subsegments.push_back(subsegment);
    }
    return subsegments;
}


void IParametric::set_parameters(std::vector<double>& params) {
    m_parameters.clear();
    m_parameters.insert(m_parameters.begin(), params.begin(), params.end());
}


void IParametric::set_parameters(std::vector<double>&& params) {
    m_parameters.clear();
    m_parameters.insert(m_parameters.begin(), params.begin(), params.end());
}


void IParametric::set_parameter_bounds(double tmin, double tmax) {
    m_parameters[0] = tmin;
    m_parameters[m_num_subsegments] = tmax;
}


void IParametric::set_initial_point(double x, double y) {
    assert(inequalities::is_close(x, m_x0_b));
    assert(inequalities::is_close(y, m_y0_b));
    m_x0_b = x;
    m_y0_b = y;
}


void IParametric::set_final_point(double x, double y) {
    assert(inequalities::is_close(x, m_x1_b));
    assert(inequalities::is_close(y, m_y1_b));
    m_x1_b = x;
    m_y1_b = y;
}


void IParametric::optimize_parameters() {
    math::optimize::SplitSegmentUniformly opt;
    opt.set_num_subsegments(m_num_subsegments);
    std::vector<double> parameters = opt.run(this);
    m_parameters.clear();
    m_parameters.insert(m_parameters.begin(), parameters.begin(), parameters.end());
}


Real2 IParametric::point_from_start(double length) {
    math::optimize::SplitSegmentAtLength opt;
    opt.set_initial_parameter(get_tmin());
    opt.set_point(evaluate_tmin());
    opt.set_parametric(this);
    opt.set_length(length);
    double t = opt.run(&opt)[0];
    double tmin = get_tmin();
    double tmax = get_tmax();
    m_parameters.resize(3);
    m_parameters[0] = tmin;
    m_parameters[1] = t;
    m_parameters[2] = tmax;
    return evaluate_model(t);
}


Real2 IParametric::point_from_end(double length) {
    math::optimize::SplitSegmentAtLength opt;
    opt.set_initial_parameter(get_tmax());
    opt.set_point(evaluate_tmax());
    opt.set_parametric(this);
    opt.set_length(length);
    double t = opt.run(&opt)[0];
    double tmin = get_tmin();
    double tmax = get_tmax();
    m_parameters.resize(3);
    m_parameters[0] = tmin;
    m_parameters[1] = t;
    m_parameters[2] = tmax;
    return evaluate_model(t);
}


bool IParametric::validate_parameter_in_bounds(double t) const {
    bool valid = inequalities::is_ge(t, get_tmin()) && inequalities::is_le(t, get_tmax());
    if (!valid) {
        std::cout << "Invlaid parameter t = " << t << '\n';
        std::cout << "Parameter bounds [" << get_tmin() << ", " << get_tmax() << "]\n";
    }
    return valid;
}


void IParametric::initialize() {
    auto [xmin, ymin] = evaluate_model(get_tmin());
    m_x0_b = xmin;
    m_y0_b = ymin;

    auto [xmax, ymax] = evaluate_model(get_tmax());
    // ensure no precision error at bounds if periodic
    bool is_periodic = inequalities::is_close(xmin, xmax) && inequalities::is_close(ymin, ymax);
    m_x1_b = is_periodic ? m_x0_b : xmax;
    m_y1_b = is_periodic ? m_y0_b : ymax;
}


void IParametric::initialize(IParametric& ip) {
    initialize();

    // ensure no precision error at bounds of the parametric
    if (inequalities::is_close(get_tmin(), ip.get_tmin())) {
        m_x0_b = ip.m_x0_b;
        m_y0_b = ip.m_y0_b;
    }

    if (inequalities::is_close(get_tmax(), ip.get_tmax())) {
        m_x1_b = ip.m_x1_b;
        m_y1_b = ip.m_y1_b;
    }
}


} // namespace parametric

} // namespace umr
