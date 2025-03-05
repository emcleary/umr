#include "parametrics_general.hpp"

#include "integrators.hpp"
#include "inequalities.hpp"


namespace umr {

namespace parametric {


Parametric::Parametric(Parametric& parametric)
        : IParametric(parametric),
          m_function(parametric.m_function),
          m_derivative(parametric.m_derivative) {}


Parametric::Parametric(Parametric& parametric, double tmin, double tmax)
        : IParametric(parametric, tmin, tmax),
          m_function(parametric.m_function), m_derivative(parametric.m_derivative) {
    initialize(parametric);
}


Parametric::Parametric(ParametricFunction function, const double tmin, const double tmax)
        : IParametric(tmin, tmax), m_function(function) {
    initialize();

    bool is_periodic = inequalities::is_close(m_x0_b, m_x1_b)
        && inequalities::is_close(m_y0_b, m_y1_b);
    if (is_periodic)
        set_num_subsegments(3);
}


Parametric::Parametric(ParametricFunction function, ParametricFunction derivative,
        const double tmin, const double tmax)
        : IParametric(tmin, tmax),
          m_function(function), m_derivative(derivative) {
    initialize();

    bool is_periodic = inequalities::is_close(m_x0_b, m_x1_b)
        && inequalities::is_close(m_y0_b, m_y1_b);
    if (is_periodic)
        set_num_subsegments(3);
}


Real2 Parametric::evaluate_model(const double t) const {
    return m_function(t);
}


Real2 Parametric::differentiate_model(double t) const {
    if (m_derivative)
        return m_derivative(t);
    return IParametric::differentiate_model(t);
}


double Parametric::get_length() {
    math::integrate::SegmentLength integrator;
    return integrator.run(this);
}


} // namespace parametric

} // namespace umr
