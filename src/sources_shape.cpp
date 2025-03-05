#include "sources_shape.hpp"

#include "inequalities.hpp"


namespace umr {

namespace source {


void Shape::add_line(mesh::Point& p0, mesh::Point& p1, size_t nseg) {
    auto line = std::make_shared<parametric::Line>(p0, p1);
    line->set_num_subsegments(nseg);
    add_segment(line);
}

void Shape::add_arc(mesh::Point& pc, double radius,
        double tmin, double tmax, size_t nseg) {
    auto arc = std::make_shared<parametric::Arc>(pc, radius, tmin, tmax);
    arc->set_num_subsegments(nseg);
    add_segment(arc);
}

void Shape::add_parametric(parametric::ParametricFunction function,
        double tmin, double tmax, size_t nseg) {
    auto parametric = std::make_shared<parametric::Parametric>(function, tmin, tmax);
    parametric->set_num_subsegments(nseg);
    add_segment(parametric);
}


} // namespace source

} // namespace umr
