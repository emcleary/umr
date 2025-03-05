#include "sources_shape_circle.hpp"

#include "parametrics.hpp"


namespace umr {

namespace source {


ShapeCircle::ShapeCircle(mesh::Point& pc, double radius) {
    add_segment(std::make_shared<parametric::Arc>(pc, radius, 0, 2 * std::numbers::pi));
}


void ShapeCircle::set_number_subsegments(std::initializer_list<size_t> nseg) {
    assert(nseg.size() == 1 &&
            "ShapeCircle::seg_number_subsegments: requires exactly 1 argument");
    assert(*nseg.begin() >= 3 &&
            "ShapeCircle::seg_number_subsegments: shape requires at least 3 subsegments");
    IShape::set_number_subsegments(nseg);
}


} // namespace source

} // namespace umr
