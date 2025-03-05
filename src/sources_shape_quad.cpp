#include "sources_shape_quad.hpp"

#include <cassert>

#include "math.hpp"
#include "parametrics.hpp"


namespace umr {

namespace source {


ShapeQuadrilateral::ShapeQuadrilateral(mesh::Point& p0, mesh::Point& p1,
        mesh::Point& p2, mesh::Point& p3) {
    validate_construction(p0, p1, p2, p3);
    add_segment(std::make_shared<parametric::Line>(p0, p1));
    add_segment(std::make_shared<parametric::Line>(p1, p2));
    add_segment(std::make_shared<parametric::Line>(p2, p3));
    add_segment(std::make_shared<parametric::Line>(p3, p0));
    assert(is_closed() && "ShapeQuadrilateral::ShapeQuadrilateral: 4 points didn't close");
}

void ShapeQuadrilateral::set_number_subsegments(std::initializer_list<size_t> nseg) {
    assert(nseg.size() == 4 &&
            "ShapeQuadrilateral::seg_number_subsegments: requires 4 arguments");
    IShape::set_number_subsegments(nseg);
}

void ShapeQuadrilateral::validate_construction(mesh::Point& p0, mesh::Point& p1,
        mesh::Point& p2, mesh::Point& p3) {
    unsigned int duplicates = p0 == p1;
    duplicates += p0 == p2;
    duplicates += p0 == p3;
    duplicates += p1 == p2;
    duplicates += p1 == p3;
    duplicates += p2 == p3;
    assert(duplicates == 0 &&
            "ShapeQuadrilateral::validate_construction: 4 points must be unique!");

    int num_ccw = math::is_ccw(p0, p1, p2);
    num_ccw += math::is_ccw(p1, p2, p3);
    num_ccw += math::is_ccw(p2, p3, p0);
    num_ccw += math::is_ccw(p3, p0, p1);
    assert(num_ccw != 2 &&
            "ShapeQuadrilateral::validate_construction: Line segments cannot cross!");
}


} // namespace source

} // namespace umr
