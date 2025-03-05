#include "sources_shape_triangle.hpp"

#include <cassert>

#include "math.hpp"
#include "parametrics.hpp"


namespace umr {

namespace source {


ShapeTriangle::ShapeTriangle(mesh::Point& p0, mesh::Point& p1, mesh::Point& p2) {
    validate_construction(p0, p1, p2);
    add_segment(std::make_shared<parametric::Line>(p0, p1));
    add_segment(std::make_shared<parametric::Line>(p1, p2));
    add_segment(std::make_shared<parametric::Line>(p2, p0));
}

void ShapeTriangle::validate_construction(
        mesh::Point& p0, mesh::Point& p1, mesh::Point& p2) {
    unsigned int duplicates = p0 == p1;
    duplicates += p0 == p2;
    duplicates += p1 == p2;
    assert(duplicates == 0 &&
            "ShapeTriangle::validate_construction: 3 points must be unique!");

    if (math::on_line(p0, p1, p2)) {
        std::cout << "ShapeTriangle::validate_construction: "
            "3 points must NOT be on the same line\n";
        exit(-1);
    }
}

void ShapeTriangle::set_number_subsegments(std::initializer_list<size_t> nseg) {
    assert(nseg.size() == 3
            && "ShapeTriangle::seg_number_subsegments: requires 4 arguments");
    IShape::set_number_subsegments(nseg);
}


} // namespace source

} // namespace umr
