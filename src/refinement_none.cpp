#include "refinement_none.hpp"

#include "delaunay.hpp"
#include "mesh.hpp"


namespace umr {

namespace mesh {

namespace algo {


void Refinement::initialize(MeshUnstructuredConstrained& mesh,
        std::vector<Real2>& points,
        std::list<parametric::IParametric*>& segments) {

    for (auto [x, y] : points)
        mesh.make_point(x, y);
}


bool Refinement::iteration(MeshUnstructuredConstrained& mesh) {
    delaunay::divide_and_conquer(mesh);
    mesh.fill_triangles();
    return true;
}


} // namespace algo

} // namespace mesh

} // namespace umr
