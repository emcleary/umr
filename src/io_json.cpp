#include "io_json.hpp"

#include <iostream>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>

#include "inequalities.hpp"
#include "point.hpp"
#include "utilities.hpp"


namespace umr {

namespace mesh {

namespace io {


void JSON::set_cells(const MeshUnstructured& mesh) {
    auto& points_set = mesh.get_points_set();
    std::vector<Point*> points_vec(points_set.begin(), points_set.end());
    std::sort(points_vec.begin(), points_vec.end(), [](Point* pa, Point* pb) {
        if (pa->x == pb->x)
            return inequalities::is_lt(pa->y, pb->y);
        return inequalities::is_lt(pa->x, pb->x);
    });
    for (int i = 0; i < (int)points_vec.size(); ++i) {
        m_output["points"][i][0] = points_vec[i]->x;
        m_output["points"][i][1] = points_vec[i]->y;
        m_point_indices[points_vec[i]] = i;
    }

    Json::Value::ArrayIndex e_idx = 0;
    for (auto edge : mesh.get_edges_set()) {
        auto it_org = m_point_indices.find(&edge->org());
        auto it_dest = m_point_indices.find(&edge->dest());
        assert(it_org != m_point_indices.end() &&
                "umr::mesh::io::JSON: edge org not in points!");
        assert(it_dest != m_point_indices.end() &&
                "umr::mesh::io::JSON: edge dest not in points!");
        m_output["edges"][e_idx][0] = (*it_org).second;
        m_output["edges"][e_idx][1] = (*it_dest).second;
        ++e_idx;
    }
}


void JSON::set_sources(const MeshUnstructuredConstrained& mesh) {
    Json::Value::ArrayIndex pt_idx = 0;
    for (auto point : mesh.get_source_points()) {
        auto it = m_point_indices.find(point);
        assert(it != m_point_indices.end() &&
                "umr::mesh::io::JSON::set_sources: source point not in points!");
        m_output["source_points"][pt_idx] = (*it).second;
        ++pt_idx;
    }

    Json::Value::ArrayIndex e_idx = 0;
    for (auto [edge, seg] : mesh.get_source_edges()) {
        auto it_org = m_point_indices.find(&edge->org());
        auto it_dest = m_point_indices.find(&edge->dest());
        assert(it_org != m_point_indices.end() &&
                "umr::mesh::io::JSON::set_sources: edge org not in points!");
        assert(it_dest != m_point_indices.end() &&
                "umr::mesh::io::JSON::set_sources: edge dest not in points!");
        m_output["source_segments"][e_idx][0] = (*it_org).second;
        m_output["source_segments"][e_idx][1] = (*it_dest).second;
        ++e_idx;
    }
}


void JSON::write(std::string& filename) {
    std::ofstream ofs(filename);
    ofs << m_output;
    ofs.close();
}


} // namespace io

} // namespace mesh

} // namespace umr
