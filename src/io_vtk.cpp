#include "io_vtk.hpp"

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


void VTK::set_cells(const MeshUnstructured& mesh) {
    auto& points = mesh.get_points_set(); // const unordered set of point pointers
    m_point_indices.clear();
    int idx = 0;
    vtkNew<vtkPoints> grid_points;
    for (auto point : points) {
        grid_points->InsertNextPoint(point->x, point->y, 0.0);
        m_point_indices[point] = idx++;
    }
    m_grid->SetPoints(grid_points);

    auto& edges = mesh.get_edges_set(); // const unordered set of edges pointers
    std::unordered_map<Triangle*, bool> already_added;
    for (auto edge : edges) {
        Triangle* tri = edge->left().get_data();
        if (tri && !already_added[tri]) {
            already_added[tri] = true;
            vtkNew<vtkGenericCell> cell;
            cell->SetCellTypeToTriangle();
            int i = 0;
            for (const Point* p : tri->get_points()) {
                int idx = m_point_indices[p];
                cell->GetPointIds()->SetId(i++, idx);
            }
            m_grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
        }

        tri = edge->right().get_data();
        if (tri && !already_added[tri]) {
            already_added[tri] = true;
            vtkNew<vtkGenericCell> cell;
            cell->SetCellTypeToTriangle();
            int i = 0;
            for (const Point* p : tri->get_points()) {
                int idx = m_point_indices[p];
                cell->GetPointIds()->SetId(i++, idx);
            }
            m_grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
        }
    }
}

void VTK::set_sources(const MeshUnstructuredConstrained& mesh) {
    vtkNew<vtkIntArray> source_points;
    source_points->SetName("SourcePoints");
    source_points->SetNumberOfComponents(1);
    for (auto p : mesh.get_source_points())
        source_points->InsertNextValue(m_point_indices[p]);

    vtkNew<vtkIntArray> source_segments;
    source_segments->SetName("SourceSegments");
    source_segments->SetNumberOfComponents(2);
    for (auto segment_param : mesh.get_source_edges()) {
        Edge* e = segment_param.first;
        source_segments->InsertNextValue(m_point_indices[&e->org()]);
        source_segments->InsertNextValue(m_point_indices[&e->dest()]);
    }

    vtkNew<vtkFieldData> fields;
    fields->AddArray(source_points);
    fields->AddArray(source_segments);
    m_grid->SetFieldData(fields);
}


void VTK::write(std::string& filename) {
    std::cout << "Writing mesh to filename " << filename << '\n';

    vtkNew<vtkUnstructuredGridWriter> writer;
    writer->SetFileName(filename.c_str());
    writer->SetInputData(m_grid);
    writer->Write();
}


} // namespace io

} // namespace mesh

} // namespace umr
