#ifndef IO_VTK_H
#define IO_VTK_H

#include <json/json.h>
#include <string>
#include <vtkUnstructuredGrid.h>

#include "io.hpp"
#include "mesh.hpp"


namespace umr {

namespace mesh {

namespace io {


/**
 * @class VTK
 * @brief Dumps mesh data to VTK files
 *
 * This class dumps mesh data to VTK files in ASCII format.  Data
 * includes all points and edges in the mesh, along with source points
 * and segments for constrained meshes.
 */
class VTK : public IO {
public:
    /** Constructor for VTK */
    VTK() {}

    /** Destructor for VTK */
    virtual ~VTK() {};

    /** Sets grid data using the given mesh */
    virtual void set_cells(const MeshUnstructured& mesh) final;
    virtual void set_sources(const MeshUnstructuredConstrained& mesh) final;

    /** Writes grid data to file */
    virtual void write(std::string& filename) final;

protected:
    vtkNew<vtkUnstructuredGrid> m_grid;
    std::unordered_map<const Point*, int> m_point_indices;
};


} // namespace io

} // namespace mesh

} // namespace umr

#endif // IO_VTK_H
