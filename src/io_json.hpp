#ifndef IO_JSON_H
#define IO_JSON_H

#include <json/json.h>
#include <string>
#include <vtkUnstructuredGrid.h>

#include "io.hpp"
#include "mesh.hpp"


namespace umr {

namespace mesh {

namespace io {


/**
 * @class JSON
 * @brief Dumps mesh data to JSON files
 *
 * This class dumps mesh data to a JSON file.  Data includes all
 * points and edges in the mesh, along with source points and segments
 * for constrained meshes.
 */
class JSON : public IO {
public:
    /** Constructor for JSON */
    JSON() {}

    /** Destructor for JSON */
    virtual ~JSON() {};

    /** Sets grid data using the given mesh */
    virtual void set_cells(const MeshUnstructured& mesh) final;
    virtual void set_sources(const MeshUnstructuredConstrained& mesh) final;

    /** Writes grid data to file */
    virtual void write(std::string& filename) final;

protected:
    Json::Value m_output;
    std::unordered_map<const Point*, int> m_point_indices;
};


} // namespace io

} // namespace mesh

} // namespace umr

#endif // IO_JSON_H
