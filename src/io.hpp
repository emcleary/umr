#ifndef IO_H
#define IO_H

#include <json/json.h>
#include <string>
#include <vtkUnstructuredGrid.h>

#include "mesh.hpp"


namespace umr {

namespace mesh {

namespace io {


class IO {
public:
    /** Constructor for IO */
    IO() {}

    /** Destructor for IO */
    virtual ~IO() {};

    /** Sets grid data using the given mesh */
    virtual void set_cells(const MeshUnstructured& mesh) = 0;
    virtual void set_sources(const MeshUnstructuredConstrained& mesh) = 0;

    /** Writes grid data to file */
    virtual void write(std::string& filename) = 0;

protected:
};


/**
 * Dumps mesh to a file. Accepted file extensions are VTK and JSON.
 *
 * @param mesh A mesh
 * @param filename Output filename
 * @param prepend Prepends filename with indexing, default true
 */
void dump(const MeshUnstructured& mesh, std::string&& filename, bool prepend = true);

/**
 * Dumps mesh to a file. Accepted file extensions are VTK and JSON.
 *
 * @param mesh A mesh
 * @param filename Output filename
 * @param prepend Prepends filename with indexing, default true
 */
void dump(const MeshUnstructured& mesh, std::string& filename, bool prepend = true);


} // namespace io

} // namespace mesh

} // namespace umr

#endif // IO_H
