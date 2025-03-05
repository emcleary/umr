#include "io.hpp"

#include <iostream>

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>

#include "io_json.hpp"
#include "io_vtk.hpp"
#include "inequalities.hpp"
#include "point.hpp"
#include "utilities.hpp"


namespace umr {

namespace mesh {

namespace io {


static int dump_counter = 0;
void prepend_filename(std::string& filename) {
    std::string c = std::to_string(dump_counter++);
    size_t target = 3;
    while (c.length() < target) {
        c = '0' + c;
    }
    filename = "data_" + c + "_" + filename;
}


void dump(const MeshUnstructured& mesh, std::string&& filename, bool prepend) {
    dump(mesh, filename, prepend);
}


void dump(const MeshUnstructured& mesh, std::string& filename, bool prepend) {
    std::unique_ptr<IO> io;
    if (is_vtk_file(filename))
        io = std::make_unique<VTK>();
    else if (is_json_file(filename))
        io = std::make_unique<JSON>();
    else {
        std::cout << "umr::mesh::io::dump: Only .vtk and .json file extensions accepted.\n";
        exit(-1);
    }

    if (prepend)
        prepend_filename(filename);

    io->set_cells(mesh);

    // could just add get_sources, etc. to MeshUnstructured and return empty sets....
    const MeshUnstructuredConstrained* m =
        dynamic_cast<const MeshUnstructuredConstrained*>(&mesh);
    if (m)
        io->set_sources(*m);

    io->write(filename);
}


} // namespace io

} // namespace mesh

} // namespace umr
