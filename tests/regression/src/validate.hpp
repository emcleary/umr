#pragma once
#include "../../../src/mesh.hpp"

bool validate_mesh(const umr::mesh::MeshUnstructured& m, unsigned int index = 0);
bool validate_points(const umr::mesh::MeshUnstructured& m, unsigned int index = 0);
