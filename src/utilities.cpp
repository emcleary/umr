#include "utilities.hpp"

#include <iostream>


bool has_suffix(std::string& str, std::string&& suffix) {
    if (str.size() > suffix.size()) {
        size_t start_index = str.size() - suffix.size();
        size_t length = suffix.size();
        return str.compare(start_index, length, suffix) == 0;
    }
    return false;
}


bool is_vtk_file(std::string& str) {
    return has_suffix(str, ".vtk");
}


bool is_json_file(std::string& str) {
    return has_suffix(str, ".json");
}
