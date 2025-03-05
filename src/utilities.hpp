#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath> // abs
#include <functional> // less, hash
#include <cassert>
#include <string>


template<typename T>
struct pointer_less_than {
    bool operator()(T* lhs, T* rhs) const {
        return (lhs && rhs)
            ? std::less<T>()(*lhs, *rhs)
            : std::less<T*>()(lhs, rhs);
    }
};


template<typename T>
struct pointer_equal {
    bool operator()(T* lhs, T* rhs) const {
        return (lhs && rhs) ? *lhs == *rhs : lhs == rhs;
    }
};


template <typename T>
struct pointer_hash {
    size_t operator()(T* var) const {
        return var ? std::hash<T>{}(*var) : std::hash<T*>{}(var);
    }
};


bool is_vtk_file(std::string& str);
bool is_json_file(std::string& str);

#endif // UTILITIES_H
