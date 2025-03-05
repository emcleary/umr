#include <cmath>

#include "point.hpp"


namespace umr {

namespace mesh {


Point::Point() : x(0), y(0) {}


Point::Point(double x, double y) : x(x), y(y) {}


Point::Point(std::pair<double, double> p) : x(p.first), y(p.second) {}


Point::Point(const Point& p) : x(p.x), y(p.y) {}


Point::Point(const Point&& p) : x(p.x), y(p.y) {}


double Point::distance_to(const Point& p) const {
    double dx = x - p.x;
    double dy = y - p.y;
    return sqrt(dx*dx + dy*dy);
}


double Point::distance_to(const double _x, const double _y) const {
    double dx = x - _x;
    double dy = y - _y;
    return sqrt(dx*dx + dy*dy);
}


bool Point::operator==(const Point& p) const {
    return this->x == p.x && this->y == p.y;
}


bool Point::operator<(const Point& p) const {
    if (this->x == p.x)
        return this->y < p.y;
    return this->x < p.x;
}


std::ostream& operator<<(std::ostream& out, const Point& p) {
    out << "(" << std::format("{}", +p.x) << ", " << std::format("{}", +p.y) << ")";
    return out;
}


size_t PointHash::operator()(const Point& p) const {
    size_t h = 17;
    h += h * 31 + std::hash<double>()(p.x);
    h += h * 31 + std::hash<double>()(p.y);
    return h;
}


} // namespace mesh

} // namespace umr


size_t std::hash<umr::mesh::Point>::operator()(const umr::mesh::Point& p) {
    size_t h = 17;
    h += h * 31 + std::hash<double>()(p.x);
    h += h * 31 + std::hash<double>()(p.y);
    return h;
}
