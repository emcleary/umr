#ifndef POINTN_H
#define POINTN_H

#include <ostream>
#include <initializer_list>
#include <vector>


/**
 * @class PointN
 * @brief N-dimensional Point in space
 *
 * This class is for N-dimensional points. It is implemented to work
 * with the ADTN class and is only intended for testing purposes.
 */
class PointN {
public:
    friend std::ostream& operator<<(std::ostream& os, PointN& p);
    PointN() {};

    PointN(std::vector<double>& x) {
        m_x.insert(m_x.begin(), x.begin(), x.end());
    }
    PointN(std::vector<double>&& x) {
        m_x.insert(m_x.begin(), x.begin(), x.end());
    }

    PointN(PointN& p) {
        m_x = p.m_x;
    }

    PointN(PointN&& p) {
        m_x = p.m_x;
    }

    PointN& operator=(PointN& p) {
        m_x = p.m_x;
        return *this;
    }

    PointN& operator=(PointN&& p) {
        m_x = p.m_x;
        return *this;
    }

    double operator[](unsigned int i) {
        return m_x[i];
    }

    unsigned int size() { return m_x.size(); }

private:
    std::vector<double> m_x;

};

std::ostream& operator<<(std::ostream& os, PointN& p);

#endif // POINTN_H
