#include "ADTN.hpp"

#include "../../../src/inequalities.hpp"


ADTN::ADTN(PointN xmin, PointN xmax)
        : m_xmin(xmin), m_xmax(xmax), m_direction(0) {
    initialize();
}


ADTN::ADTN(PointN xmin, PointN xmax, unsigned int direction)
        : m_xmin(xmin), m_xmax(xmax), m_direction(direction) {
    initialize();
}


void ADTN::add_point(PointN& p) {
    assert(p.size() == m_xmin.size());

    ADTN* node = this;
    while (node && node->m_point) {
        if (p[node->m_direction] < node->m_xmid[node->m_direction]) {
            if (!node->m_left)
                node->m_left = std::make_unique<ADTN>(node->m_xmin, node->m_xmid,
                        (node->m_direction + 1) % m_xmin.size());
            node = node->m_left.get();
        } else {
            if (!node->m_right)
                node->m_right = std::make_unique<ADTN>(node->m_xmid, node->m_xmax,
                        (node->m_direction + 1) % m_xmin.size());
            node = node->m_right.get();
        }
    }
    
    node->m_point = std::make_unique<PointN>(p);
    return;
}


bool ADTN::contains(PointN& p) {
    assert(p.size() == m_xmin.size());
    ADTN* node = this;
    static const double atol = 1e-8;
    static const double rtol = 1e-6;
    while (node && node->m_point) {
        bool close = true;
        for (size_t i = 0; i < p.size(); i++)
            close &= umr::inequalities::is_close(p[i], (*node->m_point)[i], rtol, atol);

        if (close)
            return true;

        unsigned int dir = node->m_direction;
        node = p[dir] < node->m_xmid[dir] ? node->m_left.get() : node->m_right.get();
    }
    
    return false;
}


void ADTN::initialize() {
    std::vector<double> xmid;
    xmid.resize(m_xmin.size());
    for (size_t i = 0; i < m_xmin.size(); i++)
        xmid[i] = (m_xmin[i] + m_xmax[i]) / 2;
    m_xmid = PointN(xmid);
}
