#ifndef ADTN_H
#define ADTN_H

#include <memory>
#include <vector>

#include "PointN.hpp"


/**
 * @class ADTN
 * @brief N-dimensional Alternating Digital Tree algorithm
 *
 * This class is used for searching for points in a domain using the
 * ADT algorithm. It is implemented to work in N-dimensional spaces.
 * Currently it is only used for testing purposes.
 */
class ADTN {
public:
    /** Constructor for ADTN */
    ADTN(PointN xmin, PointN xmax);

    /** Constructor for ADTN */
    ADTN(PointN xmin, PointN xmax, unsigned int direction);

    /** Add a PointN object */
    void add_point(PointN& p);

    /** Check if a PointN object exists */
    bool contains(PointN& p);

private:
    void initialize();

    PointN m_xmin, m_xmax;
    const unsigned int m_direction = 0;
    PointN m_xmid;
    std::unique_ptr<PointN> m_point;
    std::unique_ptr<ADTN> m_left, m_right;
};

#endif // ADTN_H
