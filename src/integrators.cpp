#include "integrators.hpp"

#include "parametrics.hpp"


namespace umr {

namespace math {

namespace integrate {


void segment_integrand(double t, double t_minus_tmin, double tmax_minus_t,
        double &result, void *ptr) {
    umr::parametric::IParametric* parametric =
        static_cast<umr::parametric::IParametric*>(ptr);
    auto [dxdt, dydt] = parametric->differentiate_model(t);
    result = sqrt(dxdt * dxdt + dydt * dydt);
}


void SegmentLength::initialize(void* ptr, double& tmin, double& tmax) {
    umr::parametric::IParametric* parametric =
        static_cast<umr::parametric::IParametric*>(ptr);
    tmin = parametric->get_tmin();
    tmax = parametric->get_tmax();
}


} // namespace integrate

} // namespace math

} // namespace umr
