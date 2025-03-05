#include <fstream>

#include "../ext/alglib-cpp/src/interpolation.h"

#include <cassert>
#include <bits/stdc++.h>

#include "density.hpp"
#include "inequalities.hpp"


namespace umr {

namespace density {


DensityInterpolator::DensityInterpolator(std::string& filename)
        : m_filename(filename) {
    initialize();
}


DensityInterpolator::DensityInterpolator(std::string&& filename)
        : m_filename(filename) {
    initialize();
}


double DensityInterpolator::interpolate(double x, double y) {
    return alglib::rbfcalc2(m_model, x, y);
}


std::string& DensityInterpolator::get_filename() const {
    return m_filename;
}


void DensityInterpolator::initialize() {
    load_data();
    build_model();
}


void DensityInterpolator::load_data() {
    std::ifstream file(m_filename);
    if (!file.is_open()) {
        std::cout << "DensityInterpolator::load_data(): "
            "Could not open file " << m_filename << '\n';
    }
    std::string line;
    std::vector<std::vector<double>> data;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string s;
        std::vector<double> row;
        while (std::getline(ss, s, ','))
            row.push_back(std::stod(s));
        data.push_back(row);
    }

    m_data.setlength(data.size(), data[0].size());
    for (int i = 0; i < m_data.rows(); ++i) {
        for (int j = 0; j < m_data.cols(); ++j) {
            m_data[i][j] = data[i][j];
        }
    }
}


void DensityInterpolator::build_model() {
    // 2: domain dimension
    // 1: range dimension
    alglib::rbfcreate(2, 1, m_model);
    alglib::rbfsetpoints(m_model, m_data);

    double rbase = 1.0;
    int n_layers = 3;
    double lambda_reg = 0.0; // no smoothing
    alglib::rbfsetalgohierarchical(m_model, rbase, n_layers, lambda_reg);

    alglib::rbfbuildmodel(m_model, m_rep);
    assert(m_rep.terminationtype == 1);
}


DensityManager::~DensityManager() {
    for (DensityInterpolator* interpolator : m_interpolators)
        delete interpolator;
}


void DensityManager::set_minimum_density(double d) {
    if (inequalities::is_le(d, 0)) {
        std::cout << "Minimum density must be positive!\n";
        exit(-1);
    }
    m_min_density = d;
    m_min_set = true;
}


void DensityManager::set_no_minimum_density() {
    m_min_density = std::numeric_limits<double>::max();
    m_min_set = false;
}


double DensityManager::get_minimum_density() const {
    return m_min_density;
}


void DensityManager::add_functor(DensityFunction function) {
    m_functors.push_back(function);
}


void DensityManager::add_interpolator(std::string& filename) {
    DensityInterpolator* interpolator = new DensityInterpolator(filename);
    m_interpolators.push_back(interpolator);
}


void DensityManager::add_interpolator(std::string&& filename) {
    DensityInterpolator* interpolator = new DensityInterpolator(filename);
    m_interpolators.push_back(interpolator);
}


double DensityManager::evaluate(double x, double y) const {
    double hmin = m_min_density;
    if (m_functors.size() > 0 || m_interpolators.size() > 0) {
        double hmin_functors = std::numeric_limits<double>::max();
        for (auto functor : m_functors)
            hmin_functors = std::min(hmin_functors, functor(x, y));

        double hmin_interpolators = std::numeric_limits<double>::max();
        for (DensityInterpolator* interpolator : m_interpolators)
            hmin_interpolators =
                std::min(hmin_interpolators, interpolator->interpolate(x, y));

        double calc_hmin = std::min(hmin_functors, hmin_interpolators);
        hmin = m_min_set ? std::max(hmin, calc_hmin) : calc_hmin;
    }

    if (inequalities::is_le(hmin, 0)) {
        std::cerr << "ERROR: Density must be positive! "
            "It's negative at x=" << x << ", y=" << y << '\n';
        if (!m_min_set && (m_functors.size() || m_interpolators.size()))
            std::cerr << "Minimum density not set!\n";
        exit(-1);
    }

    return hmin;
}


} // density

} // umr
