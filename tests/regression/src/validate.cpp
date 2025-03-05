#include <iostream>
#include <fstream>
#include <json/json.h>

#include "common.hpp"
#include "../../../src/utilities.hpp"
#include "../../../src/mesh.hpp"
#include "../../../src/io.hpp"
#include "PointN.hpp"
#include "ADTN.hpp"


struct Bounds {
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = 0;
    double ymax = 0;
};


class Files {
public:
    Files(const umr::mesh::MeshUnstructured& mesh, unsigned int index)
            : m_mesh(mesh), m_index(index) {
        m_filename = m_index == 0 ? "gold_standard.json"
            : std::format("gold_standard_{}.json", m_index);
        m_outfile = m_index == 0 ? "output.json"
            : std::format("output_{}.json", m_index);
        m_outfile_vtk = m_index == 0 ? "output.vtk"
            : std::format("output_{}.vtk", m_index);
    }

    bool initialize() {
        if (!load_data())
            return false;
        load_points();
        return true;
    }

    unsigned int get_points_size() {
        return m_points.size();
    }

    unsigned int get_edges_size() {
        return m_data["edges"].size();
    }

    void failed() {
        umr::mesh::io::dump(m_mesh, m_outfile, false);
        umr::mesh::io::dump(m_mesh, m_outfile_vtk, false);
    }

    ADTN get_adtn_points() {
        PointN pmin({m_bounds.xmin, m_bounds.ymin});
        PointN pmax({m_bounds.xmax, m_bounds.ymax});
        ADTN gs2(pmin, pmax);

        for (auto p : m_points)
            gs2.add_point(p);

        return gs2;
    }

    ADTN get_adtn_edges() {
        PointN pmin4({m_bounds.xmin, m_bounds.ymin, m_bounds.xmin, m_bounds.ymin});
        PointN pmax4({m_bounds.xmax, m_bounds.ymax, m_bounds.xmax, m_bounds.ymax});
        ADTN gs4(pmin4, pmax4);
    
        for (Json::Value::ArrayIndex i = 0; i != m_data["edges"].size(); ++i) {
            unsigned int j = m_data["edges"][i][0].asUInt();
            unsigned int k = m_data["edges"][i][1].asUInt();
            PointN p4({m_points[j][0], m_points[j][1], m_points[k][0], m_points[k][1]});
            PointN p4s({m_points[k][0], m_points[k][1], m_points[j][0], m_points[j][1]});
            gs4.add_point(p4);
            gs4.add_point(p4s);
        }

        return gs4;
    }
    
    
private:
    bool load_data() {
        std::ifstream input_file;
        input_file.open(m_filename);
        if (!input_file) {
            std::cout << "could not find " << m_filename << '\n';
            std::cout << "writing current mesh to file\n";
            failed();
            return false;
        }

        input_file >> m_data;
        input_file.close();
        return true;
    }

    void load_points() {
        m_points.resize(m_data["points"].size());
        for (Json::Value::ArrayIndex i = 0; i != m_data["points"].size(); ++i) {
            double x = m_data["points"][i][0].asDouble();
            double y = m_data["points"][i][1].asDouble();
            m_points[i] = PointN({x, y});
            m_bounds.xmin = std::min(m_bounds.xmin, x);
            m_bounds.ymin = std::min(m_bounds.ymin, y);
            m_bounds.xmax = std::max(m_bounds.xmax, x);
            m_bounds.ymax = std::max(m_bounds.ymax, y);
        }
    }

    const umr::mesh::MeshUnstructured& m_mesh;
    const unsigned int m_index;
    std::string m_filename, m_outfile, m_outfile_vtk;
    Json::Value m_data;
    std::vector<PointN> m_points;
    Bounds m_bounds;
};


bool test_points(const umr::mesh::MeshUnstructured& mesh, Files& files) {
    ADTN gs_points = files.get_adtn_points();

    unsigned int counter_pts = 0;
    for (auto p : mesh.get_points_set()) {
        PointN pn({p->x, p->y});
        counter_pts += !gs_points.contains(pn);
    }

    if (counter_pts == 0)
        std::cout << "All json points exist in the mesh!\n";
    else {
        std::cout << "Missing " << counter_pts << " points\n";
        files.failed();
        return false;
    }

    if (files.get_points_size() != mesh.get_points_set().size()) {
        std::cout << "Json and mesh points have different sizes!\n";
        files.failed();
        return false;
    }

    std::cout << "Valid points\n";
    return true;
}


bool test_edges(const umr::mesh::MeshUnstructured& mesh, Files& files) {
    ADTN gs_edges = files.get_adtn_edges();

    unsigned int counter_edges = 0;
    for (auto e : mesh.get_edges_set()) {
        PointN pn({e->org().x, e->org().y, e->dest().x, e->dest().y});
        counter_edges += !gs_edges.contains(pn);
    }

    if (counter_edges == 0)
        std::cout << "All json edges exist in the mesh!\n";
    else {
        std::cout << "Missing " << counter_edges << " edges\n";
        files.failed();
        return false;
    }

    if (files.get_edges_size() != mesh.get_edges_set().size()) {
        std::cout << "Json and mesh edges have different sizes!\n";
        files.failed();
        return false;
    }

    std::cout << "Valid edges\n";
    return true;
}


bool validate_mesh(const umr::mesh::MeshUnstructured& mesh, unsigned int index) {
    Files files(mesh, index);
    if (!files.initialize())
        return false;
    
    bool points_passed = test_points(mesh, files);
    bool edges_passed = test_edges(mesh, files);
    if (points_passed && edges_passed) {
        std::cout << "Valid mesh!\n";
        return true;
    }

    return false;
}


bool validate_points(const umr::mesh::MeshUnstructured& mesh, unsigned int index) {
    Files files(mesh, index);
    if (!files.initialize())
        return false;
    
    bool points_passed = test_points(mesh, files);
    if (points_passed) {
        std::cout << "Valid points\n";
        return true;
    }

    return false;
}
