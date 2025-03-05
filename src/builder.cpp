#include "builder.hpp"

#include "mesh.hpp"
#include "inequalities.hpp"
#include "refinement.hpp"


namespace umr {


Builder::Builder() : m_commandline(nullptr) {}


Builder::Builder(int argc, char* argv[])
        : m_commandline(new ArgParser(argc, argv)) {}


Builder::~Builder() {
    if (m_commandline)
        delete m_commandline;
}


void Builder::set_density(density::DensityManager* density) {
    m_density = density;
}


void Builder::set_sources(source::SourceManager* sources) {
    m_sources = sources;
}


void Builder::set_min_angle(double angle) {
    assert(inequalities::is_lt(angle, 60)
            && "Builder::set_min_angle: Minimum angle must be less than 60 degrees");
    assert(inequalities::is_gt(angle, 0)
            && "Builder::set_min_angle: Minimum angle must be positive");

    m_min_angle = angle;
}


void Builder::set_refinement_none() {
    m_refinement = RefinementAlgorithm::None;
}


void Builder::set_refinement_chew_uniform() {
    m_refinement = RefinementAlgorithm::ChewUniform;
}


void Builder::set_refinement_chew_nonuniform() {
    m_refinement = RefinementAlgorithm::ChewNonuniform;
}


void Builder::set_refinement_ruppert() {
    m_refinement = RefinementAlgorithm::Ruppert;
}


void Builder::set_terminator_on() {
    m_use_terminator = true;
}


void Builder::set_terminator_off() {
    m_use_terminator = false;
}


MeshGenerator Builder::build() {
    update_with_command_line_arguments();
    check_ready_to_build();

    m_sources->finalize();
    std::vector<Real2> point_sources = m_sources->get_points();
    std::list<parametric::IParametric*> segment_sources = m_sources->get_segments();

    mesh::algo::IRefinement* refinement;
    mesh::MeshUnstructuredConstrained* mesh;

    switch (m_refinement) {
    case RefinementAlgorithm::None: {
        assert(point_sources.size() >= 3
                && "At least 3 source points required for Delaunay triangularization.");
        refinement = new mesh::algo::Refinement;
        mesh = new mesh::MeshUnstructuredConstrained;
        break;
    }
    case RefinementAlgorithm::ChewUniform: {
        refinement = new mesh::algo::RefinementChewUniform(m_density);
        mesh::algo::TrianglePQ<mesh::algo::TriCompMaxSize>* const triangle_pq =
            new mesh::algo::TrianglePQ<mesh::algo::TriCompMaxSize>();
        mesh = new mesh::MeshUnstructuredConstrained(triangle_pq);
        break;
    }
    case RefinementAlgorithm::ChewNonuniform: {
        refinement = new mesh::algo::RefinementChewNonuniform(m_min_angle, m_density, m_use_terminator);
        mesh::algo::TrianglePQ<mesh::algo::TriCompMinAngleMin>* const triangle_pq =
            new mesh::algo::TrianglePQ<mesh::algo::TriCompMinAngleMin>();
        mesh = new mesh::MeshUnstructuredConstrained(triangle_pq);
        break;
    }
    case RefinementAlgorithm::Ruppert: {
        refinement = new mesh::algo::RefinementRuppert(m_min_angle, m_density, m_use_terminator);
        mesh::algo::TrianglePQ<mesh::algo::TriCompMinAngleMin>* const triangle_pq =
            new mesh::algo::TrianglePQ<mesh::algo::TriCompMinAngleMin>();
        mesh = new mesh::MeshUnstructuredConstrained(triangle_pq);
        break;
    }
    default:
        std::cout << "Refinement algorithm not set!\n";
        exit(1);
    }

    return MeshGenerator(mesh, std::move(point_sources), std::move(segment_sources),
            refinement, m_data_frequency);
}


void Builder::update_with_command_line_arguments() {
    if (m_commandline == nullptr)
        return;

    std::cout << "====================================================\n";
    std::cout << "Command line argument overwrites\n";

    auto [density, density_set] = m_commandline->get_density();
    if (density_set) {
        if (m_density == nullptr) {
            std::cout << "Builder::build: Constructing default density\n";
            m_density = new density::DensityManager();
        }
        std::cout << "Builder::build: Minimum density set to " << density << " from the command line\n";
        m_density->set_minimum_density(density);
    }

    auto [angle, angle_set] = m_commandline->get_angle();
    if (angle_set) {
        std::cout << "Builder::build: Minimum angle set to " << angle << " deg from the command line\n";
        m_min_angle = angle;
    }

    auto [data_frequency, data_frequency_set] = m_commandline->get_data_frequency();
    if (data_frequency_set > 0) {
        std::cout << "Builder::build: Dumping mesh data every " << data_frequency << " iterations\n";
        m_data_frequency = data_frequency;
    }

    bool skip_terminator_set = m_commandline->get_skip_terminator();
    if (skip_terminator_set) {
        std::cout << "Builder::build: Not using the Terminator algorithm\n";
        m_use_terminator = false;
    }

    auto [refinement, refinement_set] = m_commandline->get_refinement_algorithm();
    if (refinement_set) {
        std::cout << "Builder::build: Refinement algorithm set to " << refinement << " from the command line\n";
        if (refinement == "chew-uniform")
            m_refinement = RefinementAlgorithm::ChewUniform;
        else if (refinement == "chew-nonuniform")
            m_refinement = RefinementAlgorithm::ChewNonuniform;
        else if (refinement == "ruppert")
            m_refinement = RefinementAlgorithm::Ruppert;
        else {
            std::cout << "Builder::update_with_command_line_arguments: refinement algorithm option "
                      << refinement << " not setup here\n";
            exit(1);
        }
    }

    std::cout << "====================================================\n";
}


void Builder::check_ready_to_build() {
    if (m_sources == nullptr) {
        std::cout << "Builder::build: Must set sources before building.\n";
        exit(1);
    }

    if (m_density == nullptr && m_refinement != RefinementAlgorithm::None) {
        std::cout << "Builder::build: Must set density before building.\n";
        exit(1);
    }
}


} // namespace umr
