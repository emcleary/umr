#ifndef BUILDER_H
#define BUILDER_H

#include "mesh_generator.hpp"
#include "sources.hpp"
#include "density.hpp"
#include "command_line_interface.hpp"


namespace umr {


/**
 * @class RefinementAlgorithm
 *
 * An enum class used for tracking the refinement algorithm selected.
 */
enum class RefinementAlgorithm {
    None,              /**< No algorithm */
    ChewUniform,       /**< Chew Uniform algorithm */
    ChewNonuniform,    /**< Chew Nonuniform algorithm */
    Ruppert,           /**< Ruppert algorithm */
};


/**
 * @class Builder
 * @brief Builds a MeshGenerator
 *
 * The Builder takes settings and constructs a MeshGenerator object.
 * Settings can be set with the object's setters and/or through
 * command line arguments. Command line arguments will override
 * any previouly set settings.
 */
class Builder {
public:

    /**
     * @brief Creates a Builder without command line arguments.
     */
    Builder();

    /**
     * @brief Creates a Builder with command line arguments.
     *
     * @param argc Number of command line arguments.
     * @param argv Array of command line arguments.
     */
    Builder(int argc, char* argv[]);
    ~Builder();

    /** Set the density with a DensityManager. */
    void set_density(density::DensityManager* density);

    /** Set the sources with a SourceManager. */
    void set_sources(source::SourceManager* sources);

    /** Set the target minimum angle (degrees). */
    void set_min_angle(double angle);

    /** Set no refinement algorithm. */
    void set_refinement_none();

    /** Set the refinement algorithm to Chew Uniform. */
    void set_refinement_chew_uniform();

    /** Set the refinement algorithm to Chew Nonuniform. */
    void set_refinement_chew_nonuniform();

    /** Set the refinement algorithm to Ruppert. */
    void set_refinement_ruppert();

    /**
     * @brief Turn on the Terminator algorithm.
     *
     * Turns on the Terminator algorithm. Only applicable to the Chew
     * Nonuniform and Ruppert refinement algorithms.
     */
    void set_terminator_on();

    /**
     * @brief Turn off the Terminator algorithm.
     *
     * Turns off the Terminator algorithm. Only applicable to the Chew
     * Nonuniform and Ruppert refinement algorithms.
     */
    void set_terminator_off();

    /** Constructs a MeshGenerator object. */
    MeshGenerator build();

private:
    /** Overrides settings with command line arguments. */
    void update_with_command_line_arguments();

    /** Confirms mandatory settings are in place. */
    void check_ready_to_build();

    ArgParser* const m_commandline;
    RefinementAlgorithm m_refinement = RefinementAlgorithm::None;
    source::SourceManager* m_sources = nullptr;
    density::DensityManager* m_density = nullptr;
    double m_min_angle = 0;
    int m_data_frequency = 0;
    bool m_use_terminator = true;
};


} // namespace umr

#endif // BUILDER_H
