#ifndef COMMAND_LINE_INTERFACE_H
#define COMMAND_LINE_INTERFACE_H

#include <boost/program_options.hpp>


namespace umr {


namespace po = boost::program_options;


/**
 * @class ArgParser
 * @brief Command line argument parser
 *
 * This parses command line arguments. When used, (in Builder, for
 * example) run the executable with the argument "-h" or "--help" to
 * display current settings.
 */
class ArgParser {
public:

    /**
     * @brief Constructs the ArgParser
     *
     * Constructs the ArgParser taking in the standard "argc" and
     * "argv" variables in a "main" function.
     *
     * @param argc Number of command line arguments.
     * @param argv Array of command line arguments.
     */
    ArgParser(int argc, char* argv[]);
    ~ArgParser() {}

    /** Get the minimum angle plus a bool to track if set. */
    std::pair<double, bool> get_angle() const;

    /** Get the minimum density plus a bool to track if set. */
    std::pair<double, bool> get_density() const;

    /** Get the data dumping frequency plus a bool to track if set. */
    std::pair<int, bool> get_data_frequency() const;

    /** Get the refinement algorithm plus a bool to track if set. */
    std::pair<std::string, bool> get_refinement_algorithm() const;

    /** Get the skip Terminator algorithm flag (default to false). */
    bool get_skip_terminator() const;

private:
    void print_help() const;

    const std::string m_elf;
    po::options_description m_options;
    po::variables_map m_varmap;
};


} // namespace umr


#endif // COMMAND_LINE_INTERFACE_H
