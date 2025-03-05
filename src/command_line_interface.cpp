#include "command_line_interface.hpp"

#include <iostream>
#include <ranges>
#include <set>


namespace umr {


static const std::set<std::string> allowed_algorithms =
    {"chew-uniform", "chew-nonuniform", "ruppert"};


static const std::string allowed_strings = "refinement algorithms: "
    + (allowed_algorithms | std::views::join_with(std::string(", "))
            | std::ranges::to<std::string>());


ArgParser::ArgParser(int argc, char* argv[]) : m_elf(argv[0]) {
    po::options_description description("Usage: " + m_elf + " [options]\n\nAllowed options");
    description.add_options()("help,h", "produce help message");
    description.add_options()("density,d", po::value<double>(), "minimum allowed density");
    description.add_options()("angle,a", po::value<double>(), "minimum allowed angle (deg)");
    description.add_options()("data-frequency,f", po::value<int>(), "dump data every N iterations");
    description.add_options()("no-terminator", "skip the terminator algorithm");

    auto refinement_validator = [&](const std::string& value) {
        if (!allowed_algorithms.contains(value))
            throw po::error("Invalid refinement algorithm option: " + value);
    };
    description.add_options()("refinement,r",
            po::value<std::string>()->notifier(refinement_validator),
            allowed_strings.c_str());

    m_options.add(description);

    try {
        po::store(po::parse_command_line(argc, argv, m_options), m_varmap);
        po::notify(m_varmap);
    } catch (const po::error& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_help();
    }

    if (m_varmap.count("help"))
        print_help();
}


std::pair<double, bool> ArgParser::get_angle() const {
    if (m_varmap.count("angle"))
        return {m_varmap["angle"].as<double>(), true};
    return {-1, false};
}


std::pair<double, bool> ArgParser::get_density() const {
    if (m_varmap.count("density"))
        return {m_varmap["density"].as<double>(), true};
    return {-1, false};
}


std::pair<int, bool> ArgParser::get_data_frequency() const {
    if (m_varmap.count("data-frequency"))
        return {m_varmap["data-frequency"].as<int>(), true};
    return {-1, false};
}


bool ArgParser::get_skip_terminator() const {
    return m_varmap.count("no-terminator");
}


std::pair<std::string, bool> ArgParser::get_refinement_algorithm() const {
    if (m_varmap.count("refinement"))
        return {m_varmap["refinement"].as<std::string>(), true};
    return {"", false};
}


void ArgParser::print_help() const {
    std::cout << m_options << '\n';
    exit(1);
}


} // namespace umr
