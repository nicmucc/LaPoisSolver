#pragma once

#include <string>
#include <sys/stat.h>

namespace poisson {

/// Parse output directory from argv, or fall back to default.
/// Creates the directory if it doesn't exist.
///
/// Usage in main():
///   std::string out = output_dir(argc, argv);  // default: "output/"
///   field_export::to_json(grid, out + "result.json");
///
/// From command line:
///   ./parallel_plates                  -> writes to output/
///   ./parallel_plates -o results/run1  -> writes to results/run1/
inline std::string output_dir(int argc, char* argv[], const std::string& default_dir = "output") {
    std::string dir = default_dir;

    for (int i = 1; i < argc - 1; ++i) {
        std::string arg(argv[i]);
        if (arg == "-o" || arg == "--output") {
            dir = argv[i + 1];
            break;
        }
    }

    // Ensure trailing slash
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') {
        dir += '/';
    }

    // Create directory (works on POSIX and MinGW)
#ifdef _WIN32
    _mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif

    return dir;
}

} // namespace poisson
