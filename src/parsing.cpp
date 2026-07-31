#include "include/parsing.h"

frd parsing_frd(const std::string file_name) {
    frd file_info;
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open FRD file " << file_name << std::endl;
        return file_info;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // skip header lines
        size_t first_char = line.find_first_not_of(" \t\r\n");
        if (first_char == std::string::npos || line[first_char] == '*' || line[first_char] == ';') {
            continue;
        }
        
        std::stringstream ss(line);
        double f, db, deg; 

        if (ss >> f >> db) {
            file_info.freq.push_back(f);
            file_info.spl.push_back(db);

            if (ss >> deg) {
                file_info.phase.push_back(deg);
            } else {
                file_info.phase.push_back(0.0);
            }
        }
    }
    return file_info;
}

zma parsing_zma(const std::string file_name) {
    zma file_info;
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open ZMA file " << file_name << std::endl;
        return file_info;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // skip header lines
        size_t first_char = line.find_first_not_of(" \t\r\n");
        if (first_char == std::string::npos || line[first_char] == '*' || line[first_char] == ';') {
            continue;
        }

        std::stringstream ss(line);
        double f, imp, deg;

        if (ss >> f >> imp) {
            file_info.freq.push_back(f);
            file_info.imp.push_back(imp);

            if (ss >> deg) {
                file_info.phase.push_back(deg);
            } else {
                file_info.phase.push_back(0.0);
            }
        }
    }
    return file_info;
}
