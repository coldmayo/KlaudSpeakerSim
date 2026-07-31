#ifndef parsing_h_INCLUDED
#define parsing_h_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct frd {
    std::vector<double> freq;
    std::vector<double> spl;
    std::vector<double> phase;
};

struct zma {
    std::vector<double> freq;
    std::vector<double> imp;
    std::vector<double> phase;
};

frd parsing_frd(const std::string file_name);
zma parsing_zma(const std::string file_name);

#endif // parsing_h_INCLUDED
