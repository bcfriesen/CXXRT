#include <iostream>
#include <fstream>
#include <exception>

#include <yaml-cpp/yaml.h>


int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(argv[1]);
}
