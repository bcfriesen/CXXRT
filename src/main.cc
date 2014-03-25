#include <iostream>
#include <fstream>
#include <exception>

#include <yaml-cpp/yaml.h>

using namespace std;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        cout << "Usage: <executable name> <YAML control file>" << endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(argv[1]);
}
