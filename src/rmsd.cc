#include <iostream>
#include <vector>
#include <cmath>

double rmsd(const std::vector<double> vec1, const std::vector<double> vec2) {
    if (vec1.size() != vec2.size()) {
        std::cerr << "rmsd: vectors have different size!" << std::endl;
        exit(1);
    }
    double result = 0.0;
    for (unsigned int i = 0; i < vec1.size(); ++i) {
        result += std::pow(vec1.at(i) - vec2.at(i), 2);
    }
    result /= double(vec1.size());
    result = std::sqrt(result);

    return result;
}
