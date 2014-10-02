#include <iostream>
#include <vector>
#include <cmath>

#include <viennacl/vector.hpp>

double calc_rmsd(const viennacl::vector<double> vec1, const viennacl::vector<double> vec2) {
    if (vec1.size() != vec2.size()) {
        std::cerr << "rmsd: vectors have different size!" << std::endl;
        exit(1);
    }
    std::vector<double> relative_change(vec1.size());

    for (unsigned int i = 0; i < relative_change.size(); ++i) {
        relative_change.at(i) = (vec2(i) - vec1(i)) / vec1(i);
    }

    double result = 0.0;
    for (unsigned int i = 0; i < vec1.size(); ++i) {
        result += std::pow(relative_change.at(i), 2);
    }
    result /= double(vec1.size());
    result = std::sqrt(result);

    return result;
}
