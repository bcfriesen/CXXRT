#include "globals.hh"
#include "ray.hh"

void initialize_rays() {

    const int n_mu_pts = config["n_mu_pts"].as<int>();

    rays.resize(n_mu_pts);
    const double mu_min = -1.0;
    const double mu_max = +1.0;
    double mu = mu_min;
    for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
        if (std::fabs(mu) < std::numeric_limits<double>::epsilon()) {
            std::cerr << "ERROR: mu too close to zero! : " << mu << std::endl;
            exit(1);
        }
        r->bind_to_grid(mu);
        mu += (mu_max - mu_min) / double(rays.size()-1);
    }
}
