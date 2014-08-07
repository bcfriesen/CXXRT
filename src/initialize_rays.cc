#include "globals.hh"
#include "ray.hh"

void initialize_rays() {

    const double epsilon = config["epsilon"].as<double>();
    const int n_mu_pts = config["n_mu_pts"].as<int>();

    rays.resize(n_mu_pts);
    const double mu_min = -1.0;
    const double mu_max = +1.0;
    double mu = mu_min;
    for (Ray& r: rays) {
        if (std::fabs(mu) < std::numeric_limits<double>::epsilon()) {
            std::cerr << "ERROR: mu too close to zero! : " << mu << std::endl;
            exit(1);
        }
        r.bind_to_grid(mu);
        for (RayData& rd: r.raydata) {
            for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
                rwlp.epsilon = epsilon;
            }
        }
        mu += (mu_max - mu_min) / double(rays.size()-1);

        for (RayData& rd: r.raydata) {
            for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
                rwlp.calc_chi(rd.gridvoxel->rho, *(rwlp.lambda));
            }
        }

        for (auto wlv: wavelength_values) {
            r.calc_tau(wlv);
            r.calc_SC_coeffs(wlv);
        }

        for (RayData& rd: r.raydata) {
            for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
                rwlp.set_to_LTE(rd.gridvoxel->temperature);
            }
        }

        for (auto wlv: wavelength_values) {
            r.formal_soln(wlv);
        }
    }

    for (GridVoxel& gv: grid) {
        for (GridWavelengthPoint& wlp: gv.wavelength_grid) {
            gv.calc_J(*(wlp.lambda));
            gv.calc_H(*(wlp.lambda));
            gv.calc_K(*(wlp.lambda));
        }
    }
}
