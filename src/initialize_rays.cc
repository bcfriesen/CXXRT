#include "globals.hh"
#include "ray.hh"

void initialize_rays() {

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
        mu += (mu_max - mu_min) / double(rays.size()-1);

        for (RayData& rd: r.raydata) {
            for (auto &rwlp: rd.wavelength_grid) {
                rwlp.second.calc_chi(rd.gridvoxel->rho, rwlp.second.lambda);
            }
        }

        for (auto wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
            r.calc_tau(wlv->first);
            r.calc_SC_coeffs(wlv->first);
        }

        for (RayData& rd: r.raydata) {
            for (auto &rwlp: rd.wavelength_grid) {
                rwlp.second.set_to_LTE(rd.gridvoxel->temperature);
            }
        }

        for (auto wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
            r.formal_soln(wlv->first);
        }
    }

    for (GridVoxel& gv: grid) {
        for (auto &wlp: gv.wavelength_grid) {
            gv.calc_J(wlp.first);
            gv.calc_H(wlp.first);
            gv.calc_K(wlp.first);

            // In the absence of anisotropic sources/sinks, the thermalization
            // parameter epsilon will be the same for every ray at their
            // respective points where they intersect this voxel. So just use
            // the first one to set the value of epsilon for the grid.
            auto rip = gv.ray_intersection_data.at(0).ray->raydata.begin() + gv.ray_intersection_data.at(0).intersection_point;
            wlp.second.epsilon = rip->wavelength_grid[wlp.first].epsilon;
        }
    }
}
