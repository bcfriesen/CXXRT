#include "ray.hh"
#include "globals.hh"

void Lambda_iteration() {

    for (Ray& r: rays) {
      r.calc_source_fn();
      r.formal_soln();
    }

    for (GridVoxel& gv: grid) {
      gv.calc_J();
    }

}
