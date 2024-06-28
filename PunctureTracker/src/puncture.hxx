/* puncture.hxx */
/* (c) Liwei Ji 06/2024 */

#ifndef PUNCTURETRACKER_PUNCTURE_HXX
#define PUNCTURETRACKER_PUNCTURE_HXX

#include <loop_device.hxx>

#include <cctk.h>

#include <vector>

namespace PunctureTracker {

class PunctureContainer {
public:
  PunctureContainer() {}

  virtual ~PunctureContainer() = default; // Use default for trivial destructors

  void setNumPunctures() { numPunctures_ = location_[0].size(); }

  CCTK_INT getNumPunctures() { return numPunctures_; }

  std::vector<CCTK_REAL> &getTime() { return time_; }

  std::vector<CCTK_REAL> &getPreviousTime() { return previousTime_; }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getLocation() {
    return location_;
  }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getVelocity() {
    return velocity_;
  }

  std::array<std::vector<CCTK_REAL>, Loop::dim> &getBeta() { return beta_; }

  void updatePreviousTime(CCTK_ARGUMENTS);

  void interpolate(CCTK_ARGUMENTS);

  void evolve(CCTK_ARGUMENTS);

  void broadcast(CCTK_ARGUMENTS);

private:
  CCTK_INT numPunctures_;
  std::vector<CCTK_REAL> time_, previousTime_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> location_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> velocity_;
  std::array<std::vector<CCTK_REAL>, Loop::dim> beta_;
};

} // namespace PunctureTracker

#endif // #ifndef PUNCTURETRACKER_PUNCTURE_HXX
