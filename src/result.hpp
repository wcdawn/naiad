#ifndef NAIAD_RESULT_HPP
#define NAIAD_RESULT_HPP

#include <vector>

class Result
{
  public:

    Result() {}
    Result(const std::vector<std::vector<double>> & phi_, const double & keff_)
      : phi{phi_}, keff{keff_}
    {}

    std::vector<std::vector<double>> phi;
    double keff;
};

#endif
