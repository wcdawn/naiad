#ifndef NAIAD_RESULT_HPP
#define NAIAD_RESULT_HPP

#include <vector>

class Result
{
  public:

    Result() {}
    Result(const std::vector<std::vector<double>> & phi_, const double & keff_, const int pnorder_ = 0)
      : phi{phi_}, keff{keff_}, pnorder{pnorder_}
    {}

    // index:
    // [ngroup][i * (pnorder+1) + n] for spatial index i and pn index n
    std::vector<std::vector<double>> phi;
    double keff;
    int pnorder;
};

#endif
