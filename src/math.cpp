#include "math.hpp"

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

namespace naiad
{

double legendre(const int n, const double x)
{
#ifdef __STDCPP_MATH_SPEC_FUNCS__
  return std::legendre(n, x);
#else
  switch (n)
  {
    case (0):
      return 1.0;
    case (1):
      return x;
    case (2):
      return 0.5 * (3.0 * std::pow(x, 2) - 1.0);
    case (3):
      return 0.5 * (5.0 * std::pow(x, 3) - 3.0 * x);
    case (4):
      return 0.125 * (35.0 * std::pow(x, 4) - 30.0 * std::pow(x, 2) + 3.0);
    case (5):
      return 0.125 * (63.0 * std::pow(x, 5) - 70.0 * std::pow(x, 3) + 15.0 * x);
    case (6):
      return 0.0625 * (231.0 * std::pow(x, 6) - 315.0 * std::pow(x, 4) + 105.0 * std::pow(x, 2) - 5.0);
    case (7):
      return 0.0625 * (429.0 * std::pow(x, 7) - 693.0 * std::pow(x, 5) + 315.0 * std::pow(x, 3) - 35.0 * x);
    case (8):
      return (6435.0 * std::pow(x, 8) - 12012.0 * std::pow(x, 6) + 6930.0 * std::pow(x, 4) - 1260.0 * std::pow(x, 2)
              + 35.0)
             / 128.0;
    case (9):
      return (12155.0 * std::pow(x, 9) - 25740.0 * std::pow(x, 7) + 18018.0 * std::pow(x, 5) - 4620.0 * std::pow(x, 3)
              + 315.0 * x)
             / 128.0;
    case (10):
      return (46189.0 * std::pow(x, 10) - 109395.0 * std::pow(x, 8) + 90090.0 * std::pow(x, 6)
              - 30030.0 * std::pow(x, 4) + 3465.0 * std::pow(x, 2) - 63.0)
             / 256.0;
    default:
      return (2.0 * n - 1.0) * x * legendre(n - 1, x) - (n - 1.0) / n * legendre(n - 2, x);
  }
#endif
}

} // namespace naiad
