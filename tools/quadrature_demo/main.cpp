#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include "quadrature1d.hpp"
#include "quadrature_gauss_legendre.hpp"

double fun(double x)
{
  return std::cos(x);
}

int main(int argc, char * argv_in[])
{
  std::cout << "begin" << std::endl;
  std::cout << std::endl;

  // unpack argv
  std::vector<std::string> argv;
  argv.reserve(argc);
  for (int i = 0; i < argc; ++i)
    argv.emplace_back(argv_in[i]);
  std::cout << "ARGV:" << std::endl;
  for (const auto & arg : argv)
    std::cout << arg << std::endl;
  std::cout << std::endl;

  const double exact{2.0*std::sin(1.0)};

  std::unique_ptr<Quadrature1d> quad;
  for (int n = 1; n <= 64; ++n)
  {
    quad = std::make_unique<Quadrature_gauss_legendre>(n);
    const double x{quad->integrate(fun, -1.0, 1.0)};
    std::cout << std::format("n={:d} int={:.16e} err={:.16e}", n, x, exact-x) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "end" << std::endl;
  return 0;
}
