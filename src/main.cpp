#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "input.hpp"
#include "output.hpp"
#include "geometry.hpp"
#include "xslibrary.hpp"
#include "exception_handler.hpp"
#include "diffusion.hpp"
#include "result.hpp"
#include "writer.hpp"
#include "analysis.hpp"
#include "quadrature1d.hpp"
#include "quadrature_gauss_legendre.hpp"

using namespace naiad;

std::string get_stub(const std::string & fname)
{
  const auto last{fname.find_last_of('.')};
  return fname.substr(0,last); // does not include period
}

int main(int argc, char* argv[])
{

  // process arguments
  std::vector<std::string> args;
  args.reserve(argc);
  for (int i = 0; i < argc; ++i)
    args.emplace_back(argv[i]);

  if (args.size() != std::size_t{2})
  {
    std::cerr << "Expect exactly one argument (filename)." << std::endl;
    return 1;
  }

  const std::string fname_inp{args[1]};
  const std::string fname_stub{get_stub(fname_inp)};
  const std::string fname_out{fname_stub + ".out"};

  std::ofstream fout{fname_out};

  naiad::out.link_stream(std::cout);
  naiad::out.link_stream(fout);

  naiad::out << "BEGIN naiad ναϊάς" << std::endl;
  naiad::out << std::endl;

  // print arguments to terminal
  naiad::out << "=== COMMAND-LINE ARGUMENTS ===" << std::endl;
  for (const auto & x : args)
    naiad::out << x << std::endl;
  naiad::out << std::endl;

  // parse input
  const Input input{fname_inp};

  input.echo(fout);
  input.summary(naiad::out);

  const XSLibrary & xslib{input.xslibrary()};
  xslib.summary(naiad::out);

  // this is actually a copy so that I can do refinement
  Geometry geo{input.geometry()};

  naiad::out << "(before refinement)" << std::endl;
  geo.summary(naiad::out);

  for (int i = 0; i < input.refine; ++i)
    geo.refine();

  if (input.refine > 0)
  {
    naiad::out << "(after refinement)" << std::endl;
    geo.summary(naiad::out);
  }

  Result res;
  if (input.snorder == 0)
  {
    // diffusion
    const Diffusion_solver diffusion{geo, input.bc_left, input.bc_right, xslib, input.tolerance()};
    res = diffusion.solve();
  }
  else
  {
    // transport
    const std::unique_ptr<Quadrature_gauss_legendre> quadrature{std::make_unique<Quadrature_gauss_legendre>(input.snorder)};
  }

  naiad::out << "keff = " << std::format("{:.20f}", res.keff) << std::endl << std::endl;

  std::unique_ptr<Analysis> analysis{nullptr};
  switch (input.analysis_reference)
  {
    case (Analysis_reference::critical):
      analysis = std::make_unique<Analysis_critical>(geo, xslib, res);
      break;
    case (Analysis_reference::onegroup):
      analysis = std::make_unique<Analysis_onegroup>(geo, xslib, res);
      break;
    case (Analysis_reference::twogroup):
      analysis = std::make_unique<Analysis_twogroup>(geo, xslib, res);
    case (Analysis_reference::tworegion):
      analysis = std::make_unique<Analysis_tworegion>(geo, xslib, res);
    default:
      // do nothing
      break;
  }
  if (analysis)
    analysis->summary(naiad::out);

  const Writer writer{geo, res};
  const std::string fname_flux_csv{fname_stub + "_flux.csv"};
  naiad::out << "writing flux csv on " << fname_flux_csv << std::endl;
  writer.write_flux(fname_flux_csv);
  naiad::out << std::endl;

  exception.summary(naiad::out);

  naiad::out << "END naiad ναϊάς" << std::endl;
  naiad::out << "Normal termination :)" << std::endl;

  return 0;
}
