#include <iostream>
#include <vector>
#include <string>

#include "input.hpp"
#include "output.hpp"

using namespace naiad;

std::string get_stub(const std::string & fname)
{
  const auto last{fname.find_last_of('.')};
  return fname.substr(0,last+1); // includes period
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
  const std::string fname_out{fname_stub + "out"};

  naiad::out.link_stream(std::cout);
  naiad::out.link_filename(fname_out);

  naiad::out << "BEGIN naiad ναϊάς" << std::endl;
  naiad::out << std::endl;

  // print arguments to terminal
  naiad::out << "=== COMMAND-LINE ARGUMENTS ===" << std::endl;
  for (const auto & x : args)
    naiad::out << x << std::endl;
  naiad::out << std::endl;

  naiad::out << "END naiad ναϊάς" << std::endl;

  return 0;
}
