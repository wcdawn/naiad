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
    std::cout << "Expect exactly one argument (filename)." << std::endl;
    return 1;
  }

  const std::string fname_inp{args[1]};
  const std::string fname_stub{get_stub(fname_inp)};
  const std::string fname_out{fname_stub + "out"};

  naiad::out.associate_output_filename(fname_out);

  std::cout << "BEGIN naiad ναϊάς" << std::endl;
  std::cout << std::endl;

  // print arguments to terminal
  std::cout << "=== COMMAND-LINE ARGUMENTS ===" << std::endl;
  for (const auto & x : args)
    std::cout << x << std::endl;
  std::cout << std::endl;

  std::cout << "END naiad ναϊάς" << std::endl;

  return 0;
}
