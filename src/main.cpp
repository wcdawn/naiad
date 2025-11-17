#include <iostream>
#include <vector>
#include <string>

int main(int argc, char* argv[])
{

  std::cout << "BEGIN naiad ναϊάς" << std::endl;
  std::cout << std::endl;

  // process arguments
  std::vector<std::string> args;
  args.reserve(argc);
  for (int i = 0; i < argc; ++i)
    args.emplace_back(argv[i]);

  // print arguments to terminal
  std::cout << "=== COMMAND-LINE ARGUMENTS ===" << std::endl;
  for (const auto & x : args)
    std::cout << x << std::endl;
  std::cout << std::endl;

  std::cout << "END naiad ναϊάς" << std::endl;

  return 0;
}
