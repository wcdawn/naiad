#include "input.hpp"

#include <fstream>
#include <string>

#include <iostream>

std::string get_card(std::ifstream & ifs)
{
  std::string card;
  ifs >> card;
  return card;
}

namespace naiad
{

  Input::Input(const std::string & filename)
  {

    std::ifstream ifs{filename};

    std::string fname_xslib;

    while (ifs.good())
    {

      const std::string card{get_card(ifs)};
      if (card[0] == '#')
      {
        for (char x{'a'}; x != '\n'; ifs.get(x)) {}
      }

      if (card == "")
        break;

      if (card == "xslib")
      {
        ifs >> fname_xslib;
      }

    }

  }

} // namespace naiad
