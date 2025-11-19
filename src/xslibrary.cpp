#include "xslibrary.hpp"

#include <fstream>

#include <iostream>

namespace naiad
{

XSLibrary::XSLibrary(const std::string & filename_)
  : fname{filename_}
{
  std::ifstream ifs{fname};

  while (ifs.good())
  {

    std::string card;
    ifs >> card;
    if (card[0] == '#')
    {
      for (char x{'a'}; x != '\n'; ifs.get(x)) {}
      continue;
    }

    if (card == "")
      break;

    if (card == "ngroup")
    {
      ifs >> ng;
    }
    else if (card == "niso")
    {
      // this is not necessary, but I'll use it if I have it
      int niso;
      ifs >> niso;
      mat.reserve(niso);
    }
    else if (card == "nmoment")
    {
      ifs >> nmom;
    }
    else if (card == "name")
    {
      std::string name;
      ifs >> name;
      mat.emplace_back(name);
      mat.back().ngroup(ngroup());
      matid.insert({name, static_cast<int>(mat.size()-std::size_t{1})});
    }
    else if (card == "total")
    {
      mat.back().sigma_t.resize(ngroup());
      for (double & x : mat.back().sigma_t)
        ifs >> x;
    }
    else if (card == "scatter")
    {
      mat.back().scatter.resize(nmoment());
      for (Dense_matrix<double> & x : mat.back().scatter)
        x.resize(ngroup());
      int mom;
      ifs >> mom;
      for (int g = 0; g < ngroup(); ++g)
        for (int gprime = 0; gprime < ngroup(); ++gprime)
          ifs >> mat.back().scatter[mom](gprime,g);
    }
    else if (card == "nusf")
    {
      mat.back().isfis = true;
      mat.back().nusf.resize(ngroup());
      for (double & x : mat.back().nusf)
        ifs >> x;
    }
    else if (card == "chi")
    {
      mat.back().isfis = true;
      mat.back().chi.resize(ngroup());
      for (double & x : mat.back().chi)
        ifs >> x;
    }
    else
    {
      std::cout << "UNKNOWN XSLIBRARY CARD: " << card << std::endl;
    }


  }
}

const XSMaterial & XSLibrary::operator()(const std::string & name) const
{
  const auto find{matid.find(name)};
  if (find != matid.end())
    return (*this)(find->second);
  std::cerr << "Could not find material in library: " << name << std::endl;
  std::abort();
  return mat.front();
}

XSMaterial & XSLibrary::operator()(const std::string & name)
{
  const auto find{matid.find(name)};
  if (find != matid.end())
    return (*this)(find->second);
  std::cerr << "Could not find material in library: " << name << std::endl;
  std::abort();
  return mat.front();
}

int XSMaterial::ngroup(int n)
{
  ng = n;
  return ngroup();
}

} // namespace naiad
