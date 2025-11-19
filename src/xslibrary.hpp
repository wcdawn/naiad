#ifndef NAIAD_XSLIBRARY_HPP
#define NAIAD_XSLIBRARY_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>

#include "dense_matrix.hpp"

namespace naiad
{

class XSMaterial
{
  public:
    XSMaterial(){}
    XSMaterial(const std::string & n) : nam{n} {}

    int ngroup() const { return ng; }
    int ngroup(int n);
    std::string name() const { return nam; }

    void summary(std::ostream & os) const;

    std::vector<double> sigma_t;
    std::vector<Dense_matrix<double>> scatter;

    std::vector<double> nusf;
    std::vector<double> sigma_f;
    std::vector<double> chi;

    bool isfis{false};

  private:
    int ng;
    std::string nam;
};

class XSLibrary
{
  public:
    XSLibrary(){}
    XSLibrary(const std::string & filename_);

    const XSMaterial & operator()(int i) const { return mat[i]; }
    XSMaterial & operator()(int i) { return mat[i]; }

    const XSMaterial & operator()(const std::string & name) const;
    XSMaterial & operator()(const std::string & name);

    void summary(std::ostream & os) const;

    std::string filename() const { return fname; }
    int ngroup() const { return ng; }
    int nmoment() const { return nmom; }

  private:
    std::string fname;

    std::vector<XSMaterial> mat;
    std::unordered_map<std::string, int> matid;

    int ng;
    int nmom;

};

} // namespace naiad

#endif
