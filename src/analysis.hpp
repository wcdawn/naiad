#ifndef NAIAD_ANALYSIS_HPP
#define NAIAD_ANALYSIS_HPP

#include "geometry.hpp"
#include "result.hpp"
#include "xslibrary.hpp"

#include <vector>
#include <iostream>
#include <string>

namespace naiad
{

enum class Analysis_reference
{
  none,
  onegroup,
  twogroup,
  tworegion,
  critical,
};

std::string enum2str(const Analysis_reference & analysis_ref);
Analysis_reference str2enum_analysis_reference(const std::string & str);

class Analysis
{
  public:
    Analysis(const Geometry & geo_, const XSLibrary & xslib_, const Result & res_)
      : geo{geo_}, xslib{xslib_}, res{res_}
    {}

    void summary(std::ostream & os) const;

    virtual ~Analysis() {}

  protected:
    const Geometry & geo;
    const XSLibrary & xslib;
    const Result & res;

    std::vector<double> error_linf() const;
    double error_keff() const;

    virtual std::vector<std::vector<double>> exact_flux() const = 0;
    virtual double exact_keff() const = 0;

};

class Analysis_critical : public Analysis
{
  public:
    using Analysis::Analysis;
  protected:
    std::vector<std::vector<double>> exact_flux() const override;
    double exact_keff() const override;
};

class Analysis_onegroup : public Analysis
{
  public:
    using Analysis::Analysis;
  protected:
    std::vector<std::vector<double>> exact_flux() const override;
    double exact_keff() const override;
};

} // namespace naiad

#endif
