#ifndef NAIAD_OUTPUT_HPP
#define NAIAD_OUTPUT_HPP

#include <string>

namespace naiad
{

class Output
{
  public:
    Output() {}

    void associate_output_filename(const std::string & fname) {}
};

extern Output out;

} // namespace naiad

#endif
