#ifndef NAIAD_EXCEPTION_HANDLER_HPP
#define NAIAD_EXCEPTION_HANDLER_HPP

#include <string>
#include <iostream>
#include <tuple>

#include "output.hpp"

namespace naiad
{

enum class Exception_level
{
  note,
  warning,
  fatal,
};

class Exception_handler
{
  public:
    Exception_handler(){}

    void note(const std::string & msg);
    void warning(const std::string & msg);
    void fatal(const std::string & msg);

    void summary(std::ostream & os = naiad::out) const;

  private:

    static void print_msg(const std::tuple<Exception_level,std::string> & msg, std::ostream & os);

    std::vector<std::tuple<Exception_level,std::string>> arr;

};

extern Exception_handler exception;

} // namespace naiad

#endif
