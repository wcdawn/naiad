#include "exception_handler.hpp"

namespace naiad
{

Exception_handler exception;

void Exception_handler::note(const std::string & msg)
{
  arr.emplace_back(std::tuple{Exception_level::note, msg});
  print_msg(arr.back(), naiad::out);
}

void Exception_handler::warning(const std::string & msg)
{
  arr.emplace_back(std::tuple{Exception_level::warning, msg});
  print_msg(arr.back(), naiad::out);
}

void Exception_handler::fatal(const std::string & msg)
{
  arr.emplace_back(std::tuple{Exception_level::fatal, msg});
  print_msg(arr.back(), naiad::out);
  summary();
  std::abort();
}

void Exception_handler::print_msg(const std::tuple<Exception_level,std::string> & msg, std::ostream & os)
{
  switch (std::get<0>(msg))
  {
    case (Exception_level::note):
      os << "NOTE :: ";
      break;
    case (Exception_level::warning):
      os << "WARNING :: ";
      break;
    case (Exception_level::fatal):
      os << "FATAL :: ";
      break;
    default:
      std::cerr << "Failed to identify exception level." << std::endl;
      std::abort();
  }
  os << std::get<1>(msg) << std::endl;
}

void Exception_handler::summary(std::ostream & os) const
{
  os << "=== EXCEPTION SUMMARY ===" << std::endl;

  std::unordered_map<Exception_level,int> count;

  count.insert({Exception_level::note, 0});
  count.insert({Exception_level::warning, 0});
  count.insert({Exception_level::fatal, 0});

  for (const auto & e : arr)
  {
    count.at(std::get<0>(e))++;
    print_msg(e, os);
  }
  os << count.at(Exception_level::note) << " note(s)" << std::endl;
  os << count.at(Exception_level::warning) << " warning(s)" << std::endl;
  os << count.at(Exception_level::fatal) << " fatal(s)" << std::endl;
  os << std::endl;
}

} // namespace naiad
