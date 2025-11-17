#ifndef NAIAD_OUTPUT_HPP
#define NAIAD_OUTPUT_HPP

#include <string>

// Composition_stream
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

namespace naiad
{

class Composition_stream : public std::ostream
{
  public:

    Composition_stream() : std::ostream(NULL) { std::ostream::rdbuf(&my_buffer); }
    void link_stream(std::ostream & out)
    {
      out.flush();
      my_buffer.add_buffer(out.rdbuf());
    }
    void link_filename(const std::string & fname)
    {
      out_storage.emplace_back(fname);
      link_stream(out_storage.back());
    }

  private:

    class Composition_buffer : public std::streambuf
    {
      public:

        void add_buffer(std::streambuf * buf) { bufs.emplace_back(buf); }

        virtual int overflow(int c) override
        {
          for (auto & buf : bufs)
            buf->sputc(c);
          return c;
        }

      private:

        std::vector<std::streambuf *> bufs;
    };

    Composition_buffer my_buffer;
    std::vector<std::ofstream> out_storage;
};

extern Composition_stream out;

} // namespace naiad

#endif
