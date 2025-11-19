#ifndef NAIAD_DENSE_MATRIX_HPP
#define NAIAD_DENSE_MATRIX_HPP

#include <iostream>
#include <vector>

namespace naiad
{

template <typename T>
class Dense_matrix
{
  public:

    using size_type = std::size_t;
    using value_type = T;

    Dense_matrix() : m{0}, n{0} {}

    // square matrix
    Dense_matrix(size_type n_) : m{n_}, n{n_} { dat.resize(m * n); }

    // rectangular matrix
    Dense_matrix(size_type m_, size_type n_) : m{m_}, n{n_} { dat.resize(m * n); }

    Dense_matrix(size_type m_, size_type n_, const std::vector<T> & x) : m{m_}, n{n_}, dat{x} {}
    Dense_matrix(size_type m_, size_type n_, const std::vector<T> && x) : m{m_}, n{n_}, dat{std::move(x)} {}

    auto M() const { return m; }
    auto N() const { return n; }

    // necessary becasue of std::vector<bool> case
    std::vector<T>::reference operator()(size_type m_, size_type n_)
    {
      return dat[m_ + n_ * m];
    }

    const T & operator()(size_type m_, size_type n_) const
    {
      return dat[m_ + n_ * m];
    }

    Dense_matrix<T> operator=(T x)
    {
      for (auto & z : dat)
        z = x;
      return *this;
    }

    Dense_matrix<T> operator*=(T x)
    {
      for (auto & z : dat)
        z *= x;
      return *this;
    }

    Dense_matrix<T> operator/=(T x)
    {
      for (auto & z : dat)
        z /= x;
      return *this;
    }

    // take a copy so we can modify and return it
    friend Dense_matrix operator*(T x, Dense_matrix<T> z)
    {
      z *= x;
      return z;
    }

    // take a copy so we can modify and return it
    friend Dense_matrix operator/(T x, Dense_matrix<T> z)
    {
      z /= x;
      return z;
    }

    T * data() { return dat.data(); }
    const T * data() const { return dat.data(); }

    auto begin() { return dat.begin(); }
    auto end() { return dat.end(); }

    auto begin() const { return dat.begin(); }
    auto end() const { return dat.end(); }

    void print(std::ostream & os) const
    {
      for (size_type i = 0; i < n; ++i)
      {
        for (size_type j = 0; j < (m - 1); ++j)
          os << (*this)(j, i) << " , ";
        os << (*this)(m - 1, i) << "\n";
      }
    }

    void resize(const size_type m_, const size_type n_)
    {
      m = m_;
      n = n_;
      dat.resize(m * n);
    }

    void resize(const size_type square) { resize(square, square); }

  private:

    size_type m;
    size_type n;
    std::vector<T> dat;
};

} // namespace naiad

#endif
