#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>

template<class T, class A=std::allocator<T> >
struct EventMatrix {
  typedef T value_type;
  typedef std::vector<value_type, A> Container;

  EventMatrix() : _b(0) {}
  EventMatrix(int a, int b, value_type const& initial=value_type())
  : _b(0)
  {
    resize(a, b, initial);
  }
  EventMatrix(EventMatrix const& other)
  : _data(other._data), _b(other._b)
  {}

  EventMatrix& operator=(EventMatrix copy) {
    swap(*this, copy);
    return *this;
  }

  bool empty() const { return _data.empty(); }
  void clear() { _data.clear(); _b = 0; }

  int dim_a() const { return _b ? _data.size() / _b : 0; }
  int dim_b() const { return _b; }

  value_type* operator[](int a) {
    return &_data[a * _b];
  }
  value_type const* operator[](int a) const {
    return &_data[a * _b];
  }

  void resize(int a, int b, value_type const& initial=value_type()) {
    if (a == 0) {	
      b = 0;
    }
    //std::cout << a*b << " ";
    _data.resize(a * b, initial);
    //_mostrecent.resize(1, 0);
    _b = b;
    //std::cout << b << " ";
    //std::cout << _b << " " ;
    //std::cout << _data.size() << std::endl;
  }
  
  void setMostRecent(value_type const& recent=value_type()) {
	_mostrecent[0] = recent;  
  }

  value_type getMostRecent() {
	return _mostrecent[0];  
  }

  friend void swap(EventMatrix& a, EventMatrix& b) {
    using std::swap;
    swap(a._data, b._data);
    swap(a._b,    b._b);
  }

  template<class Stream>
  friend Stream& operator<<(Stream& s, EventMatrix const& value) {
    s << "<Matrix at " << &value << " dimensions "
      << value.dim_a() << 'x' << value.dim_b();
    if (!value.empty()) {
      bool first = true;
      typedef typename Container::const_iterator Iter;
      Iter i = value._data.begin(), end = value._data.end();
      while (i != end) {
        s << (first ? " [[" : "], [");
        first = false;
        s << *i;
        ++i;
        for (int b = value._b - 1; b; --b) {
          s << ", " << *i;
          ++i;
        }
      }
      s << "]]";
    }
    s << '>';
    return s;
  }

private:
  Container _data;
  int _b;
  Container _mostrecent;
};
