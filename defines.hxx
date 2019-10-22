#pragma once
#include <array>
#include <vector>
#include <set>

namespace MEXT_NS {
  template<class T, std::size_t N>
  class Tuple0 : public std::array<T, N>
  {
  public:
    enum
    {
      DIM = N
    };
    Tuple0()
    {
      memset((void*)this->data(), 0, N * sizeof(T));
    }

    Tuple0(const T* da)
    {
      memcpy(this->data(), da, N * sizeof(T));
    }

    Tuple0(const T& val)
    {
      for (size_t i = 0; i < N; i++)
        (*this)[i] = val;
    }

    Tuple0(const Tuple0<T, N>& a, const Tuple0<T, N>& b)
    {
      for (size_t i = 0; i < N; i++)
        (*this)[i] = b[i] - a[i];
    }

    int find(const T& val)const
    {
      for (size_t i = 0; i < N; i++)
      {
        if ((*this)[i] == val)
          return (int)i;
      }
      return -1;
    }

    Tuple0<T, N>& operator*=(const T& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] *= val;
      }
      return *this;
    }

    Tuple0<T, N>& operator/=(const T& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] /= val;
      }
      return *this;
    }

    Tuple0<T, N>& operator+=(const T& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] += val;
      }
      return *this;
    }

    Tuple0<T, N>& operator+=(const Tuple0<T, N>& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] += val[i];
      }
      return *this;
    }

    Tuple0<T, N>& operator-=(const Tuple0<T, N>& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] -= val[i];
      }
      return *this;
    }

    Tuple0<T, N>& operator-=(const T& val)
    {
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] -= val;
      }
      return *this;
    }

    bool operator<(const Tuple0<T, N>& val)const
    {
      for (size_t i = 0; i < N; i++)
      {
        if ((*this)[i] < val[i])
          return true;
        if (val[i] < (*this)[i])
          return false;
      }
      return false;
    }

    bool operator==(const Tuple0<T, N>& val)const
    {
      for (size_t i = 0; i < N; i++)
      {
        if ((*this)[i] != val[i])
          return false;
      }
      return true;
    }

    Tuple0<T, N>& set(const T* da)
    {
      memcpy(this->elems, da, N * sizeof(T));
      return *this;
    }

    Tuple0<T, N>& operator=(const Tuple0<T, N>&val)
    {
      if (this != &val)
      {
        for (size_t i = 0; i < N; i++)
          (*this)[i] = val[i];
      }
      return *this;
    }

    T lengthSqr()const
    {
      T rt(0);
      for (size_t i = 0; i < N; i++)
      {
        rt += ((*this)[i] * (*this)[i]);
      }
      return rt;
    }

    T distanceSqr(const Tuple0<T, N>& rhs)const
    {
      T rt(0);
      for (size_t i = 0; i < N; i++)
      {
        T diff2 = (*this)[i] - rhs[i];
        diff2 *= diff2;
        rt += diff2;
      }
      return rt;
    }

    T getMax()const
    {
      T rt((*this)[0]);
      for (size_t i = 1; i < N; i++)
      {
        rt = std::max(rt, (*this)[i]);
      }
      return rt;
    }

    T dot(const Tuple0<T, N>& rhs)const
    {
      T rt(0);
      for (size_t i = 0; i < N; i++)
      {
        rt += (*this)[i] * rhs[i];
      }
      return rt;
    }

    T norm()const
    {
      return std::sqrt(lengthSqr());
    }

    void normalize()
    {
      T len = std::sqrt(lengthSqr());
      for (size_t i = 0; i < N; i++)
      {
        (*this)[i] /= len;
      }
    }

  };


  template <class T>
  class Couple : public Tuple0<T, 2>
  {
  public:
    Couple() {}

    Couple(const T* data)
      :Tuple0<T, 2>(data)
    {}

    Couple(const T& a, const T& b)
    {
      (*this)[0] = a;
      (*this)[1] = b;
    }

    Couple<T>&   set(const T& a, const T& b)
    {
      (*this)[0] = a;
      (*this)[1] = b;
      return *this;
    }

    void          reverse()
    {
      std::swap((*this)[0], (*this)[1]);
    }

    int find(const T& val)const
    {
      for (int i = 0; i<this->size(); i++)
      {
        if (val == (*this)[i])
          return i;
      }
      return -1;
    }

    Couple<T>& operator=(const Couple<T>&val)
    {
      Tuple0<T, 2>::operator=(val);
      return *this;
    }
  };


  template <class T>
  class Triple : public Tuple0<T, 3>
  {
  public:
    Triple() {}

    Triple(const T& a, const T& b, const T& c)
    {
      (*this)[0] = a;
      (*this)[1] = b;
      (*this)[2] = c;
    }

    Triple(const T* data)
      :Tuple0<T, 3>(data)
    { }


    Triple(const Triple<T>& a, const Triple<T>& b)
    {
      set(a, b);
    }

    Triple<T>&   set(const T& a, const T& b, const T& c)
    {
      (*this)[0] = a;
      (*this)[1] = b;
      (*this)[2] = c;
      return *this;
    }


    void  set(const Triple<T>& a, const Triple<T>& b)
    {
      this->set(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    }

    Triple<T>& set(const T* da)
    {
      this->set(da[0], da[1], da[2]);
      return *this;
    }

    void  reverse()
    {
      std::swap((*this)[0], (*this)[1]);
    }

    void sort()
    {
      Sort3<T>(this->elems[0], this->elems[1], this->elems[2]);
    }

    T  dot(const Triple<T>& vec)const
    {
      return (*this)[0] * vec[0] + (*this)[1] * vec[1] + (*this)[2] * vec[2];
    }

    Triple<T> cross(const Triple<T>& vec)const
    {
      return Triple<T>((*this)[1] * vec[2] - (*this)[2] * vec[1],
        (*this)[2] * vec[0] - (*this)[0] * vec[2],
        (*this)[0] * vec[1] - (*this)[1] * vec[0]);
    }

    T len2()const
    {
      return (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] + (*this)[2] * (*this)[2];
    }
  };


  template<class T>
  inline Triple<T> operator *(const Triple<T>& a, T val) { return Triple<T>(a[0] * val, a[1] * val, a[2] * val); }

  template<class T>
  /**
  * @brief The KeyIndex3 class
  */
  class  KeyIndex3 : public  Triple<T>
  {
  public:
    KeyIndex3(const T& k0 = T(), const T& k1 = T(), const T& k2 = T())
    {
      set(k0, k1, k2);
    }

    KeyIndex3(const T* k)
    {
      set(k[0], k[1], k[2]);
    }

    KeyIndex3(const  Triple<T>& k)
    {
      set(k[0], k[1], k[2]);
    }

    void  set(const T& k0, const T& k1, const T& k2)
    {
      if (k0<k1)
      {
        if (k1<k2)
        {
          (*this)[0] = k0;
          (*this)[1] = k1;
          (*this)[2] = k2;
        }
        else
        {
          (*this)[2] = k1;
          if (k0<k2)
          {
            (*this)[0] = k0;
            (*this)[1] = k2;
          }
          else
          {
            (*this)[0] = k2;
            (*this)[1] = k0;
          }
        }
      }
      else
      {
        if (k2<k1)
        {
          (*this)[0] = k2;
          (*this)[1] = k1;
          (*this)[2] = k0;
        }
        else
        {
          (*this)[0] = k1;
          if (k2<k0)
          {
            (*this)[1] = k2;
            (*this)[2] = k0;
          }
          else
          {
            (*this)[1] = k0;
            (*this)[2] = k2;
          }
        }

      }
    }

    bool operator<(const KeyIndex3<T>& src)const
    {
      if ((*this)[0]<src[0]) return true;
      if (src[0]<(*this)[0]) return false;
      if ((*this)[1]<src[1]) return true;
      if (src[1]<(*this)[1]) return false;
      return (*this)[2]<src[2];
    }

    template<typename T1>
    bool  operator==(const KeyIndex3<T1>& src)const
    {
      if ((*this)[0] == src[0])
      {
        if ((*this)[1] == src[1])
        {
          return  (*this)[2] == src[2];
        }
      }
      return false;
    }

    bool operator!=(const KeyIndex3<T>& src)const
    {
      return !(*this == src);
    }
  };

  template<class POINT, unsigned DIM>
  class tBBox
  {
  public:
    POINT  pmin_max[2];
    typedef POINT value_type;
    typedef typename POINT::value_type FType;
    tBBox()
    {
      reset();
    }

    tBBox(const FType* coord)
    {
      pmin_max[0].set(coord);
      pmin_max[1].set(coord + DIM);
    }

    tBBox(const FType* coord0, const FType* coord1)
    {
      *this |= coord0;
      *this |= coord1;
    }

    POINT& min() { return pmin_max[0]; }
    POINT& max() { return pmin_max[1]; }

    const POINT& min()const { return pmin_max[0]; }
    const POINT& max()const { return pmin_max[1]; }

    void reset(FType imax = std::numeric_limits<FType>::max())
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        pmin_max[0][i] = imax;
        pmin_max[1][i] = -imax;
      }
    }


    tBBox<POINT, DIM>& operator |=(const FType* pt)
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (pt[i]<pmin_max[0][i]) pmin_max[0][i] = pt[i];
        if (pt[i]>pmin_max[1][i]) pmin_max[1][i] = pt[i];
      }
      return *this;
    }

    template<class POINT2>
    tBBox<POINT, DIM>& operator |=(const POINT2& pt)
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (pt[i]<pmin_max[0][i]) pmin_max[0][i] = pt[i];
        if (pt[i]>pmin_max[1][i]) pmin_max[1][i] = pt[i];
      }
      return *this;
    }

    tBBox<POINT, DIM>& operator +=(const FType* pt)
    {
      *this |= pt;
      return *this;
    }

    tBBox<POINT, DIM>& extend(const FType* pt)
    {
      *this |= pt;
      return *this;
    }

    tBBox<POINT, DIM>& extendBox(const FType* pt)
    {
      *this |= pt;
      *this |= pt + DIM;
      return *this;
    }

    tBBox<POINT, DIM>& operator |=(const tBBox<POINT, DIM>& abox)
    {
      *this |= abox.pmin_max[0];
      *this |= abox.pmin_max[1];
      return *this;
    }

    tBBox<POINT, DIM>& operator +=(const tBBox<POINT, DIM>& abox)
    {
      *this |= abox;
      return *this;
    }

    tBBox<POINT, DIM>& operator |=(const POINT& pt)
    {
      return (*this) |= pt.data();
    }

    tBBox<POINT, DIM>& operator +=(const POINT& pt)
    {
      return (*this) |= pt.data();
    }

    tBBox<POINT, DIM>& operator +=(FType thickness)
    {
      assert(thickness>0);
      for (unsigned i = 0; i<DIM; i++)
      {
        pmin_max[0][i] -= thickness;
        pmin_max[1][i] += thickness;
      }
      return *this;
    }

    tBBox<POINT, DIM>&   scale(const FType* s)
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        FType dx = pmin_max[1][i] - pmin_max[0][i];
        FType exx = (dx*(s[i] - 1)) / 2;
        pmin_max[0][i] -= exx;
        pmin_max[1][i] += exx;
      }
      return *this;
    }

    tBBox<POINT, DIM>&   scale(const FType& s)
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        FType dx = pmin_max[1][i] - pmin_max[0][i];
        FType exx = (dx*(s - 1.0)) / 2;
        pmin_max[0][i] -= exx;
        pmin_max[1][i] += exx;
      }
      return *this;
    }

    tBBox<POINT, DIM>& operator *=(FType sc)
    {
      FType s[DIM];
      for (unsigned i = 0; i<DIM; i++)
        s[i] = sc;

      this->scale(s);
      return *this;
    }

    void out_sphere_it(FType sc = 1)
    {
      POINT cen = pmin_max[0] + pmin_max[1];
      cen *= FType(0.5);
      FType offset = sqrt(pmin_max[0].distanceSqr(pmin_max[1]))*sc*FType(0.5);
      for (unsigned i = 0; i<DIM; i++)
      {
        pmin_max[0][i] = cen[i] - offset;
        pmin_max[1][i] = cen[i] + offset;
      }
    }

    void thicknessIt(FType thick)
    {
      assert(thick>0);
      for (unsigned i = 0; i<DIM; i++)
      {
        pmin_max[0][i] -= thick;
        pmin_max[1][i] += thick;
      }
    }

    bool overlaps(const tBBox<POINT, DIM>& rhs)const
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (this->pmin_max[0][i] > rhs.pmin_max[1][i] ||
          rhs.pmin_max[0][i] > this->pmin_max[1][i])
        {
          return false;
        }
      }
      return true;
    }

    void move_center(const POINT& cen)
    {
      POINT dir(pmin_max[0], pmin_max[1]);
      dir *= FType(0.5);
      pmin_max[0].set(cen);
      pmin_max[1].set(cen);
      pmin_max[0] -= dir;
      pmin_max[1] += dir;
    }

    void move_in(POINT& cen)
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (cen[i]>pmin_max[1][i])
          cen[i] = pmin_max[1][i];
        if (cen[i]<pmin_max[0][i])
          cen[i] = pmin_max[0][i];
      }
    }

    FType size(bool avg = false)const
    {
      if (avg)
      {
        double tol = this->size(false)*1.e-6;
        auto ext = extend();
        FType mass = 0.0;
        int dim = 0;
        for (unsigned i = 0; i < DIM; i++)
        {
          if (ext[i] > tol) {
            dim += 1;
            if (mass > tol)
            {
              mass *= ext[i];
            }
            else
            {
              mass = ext[i];
            }
          }

        }
        return std::pow(mass, FType(1) / dim);
      }
      return sqrt(pmin_max[0].distanceSqr(pmin_max[1]));
    }

    bool isValid()const
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (pmin_max[0][i]>pmin_max[1][i])
          return false;
      }
      return true;
    }

    bool isIn(const POINT& pt, FType eps = 0) const
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (pt[i] >(pmin_max[1][i] + eps) || (pt[i] + eps) < pmin_max[0][i])
          return false;
      }
      return true;
    }

    bool  isDisjoint(const tBBox<POINT, DIM>& other, FType eps = 0)const
    {
      for (unsigned i = 0; i<DIM; i++)
      {
        if (other.pmin_max[0][i] >(pmin_max[1][i] + eps) || (other.pmin_max[1][i] + eps) < pmin_max[0][i])
          return true;
      }
      return false;
    }

    POINT extend()const { return (pmin_max[1] - pmin_max[0]); }
    POINT center()const { return (pmin_max[1] + pmin_max[0])*0.5; }
    FType mass()const
    {
      POINT ext = this->extend();
      FType rt = ext[0];
      for (unsigned i = 1; i<DIM; i++)
        rt *= ext[i];
      return rt;
    }

    void validIt()
    {
      auto ext = this->extend();
      auto minLen = ext[0], maxLen = ext[0];
      int imin = 0, imax = 0;
      for (int i = 1; i < DIM; i++) {
        if (ext[i] < minLen) {
          imin = i;
          minLen = ext[i];
        }

        if (maxLen<ext[i]) {
          imax = i;
          maxLen = ext[i];
        }
      }
      auto eps = maxLen * 0.01;
      if (minLen < eps) {
        POINT cen = (pmin_max[0] + pmin_max[1])*FType(0.5);
        FType offset = eps;
        for (int i = 0; i < DIM; i++) {
          if (ext[i] < eps) {
            pmin_max[0][i] = cen[i] - offset;
            pmin_max[1][i] = cen[i] + offset;
          }
        }
      }
    }
  };


  template<class T, std::size_t N>
  Tuple0<int32_t, N> floor(const Tuple0<T, N>& val)
  {
    Tuple0<int32_t, N> rt;
    for (uint32_t i = 0; i < N; i++)
    {
      rt[i] = std::floor(val[i]);
    }
    return rt;
  }


  template<class T, std::size_t N>
  Tuple0<T, N> operator-(const Tuple0<T, N>& val1, const Tuple0<T, N>& val0)
  {
    Tuple0<T, N> rt;
    for (std::size_t i = 0; i < N; i++)
      rt[i] = val1[i] - val0[i];
    return rt;
  }

  template<class T, std::size_t N>
  Tuple0<T, N> operator-(const Tuple0<T, N>& val1, T val0)
  {
    Tuple0<T, N> rt;
    for (std::size_t i = 0; i < N; i++)
      rt[i] = val1[i] - val0;
    return rt;
  }

  template<class T, std::size_t N>
  Tuple0<T, N> operator+(const Tuple0<T, N>& val1, T val0)
  {
    Tuple0<T, N> rt;
    for (std::size_t i = 0; i < N; i++)
      rt[i] = val1[i] + val0;
    return rt;
  }

  template<class T, std::size_t N>
  Tuple0<T, N> operator+(const Tuple0<T, N>& val1, const Tuple0<T, N>& val0)
  {
    Tuple0<T, N> rt;
    for (std::size_t i = 0; i < N; i++)
      rt[i] = val1[i] + val0[i];
    return rt;
  }

  typedef size_t        HashValue;
  typedef Triple<int>   Int3;

  inline HashValue Hash(HashValue h11, HashValue h21)
  {
    HashValue h1(h11);
    HashValue h2(h21);
    return HashValue(((h1 << 16) | (h1 >> 16)) ^ h2);
  }

  inline HashValue Hash(const KeyIndex3<int>& val) { return Hash(Hash(val[0], val[1]), val[2]); }

  inline HashValue Hash(const Int3& val)
  {
    return Hash(KeyIndex3<int>(val[0], val[1], val[2]));
  }

  struct HashInt3
  {
    HashValue operator()(const Int3& val) const { return Hash(val); }
    bool operator()(const Int3& val0, const Int3& val1)const
    {
      KeyIndex3<int> k0(val0[0], val0[1], val0[2]);
      KeyIndex3<int>  k1(val1[0], val1[1], val1[2]);
      return k0 == k1;
    }
  };

  typedef Couple<int>          Int2;
  typedef std::vector<Int2>    DInt2Array;
  typedef std::vector<int>     DIntArray;
  typedef std::set<int>        IntSet;

  typedef boost::unordered_set<Int3, HashInt3, HashInt3> Int3Set;

};


