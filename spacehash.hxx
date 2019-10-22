#pragma once

#include <functional>
#include <vector>
#include <omp.h>

#include "defines.hxx"

namespace MEXT_NS {
  /**
   * @brief The SpaceHash class
   * space hash for collision test
   */

  template <typename T, uint32_t N> class tSpaceHash
  {
  public:
    typedef T ValueType;
    typedef Tuple0<ValueType, N> Point;
    typedef tBBox<Point, N>      Box;
    typedef std::function<void(size_t, ValueType *)> PointFunction;
    typedef std::function<void(size_t, int32_t *)> TriangleFunction;
    typedef std::function<bool(size_t, Box &)> BoxFunction;

    typedef Tuple0<int32_t, N> IntN;

    /**
     * @brief create
     * create the "SpaceHash" with n points
     * @param nPoint, number of points
     * @param pointFun, call back of getting the i's point
     * @param resolution, resolution of the SpaceHash
     */
    bool create(size_t nPoint, PointFunction pointFun, uint32_t resolution = 64, ValueType enlageFactor = 1.1, ValueType df = 0);
    bool create(const Point &min, const Point &max, uint32_t resolution = 64, ValueType enlageFactor = 1.1, ValueType df = 0);

    IntN getCellIndex(const ValueType *) const;
    IntN getCellIndex(const Point &p) const
    {
      return this->getCellIndex(p.data());
    }

    size_t cellNum() const;
    void getCellPointMap(size_t nPoint, PointFunction, std::vector<int32_t> &, std::vector<int32_t> &) const;

    size_t index(const IntN &) const;

    const Point &min() const;
    Point max() const;
    ValueType step()const { return mStep; }
    const IntN& size()const { return mSize; }
    /**
     * @brief getCollisionCellPairs
     * get the box pairs which are overlap each other
     * @param nBox, number of box
     * @param boxFun, call back of getting the i's box
     * @param pairCB, call back of returned box index pair
     */
    void getCollisionCellPairs(size_t nBox, BoxFunction boxFun, std::function<void(const Int2 &)> pairCB) const;

    void getCollisionCellPairs(size_t nBox, BoxFunction boxFun, std::function<bool(const Int2 &)>, DInt2Array&) const;


    void getCollisionCellPairs(size_t nBox0,
      BoxFunction boxFun0,
      size_t nBox1,
      BoxFunction boxFun1,
      std::function<void(const Int2 &)> pairCB) const;

    bool isInBBox(const ValueType *)const;

  protected:
    void putBox2Cells(size_t nBox, BoxFunction boxFun, DInt2Array &cellBoxPair) const;
    void putBox2Cells(size_t nBox, BoxFunction boxFun, std::vector<std::vector<size_t>> &) const;
    void getCollisionCellPairs(const DInt2Array& range, const DInt2Array& pointCellIndex, BoxFunction boxFun, std::function<void(const Int2 &)> pairCB) const;
    void getCollisionCellPairs(const DInt2Array& range, const DInt2Array& pointCellIndex, BoxFunction boxFun, std::function<bool(const Int2 &)>, DInt2Array&) const;

  protected:
    IntN mSize;
    Point mStart;
    ValueType mStep, mStepInv, mEnlageFactor;
    int32_t mResolution;
  };

  template <typename T, uint32_t N>
  /**
   * @brief The tPointHash class
   */
  class tPointHash : public tSpaceHash<T, N>
  {
  public:
    typedef typename tSpaceHash<T, N>::ValueType ValueType;
    tPointHash(size_t nPoint,
      typename tSpaceHash<T, N>::PointFunction,
      uint32_t resolution = 64,
      double enlageFactor = 1.1);
    void getIncludes(const typename tSpaceHash<T, N>::ValueType *pmin,
      const typename tSpaceHash<T, N>::ValueType *pmax,
      std::function<void(int32_t)>) const;

    /**
     * @brief getNearest
     * get the nearest point of point p
     * @param p
     * @return
     */
    int32_t getNearest(const ValueType *p, std::function<bool(size_t)>) const;

  private:
    size_t mNumPoint;
    typename tSpaceHash<T, N>::PointFunction mPointFun;
    std::vector<int32_t> mInitAddress, mCellPointMap;
  };



  template <typename T, uint32_t N>
  class tBoxHash : public tSpaceHash<T, N>
  {
  public:
    typedef typename tSpaceHash<T, N>::ValueType ValueType;
    typedef typename tSpaceHash<T, N>::BoxFunction BoxFunction;
    typedef typename tSpaceHash<T, N>::Box Box;
    typedef typename tSpaceHash<T, N>::IntN IntN;
    typedef typename tSpaceHash<T, N>::Point Point;


    tBoxHash(size_t nBox, BoxFunction boxFun,
      uint32_t resolution = 64,
      double enlageFactor = 1.1);

    void getSelfOverlap(std::function<void(const Int2 &)> pairCB) const;
    void getSelfOverlap(std::function<bool(const Int2 &)>, DInt2Array&) const;


    void getIncludes(const IntN&, std::function<bool(int32_t)>, IntSet&) const;

    void getBoxIncludes(const typename tSpaceHash<T, N>::ValueType *pmin,
      const typename tSpaceHash<T, N>::ValueType *pmax,
      std::function<bool(int32_t)>) const;

    void getIncludes(const typename tSpaceHash<T, N>::ValueType *pt,
      std::function<bool(int32_t)>) const;

    /**
    * @brief getNearest
    * get the nearest boundary boxs which includes point p
    * @param p
    * @return
    */
    int32_t getNearest(const ValueType *p, std::function<bool(size_t)>) const;
    void    getIntersected(const Point& from, const Point& to, std::function<bool(int32_t)>)const;

  private:
    BoxFunction mBoxFun;
    DInt2Array mRange, mPointCellIndex;
  };

  template <typename T, uint32_t N>
  /**
   * @brief The tPointHashUniqueInserter class
   * insert point to array and make sure there is no confilict each other wihin a tolerance
   */
  class tPointHashUniqueInserter : public tSpaceHash<T, N>
  {
  public:
    typedef typename tSpaceHash<T, N>::ValueType ValueType;
    typedef typename tSpaceHash<T, N>::Point Point;
    typedef typename tSpaceHash<T, N>::IntN IntN;

    /**
     * @brief insertPoint
     * insert point to array and make sure there is no confilict each other wihin a tolerance
     * @param point: point to be inserted
     * @param eps: tolerance
     * @return the index in point array
     */
    int32_t insertPoint(const Point &point, ValueType eps = 0);

    /**
     * @brief appendPoint: append the point without distance checking in a tolerance.
     * @param point
     * @return
     */
    int32_t appendPoint(const Point &point);


    int32_t hasPointIn(const Point &point, ValueType eps, std::function<bool(int32_t)>)const;
    int32_t nearest(const Point &point, std::function<bool(int32_t)>)const;


    const std::vector<Point>& getPointArray()const
    {
      return mPoints;
    }

    const Point& operator[](size_t i)const { return mPoints[i]; }

    void  setPointLocalIndex(size_t i, int idx);
    int   getPointLocalIndex(size_t i)const
    {
      if (i >= mLocalIndex.size())
        return -1;
      return mLocalIndex[i];
    }

  protected:
    std::vector<Point> mPoints;
    std::vector<std::vector<size_t>> mBructPointList;
    DIntArray          mLocalIndex;
  };

  template <typename T, uint32_t N>
  bool tSpaceHash<T, N>::create(size_t nPoint, PointFunction pfun, uint32_t resolution, ValueType enlageFactor, ValueType dfMin)
  {
    if (nPoint < 2 || enlageFactor <= 1.0)
      return false;
    assert(resolution >= 1.0);
    Tuple0<ValueType, N> bmin(std::numeric_limits<ValueType>::max());
    Tuple0<ValueType, N> bmax(-std::numeric_limits<ValueType>::max());
    for (size_t i = 0; i < nPoint; i++) {
      ValueType pos[N];
      pfun(i, pos);
      for (uint32_t k = 0; k < N; k++) {
        bmin[k] = std::min(bmin[k], pos[k]);
        bmax[k] = std::max(bmax[k], pos[k]);
      }
    }
    return this->create(bmin, bmax, resolution, enlageFactor, dfMin);
  }

  template <typename T, uint32_t N>
  bool
    tSpaceHash<T, N>::create(const Point &bmin, const Point &bmax,
      uint32_t resolution, ValueType enlageFactor, ValueType dfMax)
  {
    Tuple0<ValueType, N> diff(bmin, bmax);
    double lenSqrt = sqrt(diff.lengthSqr());
    if (lenSqrt < std::numeric_limits<ValueType>::epsilon())
      return false;
    mEnlageFactor = enlageFactor;
    mResolution = resolution;
    diff = Tuple0<ValueType, N>(diff.getMax());
    diff *= (ValueType)(enlageFactor - 1);

    mStart = bmin - diff;
    Tuple0<ValueType, N> end = bmax + diff;
    diff = end - mStart;
    if (N == 3) {
      mStep = std::pow(diff[0] * diff[1] * diff[2], 1.0 / 3.0) / resolution;
    }
    else if (N == 2) {
      mStep = sqrt(diff[0] * diff[1]) / resolution;
    }
    if (dfMax > std::numeric_limits<T>::epsilon())
    {
      mStep = dfMax;
    }
    mStepInv = ValueType(1) / mStep;
    Point val = (end - mStart);
    val *= mStepInv;
    for (uint32_t i = 0; i < N; i++)
    {
      mSize[i] = (IntN::value_type) std::floor(val[i]);
    }
    return true;
  }

  template <typename T, uint32_t N>
  bool tSpaceHash<T, N>::isInBBox(const ValueType * pt)const
  {
    auto pmax = this->max();
    for (uint32_t i = 0; i < N; i++) {
      if (pt[i] < mStart[i]) {
        return false;
      }

      if (pt[i] > pmax[i]) {
        return false;
      }
    }
    return true;
  }

  template <typename T, uint32_t N>
  typename tSpaceHash<T, N>::IntN
    tSpaceHash<T, N>::getCellIndex(const ValueType *pos) const
  {
    Point posI(pos);
    IntN idx;
    for (uint32_t i = 0; i < N; i++) {
      idx[i] = std::max(0, std::min(typename IntN::value_type((pos[i] - mStart[i]) * mStepInv), mSize[i] - 1));
    }
    return idx;
  }

  template <typename T, uint32_t N>
  const typename tSpaceHash<T, N>::Point &
    tSpaceHash<T, N>::min() const
  {
    return mStart;
  }

  template <typename T, uint32_t N>
  typename tSpaceHash<T, N>::Point
    tSpaceHash<T, N>::max() const
  {
    Tuple0<ValueType, N> rt;
    for (uint32_t i = 0; i < N; i++) {
      rt[i] = mStart[i] + mSize[i] * mStep;
    }
    return rt;
  }

  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCellPointMap(size_t nPoint,
      PointFunction pfun,
      std::vector<int32_t> &initAddress,
      std::vector<int32_t> &cellPointMap) const
  {
    if (nPoint < 1)
      return;
    std::vector<size_t> pointCellIndex(nPoint);
    int32_t numPt = (int32_t)nPoint;
#if defined(_USE_OMP)
#pragma omp parallel for
#endif
    for (int32_t i = 0; i < numPt; i++) {
      ValueType pos[N];
      pfun(i, pos);
      IntN idx = getCellIndex(pos);
      pointCellIndex[i] = index(idx);
    }

    size_t nCell = cellNum();
    initAddress.resize(nCell + 1);
    std::vector<int32_t> vertexCount(nCell + 1, 0);
    for (int32_t i = 0; i < nPoint; i++) {
      vertexCount[pointCellIndex[i]] += 1;
    }

    initAddress[0] = 0;
    int32_t sz_all = vertexCount[0];
    for (int32_t i = 0; i < nCell; i++) {
      initAddress[i + 1] = initAddress[i] + vertexCount[i];
      sz_all += vertexCount[i];
    }

    int32_t invalid_idx = -1;
    cellPointMap.resize(sz_all + 1, invalid_idx);
    std::vector<int32_t> node_pos(nCell + 1, 0);

    for (int32_t i = 0; i < nPoint; i++) {
      int32_t initAdd = initAddress[pointCellIndex[i]];
      int32_t pos = node_pos[pointCellIndex[i]] + initAdd;
      cellPointMap[pos] = i;
      node_pos[pointCellIndex[i]]++;
    }
  }

  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::putBox2Cells(size_t nBox, BoxFunction boxFun, std::vector<std::vector<size_t>> &pointCellIndex) const
  {
    pointCellIndex.resize(this->cellNum());
    int32_t numBox = (int32_t)nBox;
    std::vector<std::vector<Int2>> dataCollector;
#if defined(_USE_OMP)
    auto nThread = omp_get_max_threads();
    dataCollector.resize(nThread);
#else
    dataCollector.resize(1);
#endif

    if (N == 3) {
#if defined(_USE_OMP)
#pragma omp parallel for
#endif
      for (int32_t i0 = 0; i0 < numBox; i0++) {
        tBBox<Point, N> box;
        if (!boxFun(i0, box))
          continue;
        IntN iMin = getCellIndex(box.min().data());
        IntN iMax = getCellIndex(box.max().data());
        IntN idx;
#if defined(_USE_OMP)
        std::vector<Int2> &localData = dataCollector[omp_get_thread_num()];
#else
        std::vector<Int2> &localData = dataCollector[0];
#endif
        for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
          idx[0] = i;
          for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
            idx[1] = j;
            for (typename IntN::value_type k = iMin[2]; k <= iMax[2]; k++) {
              idx[2] = k;
              Int2::value_type pos = (Int2::value_type)index(idx);
              localData.push_back(Int2(pos, i0));
            }
          }
        }
      }
    }
    else if (N == 2) {
#if defined(_USE_OMP)
#pragma omp parallel for
#endif
      for (int32_t i0 = 0; i0 < numBox; i0++) {
        tBBox<Point, N> box;
        if (!boxFun(i0, box))
          continue;
        IntN iMin = getCellIndex(box.min().data());
        IntN iMax = getCellIndex(box.max().data());
#if defined(_USE_OMP)
        std::vector<Int2> &localData = dataCollector[omp_get_thread_num()];
#else
        std::vector<Int2> &localData = dataCollector[0];
#endif
        IntN idx;
        for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
          idx[0] = i;
          for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
            idx[1] = j;
            typename Int2::value_type pos = (Int2::value_type)index(idx);
            localData.push_back(Int2(pos, i0));
          }
        }
      }
    }

    for (auto &it : dataCollector) {
      for (auto &it1 : it) {
        pointCellIndex[it1[0]].push_back(it1[1]);
      }
    }
  }

  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::putBox2Cells(size_t nBox, BoxFunction boxFun, DInt2Array &pointCellIndex) const
  {
    size_t nCell = cellNum();
    pointCellIndex.clear();
    pointCellIndex.reserve(nCell);

    std::vector<std::vector<Int2>> dataCollector;
#if defined(_USE_OMP)
    auto nThread = omp_get_max_threads();
    dataCollector.resize(nThread);
#else
    dataCollector.resize(1);
#endif

    if (N == 3) {
#if defined(_USE_OMP)
#pragma omp parallel for
#endif
      for (int32_t i0 = 0; i0 < (int32_t)nBox; i0++) {
        tBBox<Point, N> box;
        if (!boxFun(i0, box)) {
          continue;
        }
        IntN iMin = getCellIndex(box.min().data());
        IntN iMax = getCellIndex(box.max().data());
        IntN idx;
#if defined(_USE_OMP)
        std::vector<Int2> &localData = dataCollector[omp_get_thread_num()];
#else
        std::vector<Int2> &localData = dataCollector[0];
#endif
        for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
          idx[0] = i;
          for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
            idx[1] = j;
            for (typename IntN::value_type k = iMin[2]; k <= iMax[2]; k++) {
              idx[2] = k;
              typename Int2::value_type pos = (Int2::value_type)(index(idx));
              localData.push_back(Int2(pos, (Int2::value_type)i0));
            }
          }
        }
      }
    }
    else if (N == 2) {
#if defined(_USE_OMP)
#pragma omp parallel for
#endif
      for (int32_t i0 = 0; i0 < (int32_t)nBox; i0++) {
        tBBox<Point, N> box;
        if (!boxFun(i0, box)) {
          continue;
        }
        IntN iMin = getCellIndex(box.min().data());
        IntN iMax = getCellIndex(box.max().data());
#if defined(_USE_OMP)
        std::vector<Int2> &localData = dataCollector[omp_get_thread_num()];
#else
        std::vector<Int2> &localData = dataCollector[0];
#endif
        IntN idx;
        for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
          idx[0] = i;
          for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
            idx[1] = j;

            Int2::value_type pos = (Int2::value_type)(index(idx));
            localData.push_back(Int2(pos, (Int2::value_type)i0));
          }
        }
      }
    }

    for (auto &it : dataCollector) {
      pointCellIndex.reserve(it.size() + pointCellIndex.size());
      for (auto &it1 : it) {
        pointCellIndex.push_back(it1);
      }
    }
    std::sort(pointCellIndex.begin(), pointCellIndex.end());
  }


  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCollisionCellPairs(const DInt2Array& range, const DInt2Array& pointCellIndex,
      BoxFunction boxFun, std::function<void(const Int2 &)> pairCB) const
  {
    int32_t num = (int32_t)range.size();

    std::vector<boost::unordered_set<Int2>> dataCollector;
#if defined(_USE_OMP)
    auto nThread = omp_get_max_threads();
    dataCollector.resize(nThread);
#pragma omp parallel for
#else
    dataCollector.resize(1);
#endif
    for (int32_t i = 0; i < num; i++) {
      if ((range[i][1] - range[i][0]) < 2)
        continue;
#if defined(_USE_OMP)
      auto &localData = dataCollector[omp_get_thread_num()];
#else
      auto &localData = dataCollector[0];
#endif
      for (int32_t j = range[i][0]; j < range[i][1]; j++) {
        for (size_t k = j + 1; k < range[i][1]; k++) {
          Int2 apair(pointCellIndex[j][1], pointCellIndex[k][1]);
          if (apair[1] < apair[0])
            std::swap(apair[0], apair[1]);
          if (localData.find(apair) == localData.end()) {
            tBBox<Point, N> box[2];
            boxFun(apair[0], box[0]);
            boxFun(apair[1], box[1]);
            if (box[0].overlaps(box[1])) {
              localData.insert(apair);
            }
          }
        }
      }
    }

    boost::unordered_set<Int2> dataSet;
    for (auto &it : dataCollector) {
      for (auto &it2 : it) {
        if (dataSet.insert(it2).second) {
          pairCB(it2);
        }
      }
    }
  }


  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCollisionCellPairs(size_t nBox, BoxFunction boxFun, std::function<void(const Int2 &)> pairCB) const
  {
    if (nBox < 1)
      return;

    auto withCellScaning = [&]() {
      size_t nCell = cellNum();
      std::vector<std::vector<size_t>> pointCellIndex(nCell);
      this->putBox2Cells(nBox, boxFun, pointCellIndex);
      std::vector<boost::unordered_set<Int2>> dataCollector;
#if defined(_USE_OMP)
      auto nThread = omp_get_max_threads();
      dataCollector.resize(nThread);
#pragma omp parallel for
#else
      dataCollector.resize(1);
#endif
      for (int32_t i = 0; i < (int32_t)nCell; i++) {
        if (pointCellIndex[i].size() < 2)
          continue;
#if defined(_USE_OMP)
        auto &localData = dataCollector[omp_get_thread_num()];
#else
        auto &localData = dataCollector[0];
#endif
        for (size_t j = 0; j < pointCellIndex[i].size(); j++) {
          for (size_t k = j + 1; k < pointCellIndex[i].size(); k++) {
            Int2 apair(pointCellIndex[i][j], pointCellIndex[i][k]);
            if (apair[1] < apair[0])
              std::swap(apair[0], apair[1]);
            if (localData.find(apair) == localData.end()) {
              tBBox<Point, N> box[2];
              boxFun(apair[0], box[0]);
              boxFun(apair[1], box[1]);
              if (box[0].overlaps(box[1])) {
                localData.insert(apair);
              }
            }
          }
        }
      }

      boost::unordered_set<Int2> dataSet;
      for (auto &it : dataCollector) {
        for (auto &it2 : it) {
          if (dataSet.insert(it2).second) {
            pairCB(it2);
          }
        }
      }
    };

    auto withCellScaningSmart = [&]() {
      DInt2Array range, pointCellIndex;

      putBox2Cells(nBox, boxFun, pointCellIndex);
      size_t nCell = cellNum();
      range.reserve(nCell);
      Int2 intv;
      intv[0] = 0;
      for (size_t i = 1; i < pointCellIndex.size(); i++) {
        if (pointCellIndex[i][0] != pointCellIndex[intv[0]][0]) {
          intv[1] = (int32_t)i;
          range.push_back(intv);
          intv[0] = intv[1];
        }
      }
      this->getCollisionCellPairs(range, pointCellIndex, boxFun, pairCB);
    };

    withCellScaning();
    //withCellScaningSmart();
    // withCellScaning2();
  }


  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCollisionCellPairs(size_t nBox, BoxFunction boxFun, std::function<bool(const Int2 &)> pairCB, DInt2Array& intersected) const
  {
    if (nBox < 1)
      return;

    auto withCellScaning = [&]() {
      size_t nCell = cellNum();
      std::vector<std::vector<size_t>> pointCellIndex(nCell);
      this->putBox2Cells(nBox, boxFun, pointCellIndex);
      std::vector<boost::unordered_set<Int2>> dataCollector;
#if 1 //defined(_USE_OMP)
      auto nThread = omp_get_max_threads();
      dataCollector.resize(nThread);
#pragma omp parallel for
#else
      dataCollector.resize(1);
#endif
      for (int32_t i = 0; i < (int32_t)nCell; i++) {
        if (pointCellIndex[i].size() < 2)
          continue;
#if 1//defined(_USE_OMP)
        auto &localData = dataCollector[omp_get_thread_num()];
#else
        auto &localData = dataCollector[0];
#endif
        for (size_t j = 0; j < pointCellIndex[i].size(); j++) {
          for (size_t k = j + 1; k < pointCellIndex[i].size(); k++) {
            Int2 apair(pointCellIndex[i][j], pointCellIndex[i][k]);
            if (apair[1] < apair[0])
              std::swap(apair[0], apair[1]);
            if (localData.find(apair) == localData.end()) {
              tBBox<Point, N> box[2];
              boxFun(apair[0], box[0]);
              boxFun(apair[1], box[1]);
              if (box[0].overlaps(box[1])) {
                localData.insert(apair);
              }
            }
          }
        }
      }

      DInt2Array pairs;
      for (auto &it : dataCollector) {
        for (auto &it2 : it) {
          pairs.push_back(it2);
        }
      }
      _Private::Equal fe;
      _Private::Less fl;
      Unique<DInt2Array, _Private::Equal, _Private::Less>(pairs, fe, fl);

      int npairs = (int)(pairs.size());
      DShortArray interFlag(npairs, 0);

      if (1) {
#pragma omp parallel for
        for (int i = 0; i < npairs; i++) {
          if (pairCB(pairs[i])) {
            interFlag[i] = 1;
          }
        }
      }
      else
      {
        for (int i = 0; i < npairs; i++) {
          if (pairCB(pairs[i])) {
            interFlag[i] = 1;
          }
        }
      }

      for (int i = 0; i < npairs; i++)
      {
        if (interFlag[i] > 0)
        {
          intersected.push_back(pairs[i]);
        }
      }
    };

    withCellScaning();
    //withCellScaningSmart();
    // withCellScaning2();
  }


  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCollisionCellPairs(const DInt2Array& range, const DInt2Array& pointCellIndex,
      BoxFunction boxFun, std::function<bool(const Int2 &)> pairCB, DInt2Array& intersected) const
  {
    int32_t num = (int32_t)range.size();

    std::vector<boost::unordered_set<Int2>> dataCollector;
#if 1//defined(_USE_OMP)
    auto nThread = omp_get_max_threads();
    dataCollector.resize(nThread);
#pragma omp parallel for
#else
    dataCollector.resize(1);
#endif
    for (int32_t i = 0; i < num; i++) {
      if ((range[i][1] - range[i][0]) < 2)
        continue;
#if 1 //defined(_USE_OMP)
      auto &localData = dataCollector[omp_get_thread_num()];
#else
      auto &localData = dataCollector[0];
#endif
      for (int32_t j = range[i][0]; j < range[i][1]; j++) {
        for (size_t k = j + 1; k < range[i][1]; k++) {
          Int2 apair(pointCellIndex[j][1], pointCellIndex[k][1]);
          if (apair[1] < apair[0])
            std::swap(apair[0], apair[1]);
          if (localData.find(apair) == localData.end()) {
            tBBox<Point, N> box[2];
            boxFun(apair[0], box[0]);
            boxFun(apair[1], box[1]);
            if (box[0].overlaps(box[1])) {
              localData.insert(apair);
            }
          }
        }
      }
    }

    DInt2Array pairs;
    for (auto &it : dataCollector) {
      for (auto &it2 : it) {
        pairs.push_back(it2);
      }
    }
    _Private::Equal fe;
    _Private::Less fl;
    Unique<DInt2Array, _Private::Equal, _Private::Less>(pairs, fe, fl);

    int npairs = (int)(pairs.size());
    DShortArray interFlag(npairs, 0);

    if (1) {
#pragma omp parallel for
      for (int i = 0; i < npairs; i++) {
        if (pairCB(pairs[i])) {
          interFlag[i] = 1;
        }
      }
    }
    else
    {
      for (int i = 0; i < npairs; i++) {
        if (pairCB(pairs[i])) {
          interFlag[i] = 1;
        }
      }
    }

    for (int i = 0; i < npairs; i++)
    {
      if (interFlag[i] > 0)
      {
        intersected.push_back(pairs[i]);
      }
    }
  }

  template <typename T, uint32_t N>
  void
    tSpaceHash<T, N>::getCollisionCellPairs(size_t nBox0,
      BoxFunction boxFun0,
      size_t nBox1,
      BoxFunction boxFun1,
      std::function<void(const Int2 &)> pairCB) const
  {
    std::vector<std::vector<size_t>> pointCellIndex;
    this->putBox2Cells(nBox0, boxFun0, pointCellIndex);

    std::vector<boost::unordered_set<Int2>> dataCollector;
#if defined(_USE_OMP)
    auto nThread = omp_get_max_threads();
    dataCollector.resize(nThread);
#pragma omp parallel for
#else
    dataCollector.resize(1);
#endif
    for (int32_t ib = 0; ib < (int32_t)nBox1; ib++) {
      tBBox<Point, N> box;
      boxFun1(ib, box);
      IntN iMin = getCellIndex(box.min().data());
      IntN iMax = getCellIndex(box.max().data());
      IntN idx;
#if defined(_USE_OMP)
      auto &localData = dataCollector[omp_get_thread_num()];
#else
      auto &localData = dataCollector[0];
#endif
      for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
        idx[0] = i;
        for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
          idx[1] = j;
          for (typename IntN::value_type k = iMin[2]; k <= iMax[2]; k++) {
            idx[2] = k;
            typename Int2::value_type pos = (Int2::value_type)(index(idx));

            for (auto &it : pointCellIndex[pos]) {
              localData.insert(Int2((Int2::value_type)it, ib));
            }
          }
        }
      }
    }
    std::vector<std::vector<bool>> checker(nBox0);
    for (size_t i = 0; i < nBox0; i++) {
      checker[i].resize(nBox1, false);
    }
    for (auto &it : dataCollector) {
      for (auto &it2 : it) {
        if (checker[it2[0]][it2[1]]) {
          continue;
        }
        checker[it2[0]][it2[1]] = true;
        pairCB(it2);
      }
    }
  }

  template <typename T, uint32_t N>
  size_t
    tSpaceHash<T, N>::cellNum() const
  {
    int32_t num = mSize[0];
    for (uint32_t i = 1; i < N; i++) {
      num *= mSize[i];
    }
    return num;
  }

  template <typename T, uint32_t N>
  size_t
    tSpaceHash<T, N>::index(const IntN &idx) const
  {
    if (N == 3) {
      return idx[0] + idx[1] * mSize[0] + idx[2] * mSize[0] * mSize[1];
    }
    else if (N == 2) {
      return idx[0] + idx[1] * mSize[0];
    }
    else {
      assert(0);
    }
    return -1;
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  template <typename T, uint32_t N>
  tPointHash<T, N>::tPointHash(size_t nPoint,
    typename tSpaceHash<T, N>::PointFunction pointFun,
    uint32_t resolution,
    double enlageFactor)
    : mNumPoint(nPoint)
    , mPointFun(pointFun)
  {
    this->create(nPoint, pointFun, resolution, enlageFactor);
    this->getCellPointMap(nPoint, pointFun, mInitAddress, mCellPointMap);
  }

  template <typename T, uint32_t N>
  void
    tPointHash<T, N>::getIncludes(const typename tSpaceHash<T, N>::ValueType *pmin0,
      const typename tSpaceHash<T, N>::ValueType *pmax0,
      std::function<void(int32_t)> pointIndexFun) const
  {
    typename tSpaceHash<T, N>::IntN iMin = this->getCellIndex(pmin0);
    typename tSpaceHash<T, N>::IntN iMax = this->getCellIndex(pmax0);
    typename tSpaceHash<T, N>::ValueType ipoint[N];
    typename tSpaceHash<T, N>::IntN idx;
    if (N == 3) {
      for (auto i = iMin[0]; i <= iMax[0]; i++) {
        idx[0] = i;
        for (auto j = iMin[1]; j <= iMax[1]; j++) {
          idx[1] = j;
          for (auto k = iMin[2]; k <= iMax[2]; k++) {
            idx[2] = k;
            auto pos = tSpaceHash<T, N>::index(idx);
            for (auto i0 = mInitAddress[pos]; i0 < mInitAddress[pos + 1]; i0++) {
              mPointFun(mCellPointMap[i0], ipoint);
              bool isIncluded = true;
              for (auto i1 = 0; i1 < N; i1++) {
                if (ipoint[i1] < pmin0[i1] || ipoint[i1] > pmax0[i1]) {
                  isIncluded = false;
                  break;
                }
              }
              if (isIncluded) {
                pointIndexFun(mCellPointMap[i0]);
              }
            }
          }
        }
      }
    }
    else {
      for (auto i = iMin[0]; i <= iMax[0]; i++) {
        idx[0] = i;
        for (auto j = iMin[1]; j <= iMax[1]; j++) {
          idx[1] = j;
          {
            auto pos = tSpaceHash<T, N>::index(idx);
            for (auto i0 = mInitAddress[pos]; i0 < mInitAddress[pos + 1]; i0++) {
              mPointFun(mCellPointMap[i0], ipoint);
              bool isIncluded = true;
              for (auto i1 = 0; i1 < N; i1++) {
                if (ipoint[i1] < pmin0[i1] || ipoint[i1] > pmax0[i1]) {
                  isIncluded = false;
                  break;
                }
              }
              if (isIncluded) {
                pointIndexFun(mCellPointMap[i0]);
              }
            }
          }
        }
      }
    }
  }

  template <typename T, uint32_t N>
  int32_t
    tPointHash<T, N>::getNearest(const ValueType *p, std::function<bool(size_t)> filter) const
  {
    if (!p)
      return -1;
    std::vector<bool> checked(this->cellNum());
    auto idx = this->getCellIndex(p);
    typename tSpaceHash<T, N>::IntN idxi;
    auto up = idx, low = idx;
    T disMin = std::numeric_limits<T>::max();
    int32_t pointIdx = -1;
    auto Check = [&](typename tSpaceHash<T, N>::IntN &idx) {
      auto bi = this->index(idx);
      if (checked[bi])
        return;
      checked[bi] = true;
      Tuple0<T, N> pi;
      for (auto i = mInitAddress[bi]; i < mInitAddress[bi + 1]; i++) {
        if (filter(mCellPointMap[i]))
          continue;
        mPointFun(mCellPointMap[i], pi.data());
        T dis2 = pi.distanceSqr(p);
        if (dis2 < disMin) {
          disMin = dis2;
          pointIdx = (int32_t)mCellPointMap[i];
        }
      }
    };

    for (;;) {
      up += 1;
      low -= 1;
      bool allChecked = true;
      for (uint32_t i = 0; i < N; i++) {
        if (up[i] < (this->mSize[i] - 1) || low[i] >= 0)
        {
          allChecked = false;
        }
        up[i] = std::min(up[i], this->mSize[i] - 1);
        low[i] = std::max(0, low[i]);
      }

      if (N == 2) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            Check(idxi);
          }
        }
      }
      else if (N == 3) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            for (auto iz = low[2]; iz <= up[2]; iz++) {
              idxi[2] = iz;
              Check(idxi);
            }
          }
        }
      }
      else {
        assert(0);
        return -1;
      }
      if (pointIdx >= 0 || allChecked)
        return pointIdx;
    }
    return pointIdx;
  }

  template <typename T, uint32_t N>
  void  tPointHashUniqueInserter<T, N>::setPointLocalIndex(size_t i, int idx)
  {
    if (i >= mPoints.size())
    {
      assert(0);
      return;
    }
    if (mLocalIndex.size() < mPoints.size())
    {
      mLocalIndex.resize(mPoints.size(), -1);
    }
    mLocalIndex[i] = idx;
  }

  template <typename T, uint32_t N>
  int32_t
    tPointHashUniqueInserter<T, N>::hasPointIn(const Point &points, ValueType eps, std::function<bool(int32_t)> filter)const
  {
    assert(eps > Gloabals::ZERO);
    if (mPoints.empty()) {
      return -1;
    }

    T eps2 = eps * eps;
    Tuple0<T, N> ptMin(points), ptMax(points);
    ptMin -= eps;
    ptMax += eps;
    IntN low = this->getCellIndex(ptMin);
    IntN up = this->getCellIndex(ptMax);
    IntN idxi;
    if (N == 2) {
      for (auto ix = low[0]; ix <= up[0]; ix++) {
        idxi[0] = ix;
        for (auto iy = low[1]; iy <= up[1]; iy++) {
          idxi[1] = iy;
          auto bi = this->index(idxi);
          for (auto it : mBructPointList[bi]) {
            if (!filter(it) && mPoints[it].distanceSqr(points) < eps2) {
              return (int32_t)it;
            }
          }
        }
      }
    }
    else if (N == 3) {
      for (auto ix = low[0]; ix <= up[0]; ix++) {
        idxi[0] = ix;
        for (auto iy = low[1]; iy <= up[1]; iy++) {
          idxi[1] = iy;
          for (auto iz = low[2]; iz <= up[2]; iz++) {
            idxi[2] = iz;
            auto bi = this->index(idxi);
            for (auto it : mBructPointList[bi]) {
              if (!filter(it) && mPoints[it].distanceSqr(points) < eps2) {
                return (int32_t)it;
              }
            }
          }
        }
      }
    }
    else {
      assert(0);
      return -1;
    }

    return -1;
  }

  template <typename T, uint32_t N>
  int32_t
    tPointHashUniqueInserter<T, N>::nearest(const Point &point, std::function<bool(int32_t)> filter)const
  {
    if (mPoints.empty()) {
      return -1;
    }

    IntN cur = this->getCellIndex(point);
    int32_t offset = 0;
    while (1)
    {
      offset += 1;
      IntN low = cur;
      IntN up = cur;
      for (uint32_t i = 0; i < N; i++)
      {
        low[i] = std::max(0, cur[i] - offset);
        up[i] = std::min(mSize[i] - 1, up[i] + offset);
      }

      double disMin = DBL_MAX;
      int32_t whi = -1;

      IntN idxi;
      if (N == 2) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            auto bi = this->index(idxi);
            for (auto it : mBructPointList[bi]) {
              if (filter(it)) {
                continue;
              }

              double dis = mPoints[it].distanceSqr(points);
              if (dis < disMin)
              {
                disMin = dis;
                whi = it;
                if (disMin < Globals::ZERO) {
                  return whi;
                }
              }
            }
          }
        }
      }
      else if (N == 3) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            for (auto iz = low[2]; iz <= up[2]; iz++) {
              idxi[2] = iz;
              auto bi = this->index(idxi);
              for (auto it : mBructPointList[bi]) {
                if (filter(it)) {
                  continue;
                }

                double dis = mPoints[it].distanceSqr(points);
                if (dis < disMin)
                {
                  disMin = dis;
                  whi = it;
                  if (disMin < Globals::ZERO) {
                    return whi;
                  }
                }
              }
            }
          }
        }
      }
      else {
        assert(0);
        return -1;
      }

      if (whi >= 0) {
        return whi;
      }
    }
    return -1;
  }

  template <typename T, uint32_t N>
  int32_t
    tPointHashUniqueInserter<T, N>::insertPoint(const Point &points, ValueType eps)
  {
    if (eps <= 0) {
      eps = this->mStep / 100;
    }

    auto pos = this->hasPointIn(points, eps, [](int32_t) {return false; });
    if (pos >= 0)
      return pos;

    return this->appendPoint(points);
  }

  template <typename T, uint32_t N>
  int32_t
    tPointHashUniqueInserter<T, N>::appendPoint(const Point &point)
  {
    if (mBructPointList.empty())
    {
      mBructPointList.resize(this->cellNum());
    }
    auto idx = this->getCellIndex(point.data());
    auto bi = this->index(idx);
    mBructPointList[bi].push_back(mPoints.size());
    mPoints.push_back(Tuple0<T, N>(point));
    return (int32_t)mPoints.size() - 1;
  }

  template <typename T, uint32_t N>
  tBoxHash<T, N>::tBoxHash(size_t nBox, BoxFunction boxFun,
    uint32_t resolution,
    double enlageFactor)
    :mBoxFun(boxFun)
  {
    if (nBox < 2 || enlageFactor <= 1.0)
      return;

    uint32_t num = (uint32_t) std::pow(nBox, 1.0 / 3.0);
    resolution = std::min(resolution, std::max(num, uint32_t(3)));
    Tuple0<ValueType, N> bmin(std::numeric_limits<ValueType>::max());
    Tuple0<ValueType, N> bmax(-std::numeric_limits<ValueType>::max());
    for (size_t i = 0; i < nBox; i++) {
      Box pos;
      if (!boxFun(i, pos))
        continue;
      for (uint32_t k = 0; k < N; k++) {
        bmin[k] = std::min(bmin[k], pos.min()[k]);
        bmax[k] = std::max(bmax[k], pos.min()[k]);
        bmin[k] = std::min(bmin[k], pos.max()[k]);
        bmax[k] = std::max(bmax[k], pos.max()[k]);
      }
    }

    this->create(bmin, bmax, resolution, enlageFactor, 0.0);
    this->putBox2Cells(nBox, boxFun, mPointCellIndex);
    size_t nCell = this->cellNum();
    mRange.reserve(nCell);
    Int2 intv;
    intv[0] = 0;
    for (size_t i = 1; i < mPointCellIndex.size(); i++) {
      if (mPointCellIndex[i][0] != mPointCellIndex[intv[0]][0]) {
        intv[1] = (int32_t)i;
        mRange.push_back(intv);
        intv[0] = intv[1];
      }
    }
  }




  template <typename T, uint32_t N>
  void tBoxHash<T, N>::getIncludes(const typename tSpaceHash<T, N>::ValueType *pt,
    std::function<bool(int32_t)> boxFun) const
  {
    IntSet boxSet;
    this->getIncludes(this->getCellIndex(pt), boxFun, boxSet);
  }

  template <typename T, uint32_t N>
  void tBoxHash<T, N>::getBoxIncludes(const typename tSpaceHash<T, N>::ValueType *pmin,
    const typename tSpaceHash<T, N>::ValueType *pmax,
    std::function<bool(int32_t)> boxFun) const
  {
    IntN iMin = this->getCellIndex(pmin);
    IntN iMax = this->getCellIndex(pmax);
    IntN idx;

    IntSet boxSet;
    if (N == 3) {
      for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
        idx[0] = i;
        for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
          idx[1] = j;
          for (typename IntN::value_type k = iMin[2]; k <= iMax[2]; k++)
          {
            idx[2] = k;
            this->getIncludes(idx, boxFun, boxSet);
          }
        }
      }
    }
    else {
      assert(N == 2);
      for (typename IntN::value_type i = iMin[0]; i <= iMax[0]; i++) {
        idx[0] = i;
        for (typename IntN::value_type j = iMin[1]; j <= iMax[1]; j++) {
          idx[1] = j;
          this->getIncludes(idx, boxFun, boxSet);
        }
      }
    }
  }


  template <typename T, uint32_t N>
  void tBoxHash<T, N>::getSelfOverlap(std::function<void(const Int2 &)> pairCB) const
  {
    this->getCollisionCellPairs(mRange, mPointCellIndex, mBoxFun, pairCB);
  }

  template <typename T, uint32_t N>
  void tBoxHash<T, N>::getSelfOverlap(std::function<bool(const Int2 &)> pairCB, DInt2Array& intersected) const
  {
    this->getCollisionCellPairs(mRange, mPointCellIndex, mBoxFun, pairCB, intersected);
  }

  template <typename T, uint32_t N>
  void tBoxHash<T, N>::getIncludes(const IntN& idx, std::function<bool(int32_t)> boxFun, IntSet& boxSet) const
  {
    Int2::value_type pos = (Int2::value_type)this->index(idx);
    auto it = std::lower_bound(mPointCellIndex.begin(), mPointCellIndex.end(), Int2(pos, -1),
      [](const Int2& a, const Int2& b) {return a[0] < b[0]; });

    if (it != mPointCellIndex.end() && it->at(0) == pos)
    {
      for (; it != mPointCellIndex.end(); it++)
      {
        if (it->at(0) != pos)
          break;
        if (boxSet.insert(it->at(1)).second)
        {
          if (!boxFun(it->at(1))) {
            return;
          }
        }
      }
    }
  }

  template <typename T, uint32_t N>
  void  tBoxHash<T, N>::getIntersected(const Point& from, const Point& to, std::function<bool(int32_t)> functor)const
  {
    int num = std::max(1.0, sqrt(from.distanceSqr(to)) / this->mStep);
    double step = 1.0 / num;
    double t = 0.0;
    auto lst = from;
    std::vector<int32_t> dataArray;
    for (int i = 0; i < num; i++)
    {
      t += step;
      auto mid = from * (1.0 - t) + to * t;
      Box abox;
      abox += lst;
      abox += mid;

      this->getBoxIncludes(abox.pmin_max[0].data(), abox.pmin_max[1].data(), [&dataArray](int32_t whi) {
        dataArray.push_back(whi);
        return true;
      });
      lst = mid;
    }

    Unique(dataArray);
    for (auto it : dataArray)
    {
      if (!functor(it)) {
        break;
      }
    }
  }

  template <typename T, uint32_t N>
  int32_t tBoxHash<T, N>::getNearest(const ValueType *p, std::function<bool(size_t)> boxFun) const
  {
    if (!p)
      return -1;
    std::vector<bool> checked(this->cellNum());
    auto idx = this->getCellIndex(p);
    typename tSpaceHash<T, N>::IntN idxi;
    auto up = idx, low = idx;
    T disMin = std::numeric_limits<T>::max();
    int32_t pointIdx = 0, iStop = 0;
    auto Check = [&](typename tSpaceHash<T, N>::IntN &idx) {
      auto pos = this->index(idx);
      if (checked[pos])
        return;
      checked[pos] = true;
      auto it = std::lower_bound(mPointCellIndex.begin(), mPointCellIndex.end(), Int2(pos, -1),
        [](const Int2& a, const Int2& b) {return a[0] < b[0]; });

      if (it != mPointCellIndex.end() && it->at(0) == pos)
      {
        for (; it != mPointCellIndex.end(); it++)
        {
          if (it->at(0) != pos)
            break;
          ++pointIdx;
          {
            if (!boxFun(it->at(1)))
            {
              iStop = true;
              return;
            }
          }
        }
      }
    };

    for (;;) {
      up += 1;
      low -= 1;
      bool allChecked = true;
      for (uint32_t i = 0; i < N; i++) {
        if (up[i] < (this->mSize[i] - 1) || low[i] >= 0)
        {
          allChecked = false;
        }
        up[i] = std::min(up[i], this->mSize[i] - 1);
        low[i] = std::max(0, low[i]);
      }

      pointIdx = 0;

      if (N == 2) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            Check(idxi);
            if (iStop > 0)
            {
              return pointIdx;
            }
          }
        }
      }
      else if (N == 3) {
        for (auto ix = low[0]; ix <= up[0]; ix++) {
          idxi[0] = ix;
          for (auto iy = low[1]; iy <= up[1]; iy++) {
            idxi[1] = iy;
            for (auto iz = low[2]; iz <= up[2]; iz++) {
              idxi[2] = iz;
              Check(idxi);
              if (iStop > 0)
              {
                return pointIdx;
              }
            }
          }
        }
      }
      else {
        assert(0);
        return -1;
      }
      if (pointIdx > 0 || allChecked)
        return pointIdx;
    }
    return pointIdx;
  }

};


