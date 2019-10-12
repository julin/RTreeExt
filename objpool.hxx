
#pragma once

template <class T,size_t __PAGE_SIZE_=1024>
class  ObjectPool
{
  struct Page
  {
    enum
    {
      SZ=__PAGE_SIZE_/CHAR_BIT+1,
    };
    T             data[__PAGE_SIZE_];
    unsigned char flag[SZ];

    Page()
    {
      memset(flag,0,SZ*sizeof(unsigned char));
    }

    void where_is(int i,int& whi,int& off)const {whi=i/CHAR_BIT;off=i%CHAR_BIT;}
    void set(int i,bool istrue=true) 
    { 
      int whi,off;
      this->where_is(i,whi,off);
      unsigned  bi=1<<off;
      if(istrue)
        flag[whi]|=bi;
      else
        flag[whi]&=~bi;
    }

    bool is(unsigned i) const 
    { 
      int whi,off;
      this->where_is(i,whi,off);
      return (flag[whi]&(1<<off))>0; 
    }

    size_t erase(size_t pos)
    {
      if (pos>=__PAGE_SIZE_) return  size_t(-1);
      if(!is(pos))
      {
        set(pos,true);
        return pos;
      }
      return  size_t(-1);
    }

    bool  valid(size_t i)const {return !is(i);}
  };


  struct Sorter0
  {
    bool operator()(const T* dat0,const T* dat1)const
    {
      return dat0<dat1;
    }
    bool operator()(const T* dat0,const Page* dat1)const
    {
      return dat0<dat1->data;
    }

    bool operator()(const Page* dat1,const T* dat0)const
    {
      return ((dat0-dat1->data)>=__PAGE_SIZE_);
    }
  };

  struct Sorter
  {
    bool operator()(size_t i,size_t j)const
    {
      return _pages->at(i)<_pages->at(j);
    }

    Sorter(const std::vector<Page*>* pages)
      :_pages(pages)
    {
    }
    const std::vector<Page*>*        _pages;
  };

  enum 
  {
    INIT_PAGES_SORTED   =1
    ,SORTEDPAGEOK       =1<<1
    ,UNDEFINED_STATE    =1<<2
  };

 
  std::vector<Page*>        _pages;
  std::vector<size_t>       _erased;
  size_t                   _size;
  
  mutable size_t*          _sortedpagesPtr;
  mutable unsigned         _flag;

private:
  size_t  getPage(size_t pos,size_t& pidx)const
  {
    size_t pn=pos/__PAGE_SIZE_;
    if(pn>=_pages.size()) return size_t(-1);
    pidx=pos%__PAGE_SIZE_;
    return pn;
  }

  size_t  getPage(const T* data,size_t& pidx)const
  {
    size_t ps=_pages.size();
    if(ps==0) return size_t(-1);
    Sorter0 compare;
    if(_flag& INIT_PAGES_SORTED)
    {
      int begin=0,end=_pages.size();
      while(begin<=end)
      {
        int mid=(end+begin)>>1;
        if(ps<=mid) return size_t(-1);
        if(compare(data,_pages[mid]))end=mid-1;
        else if(compare(_pages[mid],data)) begin=mid+1;
        else 
        {
          assert(data>=_pages[mid]->data);
          pidx=data-_pages[mid]->data;
          assert(pidx<__PAGE_SIZE_);
          return mid;
        }
      }
    }else
    {
      if(0==(_flag&SORTEDPAGEOK))
      {
        delete[] _sortedpagesPtr;
        size_t sz=_pages.size();
        _sortedpagesPtr=new size_t[sz];

        for(size_t i=0;i<sz;i++)
          _sortedpagesPtr[i]=i;

        Sorter sorter(&_pages);
        std::sort(_sortedpagesPtr,_sortedpagesPtr+sz,sorter);
        _flag|=SORTEDPAGEOK;
      }

      int begin=0,mid=0,end=_pages.size();
      while(begin<=end)
      {
        mid=(end+begin)>>1;
        int pos=_sortedpagesPtr[mid];
        if(ps<=pos) return size_t(-1);
        if(compare(data,_pages[pos]))end=mid-1;
        else if(compare(_pages[pos],data)) begin=mid+1;
        else 
        {
          assert(data>=_pages[pos]->data);
          pidx=data-_pages[pos]->data;
          assert(pidx<__PAGE_SIZE_);
          return pos;
        }
      }
    }
    return size_t(-1);
  }
public:
  ObjectPool()
    :_size(0),_sortedpagesPtr(NULL)
  {
    _flag =INIT_PAGES_SORTED;
  }

  ~ObjectPool()
  {
    clear();
  }

  void clear()
  {
    size_t i,sz=_pages.size();
    if(sz==0) return;
    for(i=0;i<sz;i++)
    {
      delete _pages[i];
    }
    _pages.clear();
    _erased.clear();
    _size=0;
    _flag |=INIT_PAGES_SORTED;
    _flag &=SORTEDPAGEOK;
    delete[] _sortedpagesPtr;
    _sortedpagesPtr=0;
  }

  T*  mAlloc(const T& obj=T())
  {
    size_t pi,off;

    if(!_erased.empty())
    {
      size_t pos=_erased[_erased.size()-1];
      _erased.pop_back();
      pi= getPage(pos,off);
    }else 
    {
      pi= getPage(_size,off);
      Sorter0 compare;
      if(pi==size_t(-1))
      {
        Page* p=new Page();
        assert(p!=NULL);
        off=0;
        pi=_pages.size();
        _pages.push_back(p);
        if((_flag&INIT_PAGES_SORTED) && pi>=1 && 
          compare(_pages[pi]->data,_pages[pi-1]->data))
          _flag &=~INIT_PAGES_SORTED;
        _flag &=~SORTEDPAGEOK;
      }
    }
    _pages[pi]->data[off]=obj;
    _pages[pi]->set(off,false);
    _size++;
    return &(_pages[pi]->data[off]);
  }

  bool  mFree(const T* data)
  {
    if(data==0) 
      return false;
    size_t pidx,pn,pos;
    pn=getPage(data,pidx);
    if(pn==size_t(-1)) 
      return false; 
    pos=pn*__PAGE_SIZE_+pidx;
    size_t rt=_pages[pn]->erase(pidx);
    if(rt==size_t(-1))
      return false;
    _size-=1;
    _erased.push_back(pos);
    return true;
  }

  bool  has(const void* data)const
  {
    if(data==0) 
      return false;
    size_t pidx,pn;
    pn=getPage((const T*)data,pidx);
    if(pn==size_t(-1))
      return false; 
    return _pages[pn]->valid(pidx);
  }

  std::vector<const T*> getValidData()const
  {
    std::vector<const T*> dataArray;
    size_t pi, off;
    for (size_t i = 0; i < _size; i++)
    {
      pi = getPage(_size, off);
      if (pi != size_t(-1))
      {
        if (_pages[pi]->valid(off))
        {
          dataArray.push_back(_pages[pi]->data + off);
        }
      }
    }
    return dataArray;
  }
};




