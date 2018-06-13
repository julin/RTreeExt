//mext_rtree.hxx
#pragma once

#include "rtree.h"

template<class DATATYPE, class ELEMTYPE, int NUMDIMS, 
class ELEMTYPEREAL = ELEMTYPE, int TMAXNODES = 6, int TMINNODES = 2>
class RTreeExt : public RTree<DATATYPE,ELEMTYPE,NUMDIMS,ELEMTYPEREAL,TMAXNODES,TMINNODES>
{
protected:
  typedef std::vector<DATATYPE>   DTArray;
  typedef RTree<DATATYPE,ELEMTYPE,NUMDIMS,ELEMTYPEREAL,TMAXNODES,TMINNODES> Base;
  typedef typename Base::Node  Node;
  typedef typename Base::Branch  Branch;
  typedef typename Base::Rect  Rect;
  struct NodePair
  {
    NodePair( Node *node0, Node *node1 ) : node_a( node0 ), node_b( node1 ){}   

    Node *node_a;
    Node *node_b;
  };

  /// Minimal bounding rectangle (n-dimensional)

  ELEMTYPE SquaredDistSup(const Rect& rect, const ELEMTYPE *point )const
  {
    int k = 0;
    ELEMTYPE mean[NUMDIMS];

    for ( k = 0; k < NUMDIMS; ++k )
      mean[k] = ( rect.m_min[k] + rect.m_max[k] ) * 0.5;

    ELEMTYPE minDist = 0.0;

    for ( k = 0; k < NUMDIMS; ++k )
    {
      ELEMTYPE rmkDiff = ( point[k] <= mean[k] ) ? ( point[k] - rect.m_min[k] ) : ( point[k] - rect.m_max[k] ); 

      ELEMTYPE value1st = rmkDiff * rmkDiff;

      ELEMTYPE value2nd = 0.0;

      for ( int j = 0; j < NUMDIMS; ++j )
      {
        if ( j == k )
          continue;

        ELEMTYPE rMjDiff = ( point[j] >= mean[j] ) ? ( point[j] - rect.m_min[j] ) : ( point[j] - rect.m_max[j] ); 
        value2nd += rMjDiff * rMjDiff;
      }

      ELEMTYPE currValue = value1st + value2nd;  

      if ( k==0 || currValue < minDist)
        minDist = currValue;
    }

    return minDist;
  }


  ELEMTYPE SquaredDistSup(const Rect& rect, const Rect& rect1 )const
  {
    int k = 0;
    ELEMTYPE pmin[NUMDIMS],pmax[NUMDIMS];
    memcpy(pmin,rect.m_min,NUMDIMS*sizeof(ELEMTYPE));
    memcpy(pmax,rect.m_max,NUMDIMS*sizeof(ELEMTYPE));

    for ( k = 0; k < NUMDIMS; ++k )
    {
      if(rect1.m_min[k]<pmin[k])
        pmin[k]=rect1.m_min[k];
      if(rect1.m_min[k]>pmax[k])
        pmax[k]=rect1.m_min[k];

      if(rect1.m_max[k]<pmin[k])
        pmin[k]=rect1.m_max[k];
      if(rect1.m_max[k]>pmax[k])
        pmax[k]=rect1.m_max[k];

    }

    ELEMTYPE dist = 0.0;

    for ( k = 0; k < NUMDIMS; ++k )
    {
      ELEMTYPE rmkDiff = pmax[k]- pmin[k];

      dist += (rmkDiff * rmkDiff);
    }

    return dist;
  }

  ELEMTYPE SquaredDistInf(const Rect& rect, const ELEMTYPE *point )const
  {
    ELEMTYPE sumDist = 0.0;

    for ( int index = 0; index < NUMDIMS; ++index )
    {
      ELEMTYPE currDiff = 0.0;

      if ( point[index] < rect.m_min[index] )
      {
        currDiff = point[index] - rect.m_min[index];

        sumDist += currDiff * currDiff;
      }
      else if ( point[index] > rect.m_max[index] )
      {
        currDiff = point[index] - rect.m_max[index];
        sumDist += currDiff * currDiff;
      }
    }

    return sumDist;
  }


  ELEMTYPE SquaredDistInf(const Rect& rect, const Rect &rect1 )const
  {
    ELEMTYPE sumDist = 0.0;

    for ( int i = 0; i < NUMDIMS; ++i )
    {
      ELEMTYPE currDiff = 0.0;

      if ( rect1.m_max[i] < rect.m_min[i] )
      {
        currDiff = rect1.m_max[i] - rect.m_min[i];

        sumDist += currDiff * currDiff;
      }
      else if ( rect1.m_min[i] > rect.m_max[i] )
      {
        currDiff = rect1.m_min[i] - rect.m_max[i];
        sumDist += currDiff * currDiff;
      }
    }

    return sumDist;
  }


  struct BranchInfo
  {
    BranchInfo() : index( -1 ), sqDistInf(0){}
    bool operator<( const BranchInfo &other ){ return ( this->sqDistInf < other.sqDistInf ); }
    bool operator<( const BranchInfo &other ) const { return ( this->sqDistInf < other.sqDistInf ); }

    int index;
    ELEMTYPE sqDistInf;
  };

  void GetRect(Node& node, Rect &rect )
  {
    int idx = 0;

    if ( !node.m_count )
    {
      for ( idx = 0; idx < NUMDIMS; idx++ )
        rect.m_min[idx] = rect.m_max[idx] = 0.0;

      return;
    }

    for ( int index = 0; index < node.m_count; index++ )
    {
      Rect &currRect = node.m_branch[index].m_rect;

      if ( !index )
      {
        rect = currRect;
      }
      else
      {
        for ( idx = 0; idx < NUMDIMS; idx++ )
        {
          if ( currRect.m_min[idx] < rect.m_min[idx] )
            rect.m_min[idx] = currRect.m_min[idx];
          else if ( currRect.m_max[idx] > rect.m_max[idx] )
            rect.m_max[idx] = currRect.m_max[idx]; 
        }
      }
    }
  }
public:
  int Search(std::function<bool(const ELEMTYPE*, const ELEMTYPE*)> fun, DTArray& found)
  {
    Search(this->m_root, fun,found);
    return found.size();
  }

  int Search(const ELEMTYPE* a_min, const ELEMTYPE* a_max, DTArray& data)
  {
#if _DEBUG
    for(int index=0; index<NUMDIMS; ++index)
    {
      assert(a_min[index] <= a_max[index]);
    }
#endif //_DEBUG

    Rect rect;

    memcpy(rect.m_min,a_min,NUMDIMS*sizeof(ELEMTYPE));
    memcpy(rect.m_max,a_max,NUMDIMS*sizeof(ELEMTYPE));

    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

    int foundCount = 0;
    Search(this->m_root, &rect, foundCount, data);

    return foundCount;
  }

  bool Search(Node* a_node, Rect* a_rect, int& a_foundCount, DTArray& data)
  {
    assert(a_node);
    assert(a_node->m_level >= 0);
    assert(a_rect);

    if(a_node->IsInternalNode()) // This is an internal node in the tree
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(Base::Overlap(a_rect, &a_node->m_branch[index].m_rect))
        {
          if(!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, data))
          {
            return false; // Don't continue searching
          }
        }
      }
    }
    else // This is a leaf node
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(Base::Overlap(a_rect, &a_node->m_branch[index].m_rect))
        {
          ++a_foundCount;
          data.push_back(a_node->m_branch[index].m_data);
        }
      }
    }

    return true; // Continue searching
  }

  bool Overlap(const Rect& a_rectA, const ELEMTYPE* start,const ELEMTYPE* end)
  {
    return BoxLineIntersect(NUMDIMS,a_rectA.m_min,a_rectA.m_max,start,end);
  }

  bool Search(Node* a_node, const ELEMTYPE* start,const ELEMTYPE* end, int& a_foundCount, DTArray& data)
  {
    assert(a_node);
    assert(a_node->m_level >= 0);

    if(a_node->IsInternalNode()) // This is an internal node in the tree
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(Overlap(a_node->m_branch[index].m_rect,start,end))
        {
          if(!Search(a_node->m_branch[index].m_child, start, end,a_foundCount, data))
          {
            return false; // Don't continue searching
          }
        }
      }
    }
    else // This is a leaf node
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(Overlap(a_node->m_branch[index].m_rect,start,end))
        {
          ++a_foundCount;
          data.push_back(a_node->m_branch[index].m_data);
        }
      }
    }

    return true; // Continue searching
  }

  int SearchLine(const ELEMTYPE* start, const ELEMTYPE* end, DTArray& data)
  {
    // NOTE: May want to return search result another way, perhaps returning the number of found elements here.
    int foundCount = 0;
    Search(this->m_root, start,end, foundCount, data);
    return foundCount;
  }

  bool Search(Node* a_node, std::function<bool(const ELEMTYPE*, const ELEMTYPE*)> overlap_fun, DTArray& found)
  {
    assert(a_node);
    assert(a_node->m_level >= 0);

    if(a_node->IsInternalNode()) // This is an internal node in the tree
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(overlap_fun(a_node->m_branch[index].m_rect.m_min,a_node->m_branch[index].m_rect.m_max))
        {
          if(!Search(a_node->m_branch[index].m_child, overlap_fun, found))
          {
            return false; // Don't continue searching
          }
        }
      }
    }
    else // This is a leaf node
    {
      for(int index=0; index < a_node->m_count; ++index)
      {
        if(overlap_fun(a_node->m_branch[index].m_rect.m_min,a_node->m_branch[index].m_rect.m_max))
        {
          found.push_back(a_node->m_branch[index].m_data);
        }
      }
    }

    return true; // Continue searching
  }



  template<class DistFunction>
  bool Search1nn( const ELEMTYPE *point, 
    DistFunction& fun, 
    DATATYPE &nn1Data, ELEMTYPE &minDist )
  {
    return Search1nnVLamda(point, [&fun,&point](DATATYPE val) {return fun(point,val); }, nn1Data, minDist);
  }

  template<class DistFunction>
  bool Search1nn(const ELEMTYPE* a_min, const ELEMTYPE* a_max, 
                 DistFunction& fun, DATATYPE &nn1Data, ELEMTYPE &minDist )
  {
    return Search1nnVLamda(a_min, a_max, [&fun](DATATYPE val) {return fun(val); }, nn1Data, minDist);
  }

  int Intersect( Node* node_a, Node* node_b, std::stack< NodePair > &nodePairStack, std::function<bool(const DATATYPE&, const DATATYPE&)> collectIntersects )
  {
    assert( node_a );
    assert( node_a->m_level >= 0 );

    assert( node_b );
    assert( node_b->m_level >= 0 );

    Rect rect_a;
    Rect rect_b;

    GetRect(*node_a, rect_a );
    GetRect(*node_b, rect_b );

    if ( !Base::Overlap( &rect_a, &rect_b ) )
      return 0;

    int index = 0;

    if ( !node_a->IsInternalNode() && !node_b->IsInternalNode() )
    {
      for ( int idx = 0; idx < node_a->m_count; ++idx )
      {
        for ( int jdx = 0; jdx < node_b->m_count; ++jdx )
        {
          if ( !Base::Overlap( &node_a->m_branch[idx].m_rect, &node_b->m_branch[jdx].m_rect ) )
            continue;

          if ( node_a->m_branch[idx].m_data == node_b->m_branch[jdx].m_data )
            continue;

          if(!collectIntersects( node_a->m_branch[idx].m_data, node_b->m_branch[jdx].m_data ))
            return -1; // user break
        }
      }
    }else if ( !node_a->IsInternalNode() && node_b->IsInternalNode() )
    {
      for ( index = 0; index < node_b->m_count; ++index )
      {
        Branch &currBranch = node_b->m_branch[index];

        if ( !Base::Overlap( &rect_a, &currBranch.m_rect ) )
          continue;
        nodePairStack.push(NodePair( node_a, currBranch.m_child) );
      }
    }else if ( node_a->IsInternalNode() && !node_b->IsInternalNode() )
    {
      for ( index = 0; index < node_a->m_count; ++index )
      {
        Branch &currBranch = node_a->m_branch[index];

        if ( !Base::Overlap( &currBranch.m_rect, &rect_b ) )
          continue;

        nodePairStack.push(NodePair(currBranch.m_child, node_b) );
      }
    }else
    {
      ELEMTYPEREAL volume_a = Base::CalcRectVolume( &rect_a );
      ELEMTYPEREAL volume_b = Base::CalcRectVolume( &rect_b );

      if ( volume_a <= volume_b )
      {
        for ( index = 0; index < node_b->m_count; ++index )
        {
          Branch &currBranch = node_b->m_branch[index];

          if ( !Base::Overlap( &rect_a, &currBranch.m_rect ) )
            continue;
          nodePairStack.push(NodePair( node_a, currBranch.m_child) );
        }
      }
      else
      {
        for ( index = 0; index < node_a->m_count; ++index )
        {
          Branch &currBranch = node_a->m_branch[index];

          if ( !Base::Overlap( &currBranch.m_rect, &rect_b ) )
            continue;
          nodePairStack.push(NodePair( currBranch.m_child, node_b) );
        }
      }
    }
    return 0;
  }

  void GetOverlap(const Base& atree, std::function<bool(const DATATYPE&,const DATATYPE&)> rec)
  {
    std::stack< NodePair > nodePairStack;
    try
    {
      nodePairStack.push(NodePair( this->m_root, atree.m_root) );
    }
    catch ( ... )
    {
      return;
    }

    while ( !nodePairStack.empty() )
    {
      NodePair currNodePair = nodePairStack.top();
      nodePairStack.pop();

      if(-1==Intersect( currNodePair.node_a, currNodePair.node_b, nodePairStack, rec ))
        break;//stop search
    }
  }

  void GetSelfOverlap(std::function<bool(const DATATYPE&, const DATATYPE&)> rec)
  {
    try
    {
      std::stack< NodePair > nodePairStack;
      nodePairStack.push(NodePair( this->m_root, this->m_root) );
      while ( !nodePairStack.empty() )
      {
        NodePair currNodePair = nodePairStack.top();
        nodePairStack.pop();
        if(-1==Intersect( currNodePair.node_a, currNodePair.node_b, nodePairStack, rec ))
          return;// user break;
      }
    }catch ( ... )
    {
      return;
    }
  }

  template<class T>
  struct VisitorDummy
  {
    bool operator()(T,T)
    {
      return true;
    }
  };

  template<class DistFunction,class Visitor>
  bool Search2nn( Node* a_node, const Branch &branch,DistFunction& fun,DATATYPE *nn2Data, ELEMTYPE &minDist,Visitor* viPtr )
  {
    assert( a_node );
    assert( a_node->m_level >= 0 );
    assert( a_node->m_count != 0 );

    bool success = true;
    int  index = 0;

    if ( a_node->IsInternalNode() )
    {
      std::vector< BranchInfo > branchInfos;
      branchInfos.resize( a_node->m_count );

      for ( index = 0; index < a_node->m_count; ++index )
      {
        branchInfos[index].index     = index;
        branchInfos[index].sqDistInf = SquaredDistInf(a_node->m_branch[index].m_rect, branch.m_rect );
      }

      std::sort( branchInfos.begin(), branchInfos.end() );

      if(!viPtr)
      {
        for (unsigned idx = 0; idx < branchInfos.size(); ++idx )
        {
          Branch &currBranch = a_node->m_branch[branchInfos[idx].index];

          ELEMTYPE sqDistSup = SquaredDistSup(currBranch.m_rect, branch.m_rect );

          if ( sqDistSup <= minDist )
          {
            minDist = sqDistSup;
          }
        }
      }

      for (unsigned idx = 0; idx < branchInfos.size(); ++idx )
      {
        if ( branchInfos[idx].sqDistInf < minDist )
        {
          Branch &currBranch = a_node->m_branch[branchInfos[idx].index]; 
          success = Search2nn( currBranch.m_child, branch, fun, nn2Data, minDist,viPtr );

          if ( !success )
            return false;
        }
      }
    } else 
    {
      DATATYPE currData0 = branch.m_data;
      for ( index = 0; index < a_node->m_count; ++index )
      {
        DATATYPE currData = a_node->m_branch[index].m_data;
        ELEMTYPE currDist = fun(currData,currData0 );
        if ( currDist <= minDist )
        {
          if(viPtr)
          {
            (*viPtr)(currData,currData0);
          }else
          {
            minDist = currDist;
            nn2Data[0] = currData;
            nn2Data[1] = currData0;
          }
        }
      }
    }

    return true;
  }

  template<class DistFunction,class Visitor>
  bool Search2nn( Node* a_node, Node* b_node,DistFunction& fun,DATATYPE *nn2Data, ELEMTYPE &minDist,Visitor*  viPtr )
  {
    assert( a_node && b_node);
    assert( a_node->m_level >= 0 && b_node->m_level>=0);
    assert( a_node->m_count != 0 && b_node->m_count != 0);

    try
    {
      int  index = 0;
      if (b_node->IsInternalNode())
      {
        for ( index = 0; index < b_node->m_count; ++index )
        {
          Search2nn(a_node,b_node->m_branch[index].m_child,fun,nn2Data,minDist,viPtr);
        }
      }else
      {
        for ( index = 0; index < b_node->m_count; ++index )
        {
          Search2nn(a_node,b_node->m_branch[index],fun,nn2Data,minDist,viPtr);
        }
      }  
    }catch(...)
    {
      return false;
    }
    return true;
  }

  Node* Root(){return this->m_root;}

  template<class DistFunction>
  bool Search2nn(Node* another,DistFunction& fun,DATATYPE *nn1Data, ELEMTYPE &minDist )
  {
    if ( minDist < INT_MAX )
      minDist = INT_MAX;
     
    bool success = Search2nn<DistFunction,VisitorDummy<DATATYPE> >( this->m_root, another, fun, nn1Data, minDist,NULL );
    return success;
  }

  template<class DistFunction,class Visitor>
  bool Search2nn(Node* another,DistFunction& fun, ELEMTYPE minDist,Visitor& vi )
  {
    DATATYPE nn1Data[2];
    bool success = Search2nn( this->m_root, another, fun, nn1Data, minDist,&vi );
    return success;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool Search1nnVLamda(const ELEMTYPE *point,
    std::function<double(DATATYPE i)> fun, DATATYPE &nn1Data, ELEMTYPE &minDist)const
  {
    nn1Data = 0;

    if (minDist < INT_MAX)
      minDist = INT_MAX;
    return Search1nnVLamda(this->m_root, point, fun, nn1Data, minDist);
  }

  bool Search1nnVLamda(const ELEMTYPE *point,
    std::function<double(DATATYPE i)> fun)const
  {
    DATATYPE nn1Data(0);
    ELEMTYPE minDist(INT_MAX);
    return Search1nnVLamda(this->m_root, point, fun, nn1Data, minDist);
  }

  bool Search1nnVLamda(Node* a_node, const ELEMTYPE *point,
    std::function<double(DATATYPE)> fun,
    DATATYPE &nn1Data, ELEMTYPE &minDist)const
  {
    assert(a_node);
    assert(a_node->m_level >= 0);
    assert(point);
    assert(a_node->m_count != 0);

    bool success = true;
    int  index = 0;

    if (a_node->IsInternalNode())
    {
      std::vector< BranchInfo > branchInfos;

      try
      {
        branchInfos.resize(a_node->m_count);

        for (index = 0; index < a_node->m_count; ++index)
        {
          branchInfos[index].index = index;
          branchInfos[index].sqDistInf = SquaredDistInf(a_node->m_branch[index].m_rect, point);
        }

        std::sort(branchInfos.begin(), branchInfos.end());
      }
      catch (...)
      {
        return false;
      }

      unsigned int idx = 0;

      for (idx = 0; idx < branchInfos.size(); ++idx)
      {
        Branch &currBranch = a_node->m_branch[branchInfos[idx].index];

        ELEMTYPE sqDistSup = SquaredDistSup(currBranch.m_rect, point);

        if (sqDistSup <= minDist)
        {
          minDist = sqDistSup;
          nn1Data = 0;
        }
      }

      for (idx = 0; idx < branchInfos.size(); ++idx)
      {
        if (branchInfos[idx].sqDistInf < minDist)
        {
          Branch &currBranch = a_node->m_branch[branchInfos[idx].index];
          success = Search1nnVLamda(currBranch.m_child, point, fun, nn1Data, minDist);

          if (!success)
            return false;
        }
      }
    }
    else // This is a leaf node
    {
      for (index = 0; index < a_node->m_count; ++index)
      {
        DATATYPE currData = a_node->m_branch[index].m_data;

        ELEMTYPE currDist = fun(currData);

        if (currDist <= minDist)
        {
          minDist = currDist;
          nn1Data = currData;
        }
      }
    }

    return true;
  }

  bool Search1nnVLamda(Node* a_node, const Rect &rect,
    std::function<double(DATATYPE)> fun,
    DATATYPE &nn1Data, ELEMTYPE &minDist)const
  {
    assert(a_node);
    assert(a_node->m_level >= 0);
    assert(a_node->m_count != 0);

    bool success = true;
    int  index = 0;

    if (a_node->IsInternalNode())
    {
      std::vector< BranchInfo > branchInfos;

      try
      {
        branchInfos.resize(a_node->m_count);

        for (index = 0; index < a_node->m_count; ++index)
        {
          branchInfos[index].index = index;
          branchInfos[index].sqDistInf = SquaredDistInf(a_node->m_branch[index].m_rect, rect);
        }

        std::sort(branchInfos.begin(), branchInfos.end());
      }
      catch (...)
      {
        return false;
      }

      unsigned int idx = 0;

      for (idx = 0; idx < branchInfos.size(); ++idx)
      {
        Branch &currBranch = a_node->m_branch[branchInfos[idx].index];

        ELEMTYPE sqDistSup = SquaredDistSup(currBranch.m_rect, rect);

        if (sqDistSup <= minDist)
        {
          minDist = sqDistSup;
          nn1Data = 0;
        }
      }

      for (idx = 0; idx < branchInfos.size(); ++idx)
      {
        if (branchInfos[idx].sqDistInf < minDist)
        {
          Branch &currBranch = a_node->m_branch[branchInfos[idx].index];
          success = Search1nnVLamda(currBranch.m_child, rect, fun, nn1Data, minDist);

          if (!success)
            return false;
        }
      }
    }
    else // This is a leaf node
    {
      //if ( &aSqDistToObjectCallback )
      {
        for (index = 0; index < a_node->m_count; ++index)
        {
          DATATYPE currData = a_node->m_branch[index].m_data;

          ELEMTYPE currDist = fun(currData);

          if (currDist <= minDist)
          {
            minDist = currDist;
            nn1Data = currData;
          }
        }
      }
    }

    return true;
  }


  bool Search1nnVLamda(const ELEMTYPE* a_min, const ELEMTYPE* a_max,
    std::function<double(DATATYPE)> fun, DATATYPE &nn1Data, ELEMTYPE &minDist)const
  {
#if _DEBUG
    for (int index = 0; index<NUMDIMS; ++index)
    {
      assert(a_min[index] <= a_max[index]);
    }
#endif //_DEBUG

    Rect rect;

    memcpy(rect.m_min, a_min, NUMDIMS*sizeof(ELEMTYPE));
    memcpy(rect.m_max, a_max, NUMDIMS*sizeof(ELEMTYPE));
    nn1Data = 0;

    if (minDist < INT_MAX)
      minDist = INT_MAX;

    return Search1nnVLamda(this->m_root, rect, fun, nn1Data, minDist);
  }
};

typedef RTreeExt<int,double,3> RTree3DInt;
typedef RTreeExt<int,double,2> RTree2DInt;

