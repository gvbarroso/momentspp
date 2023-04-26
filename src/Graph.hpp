/*
 * Authors: Gustavo V. Barroso
 * Created: 27/03/2023
 * Last modified: 26/04/2023

 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#define NIL -1

#include <iostream>
#include <list>
#include <stack>

// based on https://www.geeksforgeeks.org/tarjan-algorithm-find-strongly-connected-components/
template<typename ValueType>
class Graph
{

private:
  size_t numVertices_;
  std::list<ValueType>* adj_; // A dynamic array of adjacency lists

public:
  Graph(size_t numVertices):
  numVertices_(numVertices),
  adj_()
  {
    adj_ = new std::list<ValueType>[numVertices];
  }

  Graph(const Graph&) = delete;
  Graph(Graph&&) = delete;
  Graph& operator=(const Graph&) = delete;
  Graph& operator=(Graph&&) = delete;

  ~Graph()
  {
    delete [] adj_;
  }

public:
  size_t getNumVertices()
  {
    return numVertices_;
  }

  void addEdge(size_t v, ValueType w)
  {
    adj_[v].push_back(w);
  }

  bool isReachable(ValueType source, ValueType destination);

};

#endif
