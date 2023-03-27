/*
 * Authors: Gustavo V. Barroso
 * Created: 27/03/2023
 * Last modified: 27/03/2023

 */

#ifndef _GRAPH_H_
#define _GRAPH_H_

#define NIL -1

#include <iostream>
#include <list>
#include <stack>

// based on https://www.geeksforgeeks.org/tarjan-algorithm-find-strongly-connected-components/
class Graph
{

private:
  int numVertices_;
  std::list<int>* adj_; // A dynamic array of adjacency lists

public:
  Graph(int numVertices):
  numVertices_(numVertices),
  adj_()
  {
    adj_ = new std::list<int>[numVertices];
  }

  Graph(const Graph&) = delete;
  Graph(Graph&&) = delete;
  Graph& operator=(const Graph&) = delete;
  Graph& operator=(Graph&&) = delete;

  ~Graph() = default;

public:
  int getNumVertices()
  {
    return numVertices_;
  }

  void addEdge(int v, int w)
  {
    adj_[v].push_back(w);
  }

  bool isReachable(int source, int destination);

};

#endif
