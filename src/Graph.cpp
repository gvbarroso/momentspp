/*
 * Authors: Gustavo V. Barroso
 * Created: 27/03/2023
 * Last modified: 11/04/2023
 *
 */

#include "Graph.hpp"

bool Graph::isReachable(int source, int dest)
{
  if(source == dest)
    return true;

  // mark all the vertices as not visited
  bool *visited = new bool[numVertices_];
  for(int i = 0; i < numVertices_; i++)
    visited[i] = false;

  std::list<int> queue;
  // mark the current node as visited and enqueue it
  visited[source] = true;
  queue.push_back(source);

  // it will be used to get all adjacent vertices of a vertex
  std::list<int>::iterator it;

  while(!queue.empty())
  {
    // dequeue a vertex from queue and print it
    source = queue.front();
    queue.pop_front();

    // get all adjacent vertices of the dequeued vertex source
    // if a adjacent has not been visited, then mark it visited
    // and enqueue it
    for(it = std::begin(adj_[source]); it != std::end(adj_[source]); ++it)
    {
      // if this adjacent node is the destination node, then
      if(*it == dest)
      {
        delete [] visited;
        return true;
      }

      // else, continue to do BFS
      if(!visited[*it])
      {
        visited[*it] = true;
        queue.push_back(*it);
      }
    }
  }

  delete [] visited;

  // If BFS is complete without visiting dest
  return false;
}
