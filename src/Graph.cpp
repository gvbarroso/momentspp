/*
 * Authors: Gustavo numVertices_. Barroso
 * Created: 27/03/2023
 * Last modified: 27/03/2023
 *
 */

#include "Graph.hpp"

bool Graph::isReachable(int source, int dest)
{
  if(source == dest)
    return true;

  // Mark all the vertices as not visited
  bool *visited = new bool[numVertices_];
  for(int i = 0; i < numVertices_; i++)
    visited[i] = false;

  std::list<int> queue;
  // Mark the current node as visited and enqueue it
  visited[source] = true;
  queue.push_back(source);

  // it will be used to get all adjacent vertices of a vertex
  std::list<int>::iterator it;

  while(!queue.empty())
  {
    // Dequeue a vertex from queue and print it
    source = queue.front();
    queue.pop_front();

    // Get all adjacent vertices of the dequeued vertex source
    // If a adjacent has not been visited, then mark it visited
    // and enqueue it
    for(it = std::begin(adj_[source]); it != std::end(adj_[source]); ++it)
    {
      // If this adjacent node is the destination node, then
      if (*it == dest)
        return true;

      // Else, continue to do BFS
      if (!visited[*it])
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
