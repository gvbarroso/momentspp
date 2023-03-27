/*
 * Authors: Gustavo numVertices_. Barroso
 * Created: 27/03/2023
 * Last modified: 27/03/2023
 *
 */

#include "Graph.hpp"

void Graph::sccUtil_(int u, int disc[], int low[], std::stack<int>* st, bool stackMember[])
{
  // A static variable is used for simplicity, we can
  // avoid use of static variable by passing a pointer.
  static int time = 0;

  // Initialize discovery time and low value
  disc[u] = low[u] = ++time;
  st->push(u);
  stackMember[u] = true;

  // Go through all vertices adj_acent to this
  std::list<int>::iterator i;
  for(auto it = std::begin(adj_[u]); it != std::end(adj_[u]); ++it)
  {
    int v = *it; // v is current adj_acent of 'u'

    // If v is not visited yet, then recur for it
    if(disc[v] == -1)
    {
      sccUtil_(v, disc, low, st, stackMember);

      // Check if the subtree rooted with 'v' has a
      // connection to one of the ancestors of 'u'
      // Case 1 (per above discussion on Disc and Low value)
      low[u] = std::min(low[u], low[v]);
    }

    // Update low value of 'u' only of 'v' is still in
    // stack (i.e. it's a back edge, not cross edge).
    // Case 2 (per above discussion on Disc and Low value)
    else if(stackMember[v] == true)
      low[u] = std::min(low[u], disc[v]);


    // head node found, pop the stack and print an scc
    int w = 0; // To store stack extracted vertices
    if(low[u] == disc[u])
    {
      while(st->top() != u)
      {
        w = (int)st->top();
        std::cout << w << " ";
        stackMember[w] = false;
        st->pop();
      }
    }

    w = (int)st->top();
    std::cout << w << "\n";
    stackMember[w] = false;
    st->pop();
  }
}

// DFS traversal
void Graph::scc()
{
  int* disc = new int[numVertices_];
  int* low = new int[numVertices_];
  bool* stackMember = new bool[numVertices_];
  std::stack<int>* st = new std::stack<int>();

  // Initialize disc and low, and stackMember arrays
  for(int i = 0; i < numVertices_; i++)
  {
    disc[i] = NIL;
    low[i] = NIL;
    stackMember[i] = false;
  }

  // Call the recursive helper function to find strongly
  // connected components in DFS tree with vertex 'i'
  for(int i = 0; i < numVertices_; i++)
    if(disc[i] == NIL)
      sccUtil_(i, disc, low, st, stackMember);
}
