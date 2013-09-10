/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <limits>
#include "MGraph.h"
#include "GraphNode.h"

template <class NodeData, class EdgeWeightType> class MGraph;

template <class NodeData, class T>
T dijkstra (GraphNode<NodeData> *s, GraphNode<NodeData> *t,
	    MGraph<NodeData, T> *graph)
{   
  GraphNode<NodeData>* act_node(s); // init act_node with start node s

  // create arrays to store costs and the previous node
  unsigned n_nodes = graph->getNumberOfNodes ();
  bool* removed = new bool [n_nodes];
  T* costs = new T [n_nodes];
  GraphNode<NodeData>** pre  = new GraphNode<NodeData>* [n_nodes];
  GraphNode<NodeData>** nodes = graph->getNodes ();
  // init costs and pre
  for (unsigned k = 0; k < n_nodes; k++) {
    costs[k] = std::numeric_limits<T>::max ();
    pre[k] = s;
    removed[k] = false;
  }
  
  unsigned k = 0;
  while (s->getIndex () != nodes[k]->getIndex()) k++;
  removed[k] = true;
  costs[k] = 0;
  unsigned min_idx = k;
  
  for (unsigned k = 0; k < n_nodes; k++) {
    // compute costs for step k
    T act_node_dist (costs[min_idx]);
    for (unsigned j = 0; j < n_nodes; j++) {
      if (! removed[j]) {
	if (graph->isEdge (act_node, nodes[j] )) {
	  unsigned dist = graph->getEdgeWeight (act_node, nodes[j]);
	  if (costs[j] > act_node_dist + dist) {
	    costs[j] = act_node_dist + dist;
	    pre[j] = act_node;
	  }
	}
      }
    }
    
    // compute min from costs vector
    T min_costs; // init
    unsigned i = 0;
    while (removed[i] && i < n_nodes) i++;
    min_idx = i;
    min_costs = costs[i];
    for (unsigned j = 0; j < n_nodes; j++ ) {
      if (! removed[j]) {
	if (costs[j] < min_costs) {
	  min_idx = j;
	  min_costs = costs[j];
	}
      }
    }

    // remove
    removed[min_idx] = true;
    if (t->getIndex () == nodes[min_idx]->getIndex()) k = n_nodes;
    
    act_node = nodes[min_idx];
  }
  T dist = costs[min_idx];
  delete [] pre; delete [] costs; delete [] removed; 
  return dist;
}

#endif
