/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef GRAPHNODE_H
#define GRAPHNODE_H

template <class NodeDataType>
class GraphNode
{
  unsigned idx;
  NodeDataType node_data;

 public:
  GraphNode(unsigned i) : idx(i) {}
  GraphNode(unsigned i, NodeDataType data) : idx(i), node_data(data) {}
  unsigned getIndex() const { return idx; }
  NodeDataType& getNodeData() { return node_data; }
};

#endif
