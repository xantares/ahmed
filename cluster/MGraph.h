/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  MGraph.h
    \brief This class implements a Graph based on adjacency matrix.
*/
#ifndef MATRIXGRAPH_H
#define MATRIXGRAPH_H

#include "GraphNode.h"
#include "AdjMatrix.h"
#include "dijkstra.h"
#include "CS_getPathAndDist.h"

/** The class Graph implements a Graph G = (V, E).
 
 One can choose the NodeType as the type for the node-set and the
 EdgeType as the type for the edge-set. The EdgeType should be a
 numeric type, so that the algorithm of dijkstra for instance can be used with
 this class. The class holds the nodes/vertices of the graph into the array
 nodes. The edge-set is stored in a adjacency matrix.
 
 The algorithm of dijkstra computes shortest paths between a start node s and
 all other nodes of the graph.
 */
template <class NodeType, class EdgeType> class MGraph
{
public:
   /** Compute the distance between one and two with dijkstra algorithm
    and return this distance. 
   \param one start node of graph
   \param two target node
   \returns the dist between node one and node two in graph
   */
   	EdgeType getDist (const NodeType &one, const NodeType &two);

	EdgeType getPathsAndDist (const NodeType &one, const NodeType &two,
		std::list<NodeType>&, std::list<NodeType>&, std::list<NodeType>& path2);

   /** Returns the number of nodes of a Graph. (used by dijkstra)
   \returns number of nodes */
   unsigned getNumberOfNodes () const { return n_nodes; }
   unsigned getNumberOfEdges () const { return adj->getRowIdxArray()[adj->getRows ()] / 2; }
   /** Returns the array of nodes.
    \returns array of graph nodes */
   GraphNode<NodeType>** getNodes () const { return nodes; }

   /** Method checks existence of an edge between node one and node two.
    If there is an edge, the edge weight is returned in weight.
    \param one first node
    \param two second node
    \return true if edge exists else false
    */
   bool isEdge ( GraphNode<NodeType>* &one, GraphNode<NodeType>* &two ) const
   {
      unsigned row = one->getIndex (), col = two->getIndex ();
      if (row != col && (adj->isElement (row, col) >= 0)) return true;
      return false;
   }

   unsigned getEdgeWeight ( GraphNode<NodeType>* &one, GraphNode<NodeType>* &two ) const
   {
      unsigned row = one->getIndex (), col = two->getIndex ();
      if (row != col && (adj->isElement (row, col) >= 0)) 
      	return adj->getElement(row, col);
      return std::numeric_limits<unsigned>::max();
   }

   /** This is the constructor of the graph G = (V, E) = (in_nodes, in_adj).
    Initializes the attributes. It sets the entries of dist_matrix.
   \param in_nodes array of nodes (length: n_nodes = in_adj->getRows())
   \param in_adj adjacency matrix (set of edges)
   */
   MGraph ( NodeType* in_nodes, AdjMat* in_adj )
         :	n_nodes (in_adj->getRows()),
         nodes (new GraphNode<NodeType>* [n_nodes]),
         adj (in_adj)
   {
      for (unsigned k=0; k<n_nodes; k++)
         nodes[k] = new GraphNode<NodeType> (k, in_nodes[k]);
   };

   bool isGraphConnected () const
   {
      unsigned diam;
      return adj->isGraphConnected (diam);
   }

   void printGraph ( bool long_version = true )
   {
      std::cout << "(" << n_nodes << ", " << adj->getNNZ() << ") " << std::endl;
      if (long_version) {
         std::cout << "\t nodes: " << std::endl;
         for (unsigned k=0; k<n_nodes; k++) {
            if (nodes[k]->getNodeData()->isSeparator()) std::cout << "*";
            std::cout << "[" << nodes[k]->getNodeData()->getnbeg() << ",";
            std::cout << nodes[k]->getNodeData()->getnend() << ") ";
         }
         std::cout << std::endl;
         unsigned rows (adj->getRows());
         unsigned *data (adj->getDataArray());
         for (unsigned r = 0; r < rows; ++r) {
         	std::cout << r << ": " << std::flush;
            unsigned i0 (adj->getBeginRow(r)), i1 (adj->getBeginRow(r+1));
            for (unsigned k = i0; k < i1; ++k) {
               std::cout << adj->getCol(k) << "(" << data[k] << ") ";
            }
            std::cout << std::endl;
         }
      }
   }
   
   //! destructor
   ~MGraph ()
   {
   	  for (unsigned k=0; k<n_nodes; k++) delete nodes[k];
      delete [] nodes;
      delete adj;
   }

private: // attributes *********
   //! number of nodes
   unsigned n_nodes;
   //! nodes array
   GraphNode<NodeType>** nodes;
   //! edges from graph, saved in adjacency matrix
   AdjMat *adj;

private: // methods **********
   /** The default constructor is only declared to prevent the compiler to create
     a default constructor automaticaly. */
   MGraph ()
   {}

   /** The copy constructor is only declared to prevent the compiler to create
     a copy constructor automaticaly. The user can not make a copy of Graph. */
   MGraph ( const MGraph <NodeType, EdgeType>& )
   {}

   /** The assignment operator is hidden because it is private. It prevents
     the compiler to create one automaticaly. The user can not assign instances of
     Graphs. */
   MGraph<NodeType, EdgeType>& operator=
   (const MGraph<NodeType, EdgeType> &)
   {}
};

template<class NodeType, class EdgeType>
EdgeType MGraph<NodeType, EdgeType>::getDist
(const NodeType &one, const NodeType &two)
{
   if (adj->getNNZ() == 0) {
      std::cout << "\033[33m -> MGraph::getDist () error: no edges ... \033[0m" << std::flush;
      std::cout << "[" << one->getnbeg() << "," << one->getnend() << "] (";
      std::cout << one->getDepth() << "), ";
      std::cout << "[" << two->getnbeg() << "," << two->getnend() << "] (" << std::flush;
      std::cout << two->getDepth() << ") " << std::endl;
      printGraph();
      return 0;
   }

   unsigned idx1 = n_nodes+1, idx2 = n_nodes+1;
   unsigned k = 0;

   while ((idx1 == n_nodes+1 || idx2 == n_nodes+1) && k < n_nodes) {
      if (one == nodes[k]->getNodeData ()) idx1 = nodes[k]->getIndex ();
      if (two == nodes[k]->getNodeData ()) idx2 = nodes[k]->getIndex ();
      k++;
   }

   unsigned dist (dijkstra (nodes[idx2], nodes[idx1], this));
   //unsigned dist (dijkstra (idx2, idx1, adj));   

   return dist;
}

template<class NodeType, class EdgeType>
EdgeType MGraph<NodeType, EdgeType>::getPathsAndDist 
	(const NodeType &t1, const NodeType &t2, std::list<NodeType>& path0,
	std::list<NodeType>& path1, std::list<NodeType>& path2)
{
	unsigned src(n_nodes+1), target(n_nodes+1), k(0);
	EdgeType dist0(0), dist1(0);
	int *tree (new int[n_nodes]);
	unsigned *iA(adj->getRowIdxArray()), *jA(adj->getColIdxArray()) ;

	while ((src == n_nodes+1 || target == n_nodes+1) && k < n_nodes) {
		if (t1 == nodes[k]->getNodeData ()) src = nodes[k]->getIndex ();
		if (t2 == nodes[k]->getNodeData ()) target = nodes[k]->getIndex ();
		k++;
	}
	
//	if (src == n_nodes+1 || target == n_nodes+1) std::cout << "src or target not found" << std::endl;
	
	// init tree
	for (k=0; k<n_nodes; ++k) tree[k] = -1;
	tree[src] = src;
	
	dist0 = getPathAndDist (src, target, tree, n_nodes, iA, jA);
	// save the first path
	int node (target);
	path0.push_back (nodes[node]->getNodeData());
	while (tree[node] != -1 && tree[node] != node) {
		path0.push_back (nodes[tree[node]]->getNodeData()); 
		node = tree[node];
	}
	 	
	// *** second path ********** 
	// init tree
	for (k=0; k<n_nodes; ++k) tree[k] = -1;
	tree[src] = src;
	
	// do not allow to use this node for the path
	typename std::list<NodeType>::iterator it1 (++path0.begin());
	typename std::list<NodeType>::iterator it2 (-- path0.end());
	it2--;
	// std::list<NodeType> is a dependent name. The compiler cannot know that
	// std::list<NodeType>::iterator is a type
	unsigned node2avoid1 (n_nodes+1), node2avoid2 (n_nodes+1);
	k=0;
	while ((node2avoid1 == n_nodes+1 || (node2avoid2 == n_nodes+1)) && k < n_nodes) {
		if ((*it1) == nodes[k]->getNodeData ()) node2avoid1 = nodes[k]->getIndex ();
		if ((*it2) == nodes[k]->getNodeData ()) node2avoid2 = nodes[k]->getIndex ();
		k++;
	}
	tree[node2avoid1] = n_nodes+1;
	
	dist1 = getPathAndDist (src, target, tree, n_nodes, iA, jA);
	// save the second path
	node = target;
	path1.push_back (nodes[node]->getNodeData());
	while (tree[node] != -1 && tree[node] != node) {
		path1.push_back (nodes[tree[node]]->getNodeData()); 
		node = tree[node];
	}
	if (tree[node] == -1) path1.clear ();

	// *** third path **********
	// init tree
	for (k=0; k<n_nodes; ++k) tree[k] = -1;
	tree[src] = src;
	tree[node2avoid2] = n_nodes+1;
	
	dist1 = getPathAndDist (src, target, tree, n_nodes, iA, jA);
	// save the third path
	node = target;
	path2.push_back (nodes[node]->getNodeData());
	while (tree[node] != -1 && tree[node] != node) {
		path2.push_back (nodes[tree[node]]->getNodeData()); 
		node = tree[node];
	}
	if (tree[node] == -1) path2.clear ();

	delete [] tree;
	
	return dist0;
}

#endif
