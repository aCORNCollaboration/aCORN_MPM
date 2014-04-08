#ifndef DIRECTEDGRAPH_HH
#define DIRECTEDGRAPH_HH

#include <set>
#include <vector>

/// Base class for a node on a directed graph
class DirectedGraphNode {
public:
	/// constructor
	DirectedGraphNode() {}
	/// destructor
	virtual ~DirectedGraphNode() {}
	
	std::set<DirectedGraphNode*> links;	//< other connected nodes
};

/// Base class for a directed graph with an enumerated list of nodes
class DirectedGraph {
public:
	/// constructor
	DirectedGraph() {}
	/// destructor
	virtual ~DirectedGraph() {}
	
	std::vector<DirectedGraphNode*> nodes;	//< enumerated nodes in graph
};

#endif
