/* ============================================================================
 * I B E X - ibex_CellBufferNeighborhood.h
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Marendet, Gilles Chabert
 * Created     : May 4, 2018
 * ---------------------------------------------------------------------------- */
 
#ifndef __SIP_IBEX_CELLBUFFERNEIGHBORHOOD_H__
#define __SIP_IBEX_CELLBUFFERNEIGHBORHOOD_H__

#include "ibex_Cell.h"
#include "ibex_CellBuffer.h"
#include "ibex_IntervalVector.h"
#include "ibex_Vector.h"
#include "ibex_CellList.h"

#include <map>
#include <set>
#include <vector>

namespace ibex {

class CellBufferNeighborhood: public CellBuffer {
public:
	enum Heuristic { DIJKSTRA, A_STAR_DISTANCE };

	CellBufferNeighborhood(IntervalVector start, IntervalVector goal, Heuristic heuristic);
	virtual ~CellBufferNeighborhood();

	void flush();
	unsigned int size() const;
	bool empty() const;
	void push(Cell* cell);
	Cell* pop();
	Cell* top() const;
	void pushInner(Cell* cell);

//private:

	mutable std::vector<IntervalVector> pathFound;
	mutable bool isPathFound = false;

	class GraphNode;
	class Edge {
	public:
		GraphNode* node = nullptr;
		double weight = 0;
		Edge() {}
		Edge(GraphNode* node, double weight) : node(node), weight(weight) {}
		Edge(const Edge& other) : node(other.node), weight(other.weight) {}
		bool operator <(const Edge& other) const {
			return node < other.node;
		}
		bool operator==(const Edge& other) const {
			return node == other.node;
		}
	};
	class GraphNode {
	public:
		enum CellType {INNER, UNKNOWN};
		Cell* cell;
		Vector mid;
		std::set<Edge> neighborsWeight;
		CellType type;
		GraphNode(Cell* cell, CellType type) : cell(cell), mid(cell->box.mid()), type(type) {
		}

		std::set<Cell*> connectedComponent() const;
	};
	double heuristic_distance(const GraphNode& v1, const GraphNode& v2) const;
	double distance(const GraphNode& n1, const GraphNode& n2) const;
	Edge topGraphNode() const;
	std::vector<Edge> reconstructPath(const std::map<GraphNode*, Edge>& cameFrom, const Edge& current) const;
	std::vector<Edge> shortestPath(GraphNode* start, GraphNode* goal) const;
	std::set<GraphNode*> stack_;
	IntervalVector start_;
	GraphNode* start_node_ = nullptr;
	IntervalVector goal_;
	GraphNode* goal_node_ = nullptr;
	Heuristic heuristic_;
	mutable bool init_push_phase_ = true;
	mutable Edge last_top_{nullptr, 0.0};
	
	mutable CellList cell_list_start;
	mutable CellList cell_list_goal;

};

} // end namespace ibex

#endif // __SIP_IBEX_CELLBUFFERNEIGHBORHOOD_H__
