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

	CellBufferNeighborhood(Heuristic heuristic);
	virtual ~CellBufferNeighborhood();

	void init(const IntervalVector& start, const IntervalVector& goal);
	void flush();
	unsigned int size() const;
	bool empty() const;
	void push(Cell* cell);
	Cell* pop();
	Cell* top() const;
	void push_inner(Cell* cell);
	bool feasible_path_found();
	std::vector<IntervalVector> path() const;

private:
	struct Node {
		Cell* cell;
		bool inner;
	};

	struct Edge {
		const Node* from;
		const Node* to;
		double weight;
	};

	IntervalVector start_;
	IntervalVector goal_;
	const Node* start_node_;
	const Node* goal_node_;
	Heuristic heuristic_;
	int size_ = 0;

	mutable bool update_top_ = true;
	mutable const Node* top_ = nullptr;
	mutable std::vector<const Node*> shortest_path_;
	mutable std::vector<IntervalVector> shortest_path_iv_;
	mutable bool feasible_path_found_ = false;
	mutable bool impossible_ = false;
	mutable bool alternate_ = false;

	typedef std::map<const Node*, Edge> map_edge_by_node;
	std::map<const Node*, map_edge_by_node> graph_list_;

	void update_top() const;
	void insert_into_graph(const Node* node);
	void delete_from_graph(const Node* node);
	double weight(const Node* from, const Node* to) const;
	void update_path_a_star(const Node* start_node, const Node* goal_node) const;
	void reconstruct_path(const std::map<const Node*, const Node*>& came_from, const Node* last) const;
	const Node* first_unknown_node_on_path() const;
	
};

} // end namespace ibex

#endif // __SIP_IBEX_CELLBUFFERNEIGHBORHOOD_H__
