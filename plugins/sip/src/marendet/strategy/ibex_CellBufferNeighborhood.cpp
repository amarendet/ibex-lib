/* ============================================================================
 * I B E X - ibex_CellBufferNeighborhood.cpp
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Marendet, Gilles Chabert
 * Created     : May 4, 2018
 * ---------------------------------------------------------------------------- */
 
#include "ibex_CellBufferNeighborhood.h"
#include "ibex_utils.h"

#include <limits>
#include <stack>
#include <utility>

using namespace std;

namespace ibex {

CellBufferNeighborhood::CellBufferNeighborhood(CellBufferNeighborhood::Heuristic heuristic)
    : start_(IntervalVector::empty(1))
    , goal_(IntervalVector::empty(1))
    , heuristic_(heuristic)
{
    flush();
}

CellBufferNeighborhood::~CellBufferNeighborhood() {
    flush();
}

void CellBufferNeighborhood::init(const IntervalVector& start, const IntervalVector& goal) {
    start_ = start;
    goal_ = goal;
    start_node_ = new Node{new Cell(start_), true};
    goal_node_ = new Node{new Cell(goal_), true};
    insert_into_graph(start_node_);
    insert_into_graph(goal_node_);
}

void CellBufferNeighborhood::flush() {
    for(auto pair : graph_list_) {
        const Node* n = pair.first;
        delete n->cell;
        delete n;
    }
    update_top_ = true;
    graph_list_.clear();
    size_ = 0;
    init(start_, goal_);
}

unsigned int CellBufferNeighborhood::size() const {
    if(empty()) {
        return 0;
    } else {
        return size_;
    }
}

bool CellBufferNeighborhood::empty() const {
    update_top();
    return top_ == nullptr;
}

void CellBufferNeighborhood::push(Cell* cell) {
    const Node* new_node = new Node{cell, false};   
    insert_into_graph(new_node);
    size_ += 1;
}

Cell* CellBufferNeighborhood::pop() {
    update_top();
    if(top_ == nullptr) {
        return nullptr;
    } else {
        Cell* cell = top_->cell;
        delete_from_graph(top_);
        size_ -= 1;
        return cell;
    }
}
Cell* CellBufferNeighborhood::top() const {
    update_top();
    //cout << top_->cell->box.mid() << endl;
    return top_->cell;
}

void CellBufferNeighborhood::push_inner(Cell* cell) {
    // Idea : don't push cell already fully contained in other inner cells
    const Node* new_node = new Node{cell, true};   
    insert_into_graph(new_node);
}

bool CellBufferNeighborhood::feasible_path_found() {
    return feasible_path_found_;
}

vector<IntervalVector> CellBufferNeighborhood::path() const {
    return shortest_path_iv_;
}

void CellBufferNeighborhood::insert_into_graph(const Node* node) {
    update_top_ = true;
    map_edge_by_node neighbors;
    for(auto& pair : graph_list_) {
        const Node* other = pair.first;
        map_edge_by_node& other_neighbors = pair.second;
        if(node->cell->box.intersects(other->cell->box)) {
            double w = weight(node, other);
            neighbors.emplace(other, Edge{node, other, w});
            other_neighbors.emplace(node, Edge{other, node, w});
        }
    }
    graph_list_.emplace(node, neighbors);
}

void CellBufferNeighborhood::delete_from_graph(const Node* node) {
    update_top_ = true;
    map_edge_by_node& neighbors = graph_list_[node];
    for(auto& pair : neighbors) {
        const Edge& e = pair.second;
        map_edge_by_node& other_neighbors = graph_list_[e.to];
        other_neighbors.erase(node);
    }
    graph_list_.erase(node);
}

double CellBufferNeighborhood::weight(const Node* from, const Node* to) const {
    double w = norm(from->cell->box.mid() - to->cell->box.mid());
    //  Vector v1 = from->cell->box.mid().subvector(0,1);
    //  Vector v2 = to->cell->box.mid().subvector(0,1);
    //  double w = (v1-v2)*(v1-v2);
    // cout << w << endl;
    if(from->inner || to->inner) {
        return w;
    } else {
        return 100*w;
    }
}

void CellBufferNeighborhood::update_top() const {
    if(!update_top_) {
        return;
    }
    if(alternate_)
        update_path_a_star(goal_node_, start_node_);
    else
        update_path_a_star(start_node_, goal_node_);
    alternate_ = !alternate_;
    top_ = first_unknown_node_on_path();
    update_top_ = false;
}

const CellBufferNeighborhood::Node* CellBufferNeighborhood::first_unknown_node_on_path() const {
    shortest_path_iv_.clear();
    for(const Node* n : shortest_path_) {
        if(!n->inner) {
            return n;
        }
        shortest_path_iv_.emplace(shortest_path_iv_.begin(), n->cell->box);
    }
    feasible_path_found_ = !impossible_;
    return nullptr;
}

void CellBufferNeighborhood::update_path_a_star(const Node* start_node, const Node* goal_node) const {

    map<const Node*, double> priority; // the lowest is first;
    map<const Node*, double> costs;
    map<const Node*, const Node*> came_from;
    auto cmp = [&](const Node* n1, const Node* n2) {
        return priority[n1] < priority[n2];
    };
    auto heuristic = [&](const Node* n) {
        if(heuristic_ == Heuristic::A_STAR_DISTANCE) {
            return norm(n->cell->box.mid() - goal_node->cell->box.mid());
        } else {
            return 0.0;
        }
    };
    const double infinity = std::numeric_limits<double>::max();

    set<const Node*, decltype(cmp)> open(cmp);
    set<const Node*> closed;

    priority[start_node] = 0;
    costs[start_node] = 0;
    came_from[start_node] = nullptr;
    open.insert(start_node);

    while(!open.empty() && *open.begin() != goal_node) {
        const Node* current = *open.begin();
        open.erase(*open.begin());
        closed.insert(current);
        const map_edge_by_node& neighbors = graph_list_.at(current);
        for(auto const& pair : neighbors) {
            const Edge& e = pair.second;
            double cost = costs[current] + e.weight;
            bool in_open = set_contains(open, e.to);
            bool in_closed = set_contains(closed, e.to);
            // auto it_open = find(open.begin(), open.end(), e.to);
            // auto it_closed = find(closed.begin(), closed.end(), e.to);
            bool is_better_path = cost < map_get_with_default(costs, e.to, infinity);
            if(in_open && is_better_path) {
                open.erase(e.to);
                in_open = false;    
            }
            if(in_closed && is_better_path) { // Happens only with inadmissible heuristics
                closed.erase(e.to);
                in_closed = false;
            }
            // if(it_open != open.end() && is_better_path) {
            //     open.erase(it_open);
            //     it_open = open.end();
            // }
            // if(it_closed != closed.end() && is_better_path) { // Happens only with inadmissible heuristics
            //     closed.erase(it_closed);
            //     it_closed = closed.end();
            // }
            if(!in_open && !in_closed) {
                costs[e.to] = cost;
                priority[e.to] = cost + heuristic(e.to);
                open.insert(e.to);
                came_from[e.to] = current;
            }
        }
    }
    if(*open.begin() != goal_node) {
        // No path is possible
        ibex_warning("CellBufferNeighborhood: no potential path found");
        goal_node = nullptr;
        impossible_ = true;
    }
    reconstruct_path(came_from, goal_node);
}
void CellBufferNeighborhood::reconstruct_path(const std::map<const Node*, const Node*>& came_from, const Node* last) const {
    shortest_path_.clear();
    const Node* current = last;
    while(current != nullptr) {
        shortest_path_.emplace(shortest_path_.begin(), current);
        current = came_from.at(current);
    }
}

} // end namespace ibex
