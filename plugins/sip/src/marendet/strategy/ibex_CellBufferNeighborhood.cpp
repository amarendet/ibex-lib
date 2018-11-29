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

#include <limits>
#include <stack>
#include <utility>

namespace ibex {

CellBufferNeighborhood::CellBufferNeighborhood(IntervalVector start, IntervalVector goal, CellBufferNeighborhood::Heuristic heuristic)
    : start_(start)
    , goal_(goal)
    , heuristic_(heuristic)
{
    cell_list_start.push(new Cell(start_));
    cell_list_goal.push(new Cell(goal_));
}

CellBufferNeighborhood::~CellBufferNeighborhood()
{
    for (GraphNode* node : stack_) {
        delete node;
    }
}

void CellBufferNeighborhood::flush()
{
    for (GraphNode* node : stack_) {
        delete node->cell;
        delete node;
    }
    cell_list_start.flush();
    cell_list_goal.flush();
    cell_list_start.push(new Cell(start_));
    cell_list_goal.push(new Cell(goal_));
    stack_.clear();
}

unsigned int CellBufferNeighborhood::size() const
{
    return stack_.size();
}

bool CellBufferNeighborhood::empty() const
{
    return stack_.empty();
}

void CellBufferNeighborhood::push(Cell* cell)
{
    //std::cout << "push: " << cell->box << std::endl;
    if(!init_push_phase_ && !stack_.empty()) {
        if(start_node_ == nullptr) {
            if(start_.intersects(cell->box)) {
                cell_list_start.push(cell);
            }
            return;
        }
        if(goal_node_ == nullptr) {
            if(goal_.intersects(cell->box)) {
                cell_list_goal.push(cell);
            }
            return;
        }
    }
    GraphNode* newNode = new GraphNode(cell, GraphNode::UNKNOWN);
    for (GraphNode* node : stack_) {
        if (node->cell->box.intersects(cell->box)) {
            double weight = distance(*newNode, *node);
            //double weight = node->type == GraphNode::UNKNOWN ? 1 : 0;
            newNode->neighborsWeight.insert(Edge(node, weight));
            node->neighborsWeight.insert(Edge(newNode, weight));
        }
    }
    stack_.emplace(newNode);
    /*if (cell->box.contains(start_)) {
        start_node_ = newNode;
    }
    if (cell->box.contains(goal_)) {
        goal_node_ = newNode;
    }*/
}

void CellBufferNeighborhood::pushInner(Cell* cell)
{
    //std::cout << "inner:" << cell->box << std::endl;

    GraphNode* newNode = new GraphNode(cell, GraphNode::INNER);
    if((start_node_ == nullptr && start_.intersects(cell->box))
        || (goal_node_ == nullptr && goal_.intersects(cell->box))
        || (start_node_ != nullptr && goal_node_ != nullptr)) {
        for (GraphNode* node : stack_) {
            if (node->cell->box.intersects(cell->box)) {
                double weight = distance(*newNode, *node);
                //double weight = node->type == GraphNode::UNKNOWN ? 0 : 0;
                newNode->neighborsWeight.insert(Edge(node, weight));
                node->neighborsWeight.insert(Edge(newNode, weight));
            }
        }
    }

    if (start_node_ == nullptr) {
        if(start_.intersects(cell->box)) {
            start_node_ = newNode;
            std::cout << "starting box: " << cell->box << std::endl;
            stack_.emplace(newNode);
        }
        return;
    }
    if (goal_node_ == nullptr) {
        if(goal_.intersects(cell->box)) {
            goal_node_ = newNode;
            std::cout << "goal box: " << cell->box << std::endl;
            stack_.emplace(newNode);
        }
        return;
    }
    
    stack_.emplace(newNode);
    
}

Cell* CellBufferNeighborhood::pop()
{
    //GraphNode* top = stack_.back();
    //stack_.pop_back();
    //GraphNode* top = topGraphNode();
    Edge top = last_top_;
    if (top.node == nullptr)
        return nullptr;
    for (auto neighbor : top.node->neighborsWeight) {
        neighbor.node->neighborsWeight.erase(top);
    }
    stack_.erase(top.node);
    Cell* cell = top.node->cell;
    delete top.node;
    return cell;
}

Cell* CellBufferNeighborhood::top() const
{
    init_push_phase_ = false;
    if(start_node_ == nullptr) {
        return cell_list_start.pop();
    } else if(goal_node_ == nullptr) {
        return cell_list_goal.pop();
    }
    Edge top = topGraphNode();
    last_top_ = top;
    if (top.node == nullptr) {
        return nullptr;
    } else {
        //std::cout << "top: " << top.node->cell->box << std::endl;
        return top.node->cell;
    }
    //return stack_.back()->cell;
}

CellBufferNeighborhood::Edge CellBufferNeighborhood::topGraphNode() const
{
    // We must first prove that the starting and goal points are feasible
    /*if (start_node_->type == GraphNode::UNKNOWN) {
        return Edge(start_node_, 0.0);
    } else if (goal_node_->type == GraphNode::UNKNOWN) {
        return Edge(goal_node_, 0.0);
    }*/
    
    std::vector<Edge> path = shortestPath(start_node_, goal_node_);
    pathFound.clear();
    for (const Edge& n : path) {
        if (n.node->type == GraphNode::UNKNOWN) {
            isPathFound = false;
            return n;
        } else {
            pathFound.emplace_back(n.node->cell->box);
        }
    }
    // PATH FOUND
    /*pathFound.clear();
    for (const Edge& n : path) {
        pathFound.emplace_back(n.node->cell->box);
    }*/
    isPathFound = true;
    return Edge{ nullptr, 0.0 };
    //return stack_.back()->cell;
}

std::set<Cell*> CellBufferNeighborhood::GraphNode::connectedComponent() const
{
    std::set<Cell*> component;
    std::stack<Edge> toExplore;
    while (!toExplore.empty()) {
        Edge edge = toExplore.top();
        toExplore.pop();
        for (Edge neighbor : edge.node->neighborsWeight) {
            if (component.find(neighbor.node->cell) == component.end()) {
                // not yet in component
                toExplore.push(neighbor);
            }
        }
        component.insert(edge.node->cell);
    }
    return component;
}

double CellBufferNeighborhood::heuristic_distance(const GraphNode& v1, const GraphNode& v2) const
{
    if (heuristic_ == DIJKSTRA)
        return 0;
    else if (heuristic_ == A_STAR_DISTANCE) {
        /*IntervalVector vv1 = v1.cell->box.subvector(0,1);
        IntervalVector vv2 = v2.cell->box.subvector(0,1);
        return 5*hadamard_product((vv1-vv2),(vv1-vv2)).mig().min();*/
        return norm(v1.mid.subvector(0,1) - v2.mid.subvector(0,1));
    } else
        return 0;
}

double CellBufferNeighborhood::distance(const GraphNode& v1, const GraphNode& v2) const
{
    double w = norm(v1.mid.subvector(0,1)-v2.mid.subvector(0,1));
    if(v1.type == GraphNode::UNKNOWN || v2.type == GraphNode::UNKNOWN) {
        return 2*w;
    } else {
        return 1*w;
    }
}

std::vector<CellBufferNeighborhood::Edge> CellBufferNeighborhood::reconstructPath(
    const std::map<GraphNode*, Edge>& cameFrom, const Edge& current) const
{
    std::vector<Edge> path;
    path.emplace_back(current);
    auto it = cameFrom.find(current.node);
    while (it != cameFrom.end()) {
        path.emplace(path.end(), it->second);
        it = cameFrom.find(it->second.node);
    }
    return path;
}

std::vector<CellBufferNeighborhood::Edge> CellBufferNeighborhood::shortestPath(GraphNode* start,
    GraphNode* goal) const
{
    const double infinity = 100000;
    std::set<GraphNode*> closedSet;
    std::set<Edge> openSet;
    Edge start_edge(start, 0.0);
    openSet.emplace(start_edge);
    std::map<GraphNode*, Edge> cameFrom;
    std::map<Edge, double> gScore;
    gScore.emplace(std::make_pair(start_edge, 0.0));
    std::map<Edge, double> fScore;
    fScore.emplace(std::make_pair(start_edge, heuristic_distance(*start, *goal)));

    auto defaultMapGet = [&](const std::map<Edge, double>& map, const Edge& value) {
        const auto it = map.find(value);
        if (it == map.end()) {
            return infinity;
        } else {
            return it->second;
        }
    };

    while (!openSet.empty()) {
        double min_f_score = infinity;
        Edge current = *openSet.begin();
        for(Edge e : openSet) {
            double fScore_tmp = defaultMapGet(fScore, e);
            if(fScore_tmp < min_f_score) {
                current = e;
                min_f_score = fScore_tmp;
            }
        }
        double fScore_value = defaultMapGet(fScore, current);
        for (const Edge& n : openSet) {
            double tmp = defaultMapGet(fScore, n);
            if (tmp < fScore_value) {
                fScore_value = tmp;
                current = n;
            }
        }
        /*GraphNode* current = *std::min(openSet.begin(), openSet.end(), [&,fScore](GraphNode* n1, GraphNode* n2) {
		 return defaultMapGet(fScore, n1) < defaultMapGet(fScore, n2);
		 });*/
        if (current.node == goal) {
            return reconstructPath(cameFrom, current);
        }
        openSet.erase(current);
        closedSet.emplace(current.node);
        for (const Edge& neighbor : current.node->neighborsWeight) {
            if (closedSet.find(neighbor.node) != closedSet.end()) {
                continue;
            }
            openSet.emplace(neighbor);

            double tentative_gScore = defaultMapGet(gScore, current) + current.weight;
            if (tentative_gScore >= defaultMapGet(gScore, neighbor)) {
                continue;
            }
            cameFrom[neighbor.node] = current;
            gScore[neighbor] = tentative_gScore;
            fScore[neighbor] = tentative_gScore + heuristic_distance(*neighbor.node, *goal);

        }
    }
    return std::vector<Edge>();
}

} // end namespace ibex
