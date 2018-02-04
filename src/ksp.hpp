/*
 * ksp.hpp
 *
 *  Created on: 16 Oct 2017
 *      Author: mahmoud
 */

#ifndef KSP_HPP_
#define KSP_HPP_

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "helpers.hpp"

using namespace boost;
using namespace std;

template<class Graph>
size_t shortestPath(const Graph &g, Vertex<Graph> s, Vertex<Graph> t,
		vector<Vertex<Graph>> &parentMap) {
	vector<int> distMap(num_vertices(g));
	Vertex<Graph> s_ = vertex(s, g);
	// http://www.boost.org/doc/libs/1_65_1/libs/graph/doc/dijkstra_shortest_paths.html
	dijkstra_shortest_paths(g, s_,
			predecessor_map(
					make_iterator_property_map(parentMap.begin(),
							get(vertex_index, g))).distance_map(
					make_iterator_property_map(distMap.begin(),
							get(vertex_index, g))));

	return distMap[t];
}

template<class Graph>
size_t shortestPath(Graph &g, Vertex<Graph> s, Vertex<Graph> t) {
	vector<Vertex<Graph>> p(num_vertices(g));

	return shortestPath(g, s, t, p);
}

/**
 * a k shortest path solution
 * in order to get the kth SP, remove an edge from previous SPs
 * try all remove combinations
 */
template<class Graph>
size_t bestNextEdges(Graph &g, const int d, const Vertex<Graph> s, const Vertex<Graph> t,
		vector<Edge<Graph>> &result) {

	vector<Edge<Graph>> bestCandidates; // best edge selection from recursive calls

	vector<Vertex<Graph>> parent(num_vertices(g));
	size_t bestSP = INF;
	size_t len = shortestPath(g, s, t, parent);
//	printPath(getPath(parent, s, t));
	if (len == INF) {
		return len;
	}

	Edge<Graph> bestEdge; //  best edge in this call

	if (d < 1) {
//		mylog << " \nbuttom level reached, bestSP=" << len;
		return len;

	} else { // try all edge deletes
		for (Vertex<Graph> v = t; v != s; v = parent[v]) { // for each path edge from t to s
			auto e = edge(parent[v], v, g).first;
			auto ew = get(edge_weight, g, e);
//			mylog << "\nhide edge " << edgeToStr(parent[v], v, g);
			put(edge_weight, g, e, INF); // hide the edge

//			mylog << " \n  calling bestNextEdges(),d-1=" << d - 1;

			len = bestNextEdges(g, d - 1, s, t, result); // SP not using e and the next best edges selected recursively (i.e candidates)

//			mylog << " \n  bestNextEdges() returned with len=" << len;
//			mylog << "\nuhide edge " << edgeToStr(parent[v], v, g);
			put(edge_weight, g, e, ew); // unhide the edge
//			mylog << "\nspLen=" << len << " bestSP=" << bestSP;

			if (len < bestSP) {
				bestSP = len;
//				mylog << "\nbestSP updated to " << bestSP;
				bestEdge = e;
				bestCandidates.clear();
				bestCandidates.insert(bestCandidates.end(), result.begin(),
						result.end());
			}
			result.clear();
		}
	}
// return the best deletions
	result.push_back(bestEdge);
	result.insert(result.end(), bestCandidates.begin(), bestCandidates.end());
	return bestSP;
}

template<class Graph>
bool k_SP(Graph &g, const Vertex<Graph> s, const Vertex<Graph> t, const int k,
		ParentMap &parent) {

	vector<Edge<Graph>> excludeList;

	if (k > 0) { // second  SP and so on
		if (bestNextEdges(g, k, s, t, excludeList) == INF) {
			return false;
		}
		for (auto e : excludeList) {
			put(edge_weight, g, e, INF); // hide the edge from SP algorithms
		}
	}
	size_t len = shortestPath(g, s, t, parent);
//	mylog << "shortestPath=" << len;
//	mylog << " excludeList:" << excludeList.size();
	for (auto e : excludeList) {
		put(edge_weight, g, e, 1); // unhide the excluded edges
	}
	return len < INF && parent[t] != t;
}

#endif /* KSP_HPP_ */
