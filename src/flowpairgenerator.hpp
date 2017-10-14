#ifndef FLOWPAIRGENERATOR_HPP_
#define FLOWPAIRGENERATOR_HPP_

#include <random>
#include <boost/graph/random.hpp>
#include <ctime>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "helpers.hpp"

using namespace boost;
using namespace std;

template<class Graph>
size_t shortestPath(const Graph &g, Vertex<Graph> s, Vertex<Graph> t,
		vector<Vertex<Graph>> &parentMap) {
	vector<int> distMap(num_vertices(g));
	Vertex<Graph> s_ = vertex(s, g);
	dijkstra_shortest_paths(g, s_,
			predecessor_map(
					make_iterator_property_map(parentMap.begin(),
							get(vertex_index, g))).distance_map(
					make_iterator_property_map(distMap.begin(),
							get(vertex_index, g))));
//	mylog << "\ndistance from " << s << " to " << t << " is " << distMap[t]
//			<< " in the graph:\n";
//	print_graph(g);
	return distMap[t];
}

template<class Graph>
size_t shortestPath(Graph &g, Vertex<Graph> s, Vertex<Graph> t) {
	vector<Vertex<Graph>> p(num_vertices(g));
	return shortestPath(g, s, t, p);
}

// set edge labels in the target graph for the path specified by the parent map
template<class Graph>
void addSP(Graph &g, Graph &target, int flowID, Usage_E usage, ParentMap parent,
		Vertex<Graph> s, Vertex<Graph> t) {

	Vertex<Graph> v = t;
	typename graph_traits<Graph>::vertex_iterator vi, vend;
	tie(vi, vend) = vertices(target);

	do {
		if (count_if(vi, vend,
				[v,target](Vertex<Graph> i) {return target.local_to_global(i)==v;})
				== 0) {  // check for duplicates
			mylog << "\nadding vertex " << v;
			add_vertex(v, target);
		}
		Edge<Graph> e = edge(parent[v], v, g).first;
		mylog << " ;setting  " << getFlowName(g, e, flowID) << " usage="
				<< usage << " for edge " << edgeToStr(e, g);
		g[e].flows[flowID].usage =
				(g[e].flows[flowID].usage == flow_none) ? usage : flow_both;
		v = parent[v];
	} while (v != s);
}

/**
 * a k shortest path solution
 * in order to get the kth SP, remove an edge from previous SPs
 * try all remove combinations
 */
template<class Graph>
size_t bestNextEdges(Graph &g, int d, Vertex<Graph> s, Vertex<Graph> t,
		vector<Edge<Graph>> &result) {
	Edge<Graph> bestEdge; //  best edge in this call
	vector<Edge<Graph>> bestCandidates; // best edge selection from recursive calls

	int bestSP = INF;
	vector<Vertex<Graph>> parent(num_vertices(g));
	size_t len = shortestPath(g, s, t, parent);
	if (len == INF) {
		return len;
	}
	for (Vertex<Graph> v = t; v != s; v = parent[v]) { // for each path edge from t to s

		auto pair = edge(parent[v], v, g);
		auto e = pair.first;
		auto ew = get(edge_weight, g, e);
		mylog << "\nhide edge " << edgeToStr(parent[v], v, g);
		put(edge_weight, g, e, INF); // hide the edge

		if (d == 0) {
			mylog << " \nbuttom level reached";
			len = shortestPath(g, s, t); // shortest path not using e
		} else {
			len = bestNextEdges(g, d - 1, s, t, result); // SP not using e and the next best edges selected recursively (i.e candidates)
			mylog << " \n  bestNextEdges() returned with len=" << len;
		}

		mylog << "\nuhide edge " << edgeToStr(parent[v], v, g);
		put(edge_weight, g, e, ew); // unhide the edge
		mylog << "\nspLen=" << len;

		if (len <= bestSP) {
			bestSP = len;
			bestEdge = e;
			bestCandidates.clear();
			bestCandidates.insert(bestCandidates.end(), result.begin(),
					result.end());
		}
		result.clear();
	}

	//mylog << "\nremoving best edge " << edgeToStr(bestEdge, g);
	//put(edge_weight, g, bestEdge, INF);
	//remove_edge(bestEdge, g);
	result.push_back(bestEdge);
	result.insert(result.end(), bestCandidates.begin(), bestCandidates.end());
	return bestSP;
}

template<class Graph>
bool k_SP(Graph &g, Vertex<Graph> s, Vertex<Graph> t, int k,
		ParentMap &parent) {

	vector<Edge<Graph>> excludeList;

	if (k > 0) { // second  SP and so on
		bestNextEdges(g, k - 1, s, t, excludeList);
		for (auto e : excludeList) {
			put(edge_weight, g, e, INF); // hide the edge from SP algorithms
		}
	}
	size_t len = shortestPath(g, s, t, parent);
	mylog << "shortestPath=" << len << "INF=" << INF;
	for (auto e : excludeList) {
		put(edge_weight, g, e, 1); // unhide the excluded edges
	}
	return len < INF;
}

template<class Graph>
bool generateFlowPairs(Graph &g, Graph &pair1, Graph &pair2, Vertex<Graph> s,
		Vertex<Graph> t) {

	auto k = num_vertices(g);
	// define pairs for flow1 and flow2, as parent maps
	ParentMap f1old(k);
	ParentMap f1new(k);
	ParentMap f2old(k);
	ParentMap f2new(k);

	add_vertex(s, pair1);
	add_vertex(s, pair2);

	mylog << "\nold flow for pair1:";
	bool found = k_SP(g, s, t, 0, f1old);
	if (found) {
		print_parent_map(f1old, s, t);
		addSP(g, pair1, BLUE, flow_old, f1old, s, t);
		//print_network(pair1);

	} else {
		mylog << "\ns,t is disconnected";
		return false;
	}

	mylog << "\nnew flow for pair1:";
	found = k_SP(g, s, t, 1, f1new);
	if (found) {
		print_parent_map(f1new, s, t);
		addSP(g, pair1, BLUE, flow_new, f1new, s, t);
		//print_network(pair1);

	} else { // no more s,t path is left
		mylog << "\n SP2 does not exist.";
		return false;
	}

	mylog << "\nold flow for pair2:";
	found = k_SP(g, s, t, 2, f2old);
	if (found) {
		print_parent_map(f2old, s, t);
		addSP(g, pair2, RED, flow_old, f2old, s, t);

	} else { // no more s,t path is left
		mylog << "\n SP3 does not exist, taking flow1_new path instead";
		addSP(g, pair2, RED, flow_old, f1new, s, t);
	}
	//print_network(pair2);

	mylog << "\nnew flow for pair2:";
	found = k_SP(g, s, t, 3, f2new);
	if (found) {
		mylog << "nextSPExist\n";
		print_parent_map(f2new, s, t);
		addSP(g, pair2, RED, flow_new, f2new, s, t);

	} else { // no more s,t path is left
		mylog << "\nSP4 does not exist, taking flow1_old path instead";
		addSP(g, pair2, RED, flow_new, f1old, s, t);
	}
	//print_network(pair2);

	return true;
}

template<class Graph>
bool generateFlowPairs1(Graph &g, Graph &pair1, Graph &pair2, Vertex<Graph> s,
		Vertex<Graph> t) {

	auto k = num_vertices(g);
	// define pairs for flow1 and flow2, as parent maps
	ParentMap f1old(k);
	ParentMap f1new(k);
	ParentMap f2old(k);
	ParentMap f2new(k);

	add_vertex(s, pair1);
	add_vertex(s, pair2);

	mylog << "\nold flow for pair1:";
	bool found = k_SP(g, s, t, 0, f1old);
	if (found) {
		print_parent_map(f1old, s, t);
		addSP(g, pair1, BLUE, flow_old, f1old, s, t);
		//print_network(pair1);

	} else {
		mylog << "\ns,t is disconnected";
		return false;
	}

	mylog << "\nold flow for pair2:";
	found = k_SP(g, s, t, 1, f2old);
	if (found) {
		print_parent_map(f2old, s, t);
		addSP(g, pair2, RED, flow_old, f2old, s, t);

	} else { // no more s,t path is left
		mylog << "\n SP2 does not exist.";
		return false;
	}

	mylog << "\nnew flow for pair1:";
	found = k_SP(g, s, t, 2, f1new);
	if (found) {
		print_parent_map(f1new, s, t);
		addSP(g, pair1, BLUE, flow_new, f1new, s, t);
		//print_network(pair1);

	} else { // no more s,t path is left
		mylog << "\n SP3 does not exist, taking flow2_old path instead";
		addSP(g, pair1, BLUE, flow_new, f2old, s, t);
	}

	//print_network(pair2);

	mylog << "\nnew flow for pair2:";
	found = k_SP(g, s, t, 3, f2new);
	if (found) {
		mylog << "nextSPExist\n";
		print_parent_map(f2new, s, t);
		addSP(g, pair2, RED, flow_new, f2new, s, t);

	} else { // no more s,t path is left
		mylog << "\nSP4 does not exist, taking flow1_old path instead";
		addSP(g, pair2, RED, flow_new, f1old, s, t);
	}
	//print_network(pair2);

	return true;
}

struct cycle_remover: public default_dfs_visitor {
	cycle_remover(myTypes::MyGraph &g_) :
			g_(g_) {
	}
	void back_edge(const Edge<myTypes::MyGraph> &e, const myTypes::MyGraph &g) {
		feedback_vertex_map[source(e, g)] = get(vertex_discover_time, g,
				target(e, g));
	}
	void discover_vertex(const Vertex<myTypes::MyGraph> &v,
			const myTypes::MyGraph &g) {
//		mylog << "\ndiscover_vertex " << v;
		put(vertex_discover_time, g_, v, time++); // init
	}
	void start_vertex(const Vertex<myTypes::MyGraph> &v,
			const myTypes::MyGraph &g) {
//		mylog << "\nstart_vertex " << v;
	}
	void finish_vertex(const Vertex<myTypes::MyGraph> &v,
			const myTypes::MyGraph &g) {
//		mylog << "\nfinish_vertex " << v;
	}
	void finish_edge(const Edge<myTypes::MyGraph> &e,
			const myTypes::MyGraph &g) {
	}
	myTypes::MyGraph &g_;
	myTypes::Index1Map feedback_vertex_map;
	size_t time = 0;
};

template<class Graph>
bool generateFlowPairs2(Graph &g, Graph &pair1, Graph &pair2, Vertex<Graph> s,
		Vertex<Graph> t) {

	auto k = num_vertices(g);
	// define pairs for flow1 and flow2, as parent maps
	ParentMap f1old(k);
	ParentMap f1new(k);
	ParentMap f2old(k);
	ParentMap f2new(k);

	add_vertex(s, pair1);
	add_vertex(s, pair2);
	bool found;
	/* initialize random seed: */
	srand(time(NULL));

	for (int i : { 0, 1, 2, 3 }) {

		int rnd = rand() % 4;
		switch (rnd) {
		case 0: {
			mylog << "\nold flow for pair1:";
			found = k_SP(g, s, t, 0, f1old);
			if (found) {
				print_parent_map(f1old, s, t);
				addSP(g, pair1, BLUE, flow_old, f1old, s, t);
				//print_network(pair1);

			} else {
				mylog << "\ns,t is disconnected";
				return false;
			}
			break;
		}
		case 1: {
			mylog << "\nnew flow for pair1:";
			found = k_SP(g, s, t, 1, f1new);
			if (found) {
				print_parent_map(f1new, s, t);
				addSP(g, pair1, BLUE, flow_new, f1new, s, t);
				//print_network(pair1);

			} else { // no more s,t path is left
				mylog << "\n SP2 does not exist.";
				return false;
			}
			break;
		}
		case 2: {
			mylog << "\nold flow for pair2:";
			found = k_SP(g, s, t, 2, f2old);
			if (found) {
				print_parent_map(f2old, s, t);
				addSP(g, pair2, RED, flow_old, f2old, s, t);

			} else { // no more s,t path is left
				mylog << "\n SP3 does not exist, taking flow1_new path instead";
				addSP(g, pair2, RED, flow_old, f1new, s, t);
			}
			//print_network(pair2);
			break;
		}
		case 3: {
			mylog << "\nnew flow for pair2:";
			found = k_SP(g, s, t, 3, f2new);
			if (found) {
				mylog << "nextSPExist\n";
				print_parent_map(f2new, s, t);
				addSP(g, pair2, RED, flow_new, f2new, s, t);

			} else { // no more s,t path is left
				mylog << "\nSP4 does not exist, taking flow1_old path instead";
				addSP(g, pair2, RED, flow_new, f1old, s, t);
			}
			//print_network(pair2);
		}
		}
	}
	return true;
}

void postGenerate(myTypes::MyGraph &g) {
	graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		if (g[*ei].flows[BLUE].usage == flow_none
				&& g[*ei].flows[RED].usage == flow_none) {
			g[*ei].capacity = -1;
			//remove_edge(source(*ei, g), target(*ei, g), g);

		} else if (g[*ei].flows[BLUE].usage == g[*ei].flows[RED].usage) { // i.e. both flows chose either flow_old or flow_new on this edge
			g[*ei].capacity = 2;
		}
	}
	remove_edge_if(
			[&g](Edge<myTypes::MyGraph> e) {
				return g[e].flows[BLUE].usage == flow_none&& g[e].flows[RED].usage == flow_none;
			}, g);

	auto indexMap = get(vertex_index, g);
	auto nameMap = get(vertex_name, g);
	typename graph_traits<myTypes::MyGraph>::vertex_iterator vi, vend;
	for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		put(nameMap, *vi, to_string(indexMap[*vi]));
	}
}

#endif /* FLOWPAIRGENERATOR_HPP_ */
