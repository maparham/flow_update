#ifndef FLOWPAIRGENERATOR_HPP_
#define FLOWPAIRGENERATOR_HPP_

#include <random>
#include <boost/graph/random.hpp>
#include <ctime>
#include <boost/graph/breadth_first_search.hpp>

#include "helpers.hpp"
#include "ksp.hpp"

using namespace boost;
using namespace std;

// set edge labels in the target graph for the path specified by the parent map
template<class Graph>
void addSP(Graph &g, int flowID, Usage_E usage, const ParentMap &parent,
		Vertex<Graph> s, Vertex<Graph> t) {

	Vertex<Graph> v = t;
	do {
		Edge<Graph> e = edge(parent[v], v, g).first;
		mylog << " ;setting  " << getFlowName(g, e, flowID) << " usage="
				<< usage << " for edge " << edgeToStr(e, g);
		g[e].flows[flowID].usage =
				(g[e].flows[flowID].usage == flow_none) ? usage : flow_both;
		v = parent[v];
	} while (v != s);
}

template<class Graph>
bool generateFlowPairs(Graph &g, Vertex<Graph> s, Vertex<Graph> t) {

	auto k = num_vertices(g);
	// define pairs for flow1 and flow2, as parent maps
	ParentMap f1old(k);
	ParentMap f1new(k);
	ParentMap f2old(k);
	ParentMap f2new(k);

	mylog << "\nold flow for pair1:";
	bool found = k_SP(g, s, t, 0, f1old);
	if (found) {
		print_parent_map(f1old, s, t);
		addSP(g, BLUE, flow_old, f1old, s, t);
		//print_network(pair1);

	} else {
		mylog << "\ns,t is disconnected";
		return false;
	}

	mylog << "\nnew flow for pair1:";
	found = k_SP(g, s, t, 1, f1new);
	if (found) {
		print_parent_map(f1new, s, t);
		addSP(g, BLUE, flow_new, f1new, s, t);
		//print_network(pair1);

	} else { // no more s,t path is left
		mylog << "\n SP2 does not exist.";
		return false;
	}

	mylog << "\nold flow for pair2:";
	found = k_SP(g, s, t, 2, f2old);
	if (found) {
		print_parent_map(f2old, s, t);
		addSP(g, RED, flow_old, f2old, s, t);

	} else { // no more s,t path is left
		mylog << "\n SP3 does not exist, taking flow1_new path instead";
		addSP(g, RED, flow_old, f1new, s, t);
	}
	//print_network(pair2);

	mylog << "\nnew flow for pair2:";
	found = k_SP(g, s, t, 3, f2new);
	if (found) {
		mylog << "nextSPExist\n";
		print_parent_map(f2new, s, t);
		addSP(g, RED, flow_new, f2new, s, t);

	} else { // no more s,t path is left
		mylog << "\nSP4 does not exist, taking flow1_old path instead";
		addSP(g, RED, flow_new, f1old, s, t);
	}
	//print_network(pair2);

	return true;
}

template<class Graph>
bool generateFlowPairs_rand(Graph &g, Vertex<Graph> s, Vertex<Graph> t) {

	auto n = num_vertices(g);
	// define pairs for flow1 and flow2, as parent maps
	ParentMap f1old(n, -1);
	ParentMap f1new(n, -1);
	ParentMap f2old(n, -1);
	ParentMap f2new(n, -1);

	bool found;

	vector<int> perm = { 0, 1, 2, 3 };
//	random_shuffle(perm.begin(), perm.end());

	mylog << "\nold flow for pair1 uses SP:" << perm[0];
	found = k_SP(g, s, t, perm[0], f1old);
	if (found) {
		print_parent_map(f1old, s, t);
		addSP(g, BLUE, flow_old, f1old, s, t);

	} else if (f2new[t] > 0 && f2new[t] != t && f2new[t] < n) { // can we steal from f2new?
		addSP(g, BLUE, flow_old, f2new, s, t);

	} else {
		mylog << "\ns,t is disconnected";
		return false;
	}

	mylog << "\nnew flow for pair1 uses SP:" << perm[1];
	found = k_SP(g, s, t, perm[1], f1new);
	if (found) {
		print_parent_map(f1new, s, t);
		addSP(g, BLUE, flow_new, f1new, s, t);

	} else if (f2old[t] > 0 && f2old[t] != t && f2new[t] < n) {
		addSP(g, BLUE, flow_new, f2old, s, t);

	} else { // no more s,t path is left
		mylog << "\n SP2 does not exist.";
		return false;
	}

	mylog << "\nold flow for pair2 uses SP:" << perm[2];
	found = k_SP(g, s, t, perm[2], f2old);
	if (found) {
		print_parent_map(f2old, s, t);
		addSP(g, RED, flow_old, f2old, s, t);

	} else if (f1new[t] > 0 && f1new[t] != t && f2new[t] < n) {
		addSP(g, RED, flow_old, f1new, s, t);

	} else { // no more s,t path is left
		mylog << "\n SP3 does not exist, taking flow1_new path instead";
		addSP(g, RED, flow_old, f1new, s, t);
	}

	mylog << "\nnew flow for pair2 uses SP:" << perm[3];
	found = k_SP(g, s, t, perm[3], f2new);
	if (found) {
		mylog << "nextSPExist\n";
		print_parent_map(f2new, s, t);
		addSP(g, RED, flow_new, f2new, s, t);

	} else if (f1old[t] > 0 && f1old[t] != t && f2new[t] < n) {
		addSP(g, RED, flow_new, f1old, s, t);

	} else { // no more s,t path is left
		mylog << "\nSP4 does not exist, taking flow1_old path instead";
		addSP(g, RED, flow_new, f1old, s, t);
	}

	return true;
}

void postGenerate(myTypes::MyGraph &g) {
	// set vertex names
	auto indexMap = get(vertex_index, g);
	auto nameMap = get(vertex_name, g);
	typename graph_traits<myTypes::MyGraph>::vertex_iterator vi, vend;
	for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		put(nameMap, *vi, to_string(indexMap[*vi]));
	}

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
		g[*e_it].flows[BLUE].usage = flow_none;
		g[*e_it].flows[RED].usage = flow_none;
	}
}

void setCapacity(myTypes::MyGraph &g) {
	/* initialize random seed: */
	srand(time(NULL));

	graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		if (g[*ei].flows[BLUE].usage == flow_none
				&& g[*ei].flows[RED].usage == flow_none) {
			g[*ei].capacity = -1;

		} else if (g[*ei].flows[BLUE].usage == g[*ei].flows[RED].usage) { // i.e. both flows chose either flow_old or flow_new on this edge
			g[*ei].capacity = 2;

		} else {
			g[*ei].capacity = rand() % 2 + 1;
		}
	}
	// remove non-flow edges
	remove_edge_if(
			[&g](Edge<myTypes::MyGraph> e) {
				return g[e].flows[BLUE].usage == flow_none&& g[e].flows[RED].usage == flow_none;
			}, g);
}

class randomNetwork {
	mt19937 rng;
public:
	randomNetwork() {
		rng.seed(uint32_t(time(0)));
	}
	bool generate(myTypes::MyGraph &g, size_t n_v, size_t n_e) {
		generate_random_graph(g, n_v, n_e, rng, false, false);

		postGenerate(g);

		bool success = generateFlowPairs_rand(g, vertex(0, g), vertex(1, g));
		if (!success) {
			mylog << "\npairs could not be generated";
			return false;
		}
		setCapacity(g);
		return true;
	}
};

#endif /* FLOWPAIRGENERATOR_HPP_ */
