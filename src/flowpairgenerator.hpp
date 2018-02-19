#ifndef FLOWPAIRGENERATOR_HPP_
#define FLOWPAIRGENERATOR_HPP_

#include <random>
#include <boost/graph/random.hpp>
#include <ctime>
#include <boost/graph/breadth_first_search.hpp>

#include <functional>
#include "flowutil.hpp"
#include "ksp.hpp"

using namespace boost;
using namespace std;

// set edge labels in the target graph for the path specified by the parent map
template<class Graph>
void addSP(Graph &g, int flowID, Usage_E usage, const Path<Graph> &path,
		Vertex<Graph> s, Vertex<Graph> t) {

	for (int i = path.size(); i-- > 1;) {
		Edge<Graph> e = edge(path[i - 1], path[i], g).first;
//		mylog << " ;setting  " << getFlowName(g, e, flowID) << " usage="
//				<< usage << " for edge " << edgeToStr(e, g);
		g[e].flows[flowID].usage =
				(g[e].flows[flowID].usage == flow_none) ? usage : flow_both;
	}
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

	vector<int> perm = { 0, 1, 2, 0 };
	random_shuffle(perm.begin(), perm.end());

	mylog << "\nold flow for pair1 uses SP:" << perm[0];
	found = k_SP(g, s, t, perm[0], f1old);
	if (found) {
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
		addSP(g, RED, flow_new, f2new, s, t);

	} else if (f1old[t] > 0 && f1old[t] != t && f2new[t] < n) {
		addSP(g, RED, flow_new, f1old, s, t);

	} else { // no more s,t path is left
		mylog << "\nSP4 does not exist, taking flow1_old path instead";
		addSP(g, RED, flow_new, f1old, s, t);
	}

	return true;
}

void setVertexNames(myTypes::MyGraph &g) {
	// set vertex names
	auto indexMap = get(vertex_index, g);
	auto nameMap = get(vertex_name, g);
	typename graph_traits<myTypes::MyGraph>::vertex_iterator vi, vend;
	for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
		put(nameMap, *vi, to_string(indexMap[*vi]));
	}
}

void postGenerate(myTypes::MyGraph &g) {
	setVertexNames(g);
	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
		g[*e_it].flows[BLUE].usage = flow_none;
		g[*e_it].flows[RED].usage = flow_none;
	}
}
void setCapacityRand(myTypes::MyGraph &g) {
	/* initialize random seed: */
	srand(time(NULL));

	graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		if (g[*ei].flows[BLUE].usage == flow_none
				&& g[*ei].flows[RED].usage == flow_none) {
			g[*ei].capacity = -1;

		} else if (g[*ei].flows[BLUE].usage
				== g[*ei].flows[RED].usage) { // i.e. both flows chose either flow_old or flow_new on this edge
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

		Vertex<myTypes::MyGraph> s = 0, t = 1;
		int bestS = MINUSINF, bestT=INF;
		typename graph_traits<myTypes::MyGraph>::vertex_iterator v,
				vend;
		/*
		 for (tie(v, vend) = vertices(g); v != vend; ++v) {
		 int dout = out_degree(*v, g), din = in_degree(*v, g);
		 //int m = dout - din;
		 if (bestS < dout) {
		 //			if (out_degree(*v, g) > 1) {
		 s = *v;
		 bestS = dout;
		 } else if (bestT < din) {
		 //			} else if (in_degree(*v, g) > 1) {
		 t = *v;
		 bestT = din;
		 }
		 }
		 */

		bool success = generateFlowPairs_rand(g, s, t);
		if (!success) {
			mylog << "\npairs could not be generated";
			return false;
		}
		setCapacityRand(g);
		return true;
	}
};

template<class G = myTypes::MyGraph>
class OverlapFlows {
	G &g;
	const int s;
	const int t;

	Path<G> next_path(Path<G>& basePath, const int mode = 1) {
		printf("basePath.size()=%d\n", basePath.size());
		assert(basePath.size() > 0 && basePath[0] == s && basePath.back() == t);
		vector<Vertex<G>> path;
		int s_j = s, j = 0;
		for (int i = 1; i < basePath.size(); ++i) {
			PRINTF("detour for [%d,%d]\n", s_j, basePath[i]);
			const vector<Vertex<G>>& tmp =
					k_SP(g, s_j, basePath[i], mode);
			if (tmp.size() == 0) {
				PRINTF("no detour, continue\n");
				continue;
			}
			assert(tmp.size() > 2);
//			path.insert(path.end(), basePath.begin() + j, basePath.begin() + i - 1);
			path.insert(path.end(), tmp.begin(), tmp.end() - 1);	// append the detour until one before the last node
			s_j = basePath[i];
			j = i;
		}
		if (path.size() < 1) {
			PRINTF("no next path possible\n");
			return {};
		}

		path.insert(path.end(), basePath.begin() + j, basePath.end()); // append the remaining

		assert(path[0] == s && path.back() == t);
		for (int i = 1; i < path.size(); ++i) {
			assert(path[i - 1] != path[i]);
		}

		if (prefixPath(path, basePath)) {
			PRINTF("the same as basePath, recurse with mode=%d\n", mode + 1);
			return next_path(basePath, mode + 1);
		}

		return path;
	}

public:
	OverlapFlows(G& g, const int s = 0, const int t = 1) :
			g(g), s(s), t(t) {
	}
	bool allocate();
};

template<class G>
bool OverlapFlows<G>::allocate() {
	mylog << "BLUE.old flow:\n";
	Path<G> f1old = k_SP(g, s, t, 0);
	if (f1old.size() > 1) {
		PRINTF("adding path1: ");
		printPath<>(f1old);
		addSP(g, BLUE, flow_old, f1old, s, t);
	} else {
		PRINTF("\ns=%d,t=%d is disconnected", s, t);
		return false;
	}

	mylog << "\nBLUE.new flow:\n";
	Path<G> f1new = next_path(f1old);
	if (f1new.size() > 1) {
		PRINTF("adding path2: ");
		printPath<>(f1new);
		addSP(g, BLUE, flow_new, f1new, s, t);
	} else {
		mylog << "\ns,t is disconnected, no path2";
		return false;
	}

	mylog << "\nRED.old flow:\n";
	Path<G> f2old = next_path(f1new);
	if (f2old.size() > 1) {
		assert(prefixPath(f2old, f1old) == false);
		addSP(g, BLUE, flow_old, f2old, s, t);

	} else {
		mylog << "\ns,t is disconnected(2)";
		return false;
	}

	mylog << "\nRED.new flow:\n";
	Path<G> f2new = next_path(f2old);
	if (f2new.size() > 1) {
		assert(prefixPath(f2new, f1old) == false);
		assert(prefixPath(f2new, f1new) == false);
		addSP(g, BLUE, flow_new, f2new, s, t);
	} else {
		mylog << "\ns,t is disconnected(3)";
		return false;
	}
	// auto set capacities
	setMinimalCapacities(g);
	return true;
}
#endif /* FLOWPAIRGENERATOR_HPP_ */
