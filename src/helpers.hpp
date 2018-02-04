#ifndef HELPERS_HPP_
#define HELPERS_HPP_

#define DEBUG 1

#ifdef DEBUG
#define PRINTF printf
#else
#define PRINTF(format, args...) ((void)0)
#endif

#include <iostream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <boost/graph/graphviz.hpp>

using namespace boost;
using namespace std;

#define INF numeric_limits<int>::max()
#define MINUSINF numeric_limits<int>::min()

struct Logger: std::ostream {
	template<typename T>
	Logger& operator <<(const T& x) {
#ifdef DEBUG
		std::cout << x;
		std::cout.flush();
#endif
		return *this;
	}
} mylog;

enum FlowID {
	BLUE, RED
};
const int fid[] = { BLUE, RED };

// define custom property names
enum Usage_E {
	flow_none, flow_old, flow_new, flow_both
};

string fnames[] = { "old", "new", "both", "none" };

struct Flow {
	Flow() :
			usage(flow_none) {
	}
	Flow(Usage_E usage) :
			usage(usage) {
	}
	Usage_E usage;
	Usage_E operator()() {
		return usage;
	}
};

struct Edge_label {
	Edge_label() :
			capacity(1) {
	}
	int capacity;
	Flow flows[2];
};
struct graph_label {
	int id;
};

template<class Graph>
using Edge= typename graph_traits<Graph>::edge_descriptor;

template<class Graph>
using Vertex = typename graph_traits<Graph>::vertex_descriptor;

template<class G>
using Path = vector<Vertex<G>>;

struct myTypes {
	// create a typedef for the Graph type
	typedef property<edge_index_t, size_t,
			property<edge_weight_t, size_t, Edge_label>> EdgeProp;
	typedef property<vertex_name_t, string,
			property<vertex_discover_time_t, size_t,
					property<vertex_index1_t, size_t>>> VertexProp;
	typedef property<graph_name_t, int> GraphProp;

	typedef adjacency_list<vecS, vecS, bidirectionalS,
			property<vertex_color_t, default_color_type, VertexProp>, EdgeProp,
			GraphProp> BiGraph;

	typedef subgraph<BiGraph> MyGraph;
	typedef std::vector<Vertex<MyGraph>> VertexList;

//	typedef adjacency_list<vecS, vecS, bidirectionalS,
//			property<vertex_distance_t, size_t>> DAG;

	typedef adjacency_list<vecS, vecS, directedS,
			property<vertex_distance_t, size_t>> Directed;

	typedef pair<bool, size_t> Result;

	typedef typename property_map<MyGraph, vertex_index1_t>::type Index1Map;
};

typedef vector<Vertex<myTypes::MyGraph>> ParentMap;

template<class Graph, class PrintLog>
void print_network(Graph& g, PrintLog &plog) {
	plog << "\n";
	typedef typename boost::graph_traits<Graph>::vertex_iterator Viter;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;

	Viter ui, uiend;

	typename property_map<Graph, vertex_name_t>::type name_map = get(
			vertex_name, g);
	typename property_map<Graph, edge_weight_t>::type weight_map = get(
			edge_weight, g);
	Vertex<Graph> v;
	for (tie(ui, uiend) = vertices(g); ui != uiend; ++ui) {
		OutEdgeIter out, out_end;
		v = *ui;
		plog << (name_map[v] == "" ? to_string(v) : name_map[v]) << "\t";

		tie(out, out_end) = out_edges(*ui, g);
		for (; out != out_end; ++out) {
			v = target(*out, g);
			plog << "--(" << getFlowName(g, *out, BLUE) << ","
					<< getFlowName(g, *out, RED) << "/" << g[*out].capacity
					<< "/" << ((weight_map[*out] == INF)?"INF":to_string(weight_map[*out]))
			<< ")--> " << (name_map[v]=="" ?to_string(v): name_map[v] )
			<< "\t";
		}
		InEdgeIter in, in_end;
		plog << "\n" << "\t";
		tie(in, in_end) = in_edges(*ui, g);
//		for (; in != in_end; ++in)
//			mylog << "<--(" << G[*in].flows[0].ID << "/" << G[*in].capacity
//					<< ")-- " << source(*in, G) << "\t";

		plog << "\n";
	}
}

template<class Graph>
void print_network(Graph& g) {
	print_network(g, mylog);
}

template<class Graph>
void print_network1(Graph& g, const string &title) {
	mylog << title << '\n';
	print_network(g, mylog);
}

template<class Graph>
void print_network_forced(Graph& g) {
	print_network(g, std::cout);
}

template<class Graph>
struct FlowEdgeFilter {
	FlowEdgeFilter() {
	}
	FlowEdgeFilter(int fid, Graph &g) :
			fid(fid), g(&g) {
	}
	bool operator()(const Edge<Graph> &e) const {
		return (*g)[e].flows[fid].usage != flow_none;
	}
	int fid;
	Graph *g;
};
template<class Graph>
struct FlowVertexFilter {
	FlowVertexFilter() {
	}
	FlowVertexFilter(int fid, Graph &g) :
			fid(fid), g(&g) {
	}
	bool operator()(const Vertex<Graph> &v) const {
//		// return true only if v has any in/out edge with flow fid
		auto outitr = out_edges(v, *g);
		for (; outitr.first != outitr.second; ++outitr.first) {
			if ((*g)[*outitr.first].flows[fid].usage != flow_none) {
				return true;
			}
		}
		auto initr = in_edges(v, *g);
		for (; initr.first != initr.second; ++initr.first) {
			if ((*g)[*initr.first].flows[fid].usage != flow_none) {
				return true;
			}
		}
		return false;
	}
	int fid;
	Graph *g;
};

typedef filtered_graph<myTypes::MyGraph, FlowEdgeFilter<myTypes::MyGraph>,
		FlowVertexFilter<myTypes::MyGraph>> FlowPair;

template<class Graph>
string edgeToStr(Vertex<Graph> u, Vertex<Graph> v, Graph &g) {
	stringstream sstm;
	auto nameMap = get(vertex_name, g);
//	mylog << "(" << u << "," << v << ")";
//	sstm << "(" << nameMap[u] << "," << nameMap[v] << ")";
	sstm << "(" << u << "," << v << ")";
	return sstm.str();
}

template<class Graph>
string edgeToStr(Edge<Graph> e, Graph &g) {
	stringstream sstm;
	return edgeToStr(source(e, g), target(e, g), g);
}

string edgeToStr(Edge<myTypes::MyGraph> e, myTypes::MyGraph &g) {
	stringstream sstm;
	if (g.is_root()) {
		return edgeToStr(source(e, g), target(e, g), g);

	} else {
		auto e_g = g.local_to_global(e);
		return edgeToStr(e_g, g.root());
	}
}

template<class G = myTypes::MyGraph>
void printPath(Path<G> p) {
	for (int i = 0; i < p.size(); ++i) {
		PRINTF(",%d", p[i]);
	}
	PRINTF("\n");
}

template<class Vertex>
vector<Vertex> getPath(ParentMap &p, Vertex s, Vertex t) {
	vector<Vertex> path;
	Vertex v = t;
	while (p[v] != 1 && v != p[v]) {
		path.push_back(v);
		v = p[v];
	}
	if (v == s) {
		path.push_back(s);
	}
	reverse(path.begin(), path.end());
	return path;
}

template<class Graph>
vector<Vertex<Graph>> k_SP(Graph &g, Vertex<Graph> s, Vertex<Graph> t, int k) {
	ParentMap pm(num_vertices(g), -1);
	if (k_SP(g, s, t, k, pm) == false) {
		return {};
	}
	return getPath(pm, s, t);
}

// saving dotfiles for rendering to PNG
template<typename G>
void save_dot_file(std::string const& fname, G& graph) {
	dynamic_properties dp;
	dp.property("node_id", get(vertex_index, graph));
	dp.property("label", get(vertex_name, graph));
	dp.property("weight", get(edge_weight, graph));

	std::ofstream ofs(fname);
	write_graphviz_dp(ofs, graph, dp);
}

void flowPairs(myTypes::MyGraph& g, vector<FlowPair>& flowpairs) {
	// extract pairs
	FlowEdgeFilter<myTypes::MyGraph> edgeFilter1(BLUE, g), edgeFilter2(RED, g);
	FlowVertexFilter<myTypes::MyGraph> vertexFilter1(BLUE, g), vertexFilter2(
			RED, g);
	FlowPair p_blue(g, edgeFilter1, vertexFilter1);
	FlowPair p_red(g, edgeFilter2, vertexFilter2);
	flowpairs.push_back(p_blue);
	flowpairs.push_back(p_red);
}

template<class G>
bool trivial(G& g) {
	typename graph_traits<G>::edge_iterator e, e_end;
	bool b = true;
	for (tie(e, e_end) = edges(g); e != e_end; ++e) {
		if (g[*e].flows[0].usage == flow_old || g[*e].flows[0].usage == flow_new) {
			b = false;
		}
		if (g[*e].flows[1].usage == flow_old || g[*e].flows[1].usage == flow_new) {
			b = false;
		}
	}
	return b;
}

#endif /* HELPERS_HPP_ */
