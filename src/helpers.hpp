#ifndef HELPERS_HPP_
#define HELPERS_HPP_

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
#define DEBUG 1

struct Logger {
	template<typename T>
	Logger& operator <<(const T& x) {
#ifdef DEBUG
		std::cout << x;
#endif
		return *this;
	}
} mylog;

enum FlowID {
	BLUE, RED
};

// define custom property names
enum Usage_E {
	flow_old, flow_new, flow_both, flow_none
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
			GraphProp> DiGraph;

	typedef subgraph<DiGraph> MyGraph;
	typedef std::vector<Vertex<MyGraph>> VertexList;

	typedef adjacency_list<vecS, vecS, bidirectionalS,
			property<vertex_distance_t, size_t>> DAG;

	typedef pair<bool, size_t> Result;

	typedef typename property_map<MyGraph, vertex_index1_t>::type Index1Map;
};

typedef vector<Vertex<myTypes::MyGraph>> ParentMap;

template<class Graph>
string getFlowName(Graph &g, Edge<Graph> e, int fid) {
	if (g[e].flows[fid].usage == flow_old) {
		return to_string(fid) + "_old";

	} else if (g[e].flows[fid].usage == flow_new) {
		return to_string(fid) + "_new";

	} else if (g[e].flows[fid].usage == flow_both) {
		return to_string(fid) + "_both";
	}
	return "";
}

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

// Make convenient labels for the vertices
struct Example {
	int S;
	int T;
	int U;
	int V;
	int W;
	int numNodes = 5;
} example;

void exampleNetwork(myTypes::MyGraph &g, myTypes::MyGraph &pair1,
		myTypes::MyGraph &pair2) {
	put(vertex_name, g, example.S = add_vertex(g), "S");
	put(vertex_name, g, example.T = add_vertex(g), "T");
	put(vertex_name, g, example.W = add_vertex(g), "W");
	put(vertex_name, g, example.U = add_vertex(g), "U");
	put(vertex_name, g, example.V = add_vertex(g), "V");

	Edge<myTypes::MyGraph> e;
	e = add_edge(example.U, example.W, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(example.W, example.V, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(example.S, example.W, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(example.S, example.U, g).first;
	g[e].capacity = 2;
	g[e].flows[RED] = Flow(flow_old);
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(example.W, example.T, g).first;
	g[e].capacity = 2;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(example.U, example.V, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(example.V, example.T, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);
	g[e].flows[BLUE] = Flow(flow_new);

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}

	add_vertex(example.S, pair1);
	add_vertex(example.T, pair1);
	add_vertex(example.W, pair1);
	add_vertex(example.U, pair1);
	add_vertex(example.V, pair1);

	add_vertex(example.S, pair2);
	add_vertex(example.T, pair2);
	add_vertex(example.W, pair2);
	add_vertex(example.U, pair2);
	add_vertex(example.V, pair2);

}

template<class Graph>
string edgeToStr(Vertex<Graph> u, Vertex<Graph> v, Graph g) {
	stringstream sstm;
	auto nameMap = get(vertex_name, g);
	sstm << "(" << nameMap[u] << "," << nameMap[v] << ")";
	return sstm.str();
}

template<class Graph>
string edgeToStr(Edge<Graph> e, Graph g) {
	stringstream sstm;
	return edgeToStr(source(e, g), target(e, g), g);
}

template<class Vertex>
void print_parent_map(ParentMap p, Vertex s, Vertex t) {
	Vertex v = t;
	int n = 50;
	mylog << "\n";
	do {
		mylog << v << ',';
		v = p[v];
	} while (v != s && n--);
	mylog << s;
}

// saving dotfiles for rendering to PNG
template<typename G>
void save_dot_file(std::string const& fname, G& graph) {
	dynamic_properties dp;
//	dp.property("node_id", boost::get(&VertexProperties::id, graph));
//	dp.property("label", boost::get(&VertexProperties::label, graph));
//	dp.property("weight", boost::get(&EdgeProperties::weight, graph));

	std::ofstream ofs(fname);
	write_graphviz_dp(ofs, graph, dp);
}

#endif /* HELPERS_HPP_ */
