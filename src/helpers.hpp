using namespace boost;
using namespace std;

enum FlowID {
	BLUE, RED
};

// define custom property names
enum Usage_E {
	flow_old, flow_new, flow_both, flow_none
};

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

struct myTypes {
	// create a typedef for the Graph type
	typedef property<edge_index_t, int, Edge_label> myEdgeProp_t;
	typedef property<vertex_name_t, string> myVertexProp_t;
	typedef property<graph_name_t, int> myGraphProp_t;

	typedef adjacency_list<vecS, vecS, bidirectionalS,
			property<vertex_color_t, default_color_type, myVertexProp_t>,
			myEdgeProp_t, myGraphProp_t> DiGraph;

	typedef subgraph<DiGraph> MyGraph;
	typedef graph_traits<MyGraph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<MyGraph>::vertex_descriptor Vertex;
	typedef std::vector<Vertex> VertexList;
};

template<class Graph>
void print_network(const Graph& G) {
	typedef typename boost::graph_traits<Graph>::vertex_iterator Viter;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;

	int fid = get_property(G, graph_name);
	Viter ui, uiend;
	tie(ui, uiend) = vertices(G);

	for (; ui != uiend; ++ui) {
		OutEdgeIter out, out_end;
		cout << get(vertex_name, G.root())[G.local_to_global(*ui)] << "\t";

		tie(out, out_end) = out_edges(*ui, G);
		for (; out != out_end; ++out)
			cout << "--("

//					<< count_if(begin(G[*out].flows), end(G[*out].flows),
//							[](Flow f) {return f.usage!=flow_none;}) << ','
					<< ((G[*out].flows[fid].usage == flow_old) ? "old" : "new")
					<< "/" << G[*out].capacity << ")--> "
					<< get(vertex_name, G.root())[G.local_to_global(
							target(*out, G))] << "\t";

		InEdgeIter in, in_end;
		cout << endl << "\t";
		tie(in, in_end) = in_edges(*ui, G);
//		for (; in != in_end; ++in)
//			cout << "<--(" << G[*in].flows[0].ID << "/" << G[*in].capacity
//					<< ")-- " << source(*in, G) << "\t";

		cout << endl;
	}
}
