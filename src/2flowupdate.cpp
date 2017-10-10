#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <random>
#include <boost/graph/random.hpp>
#include <ctime>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "helpers.hpp"

using namespace boost;
using namespace std;

bool isForkNode(const Vertex<myTypes::MyGraph> &v, const myTypes::MyGraph &g,
		FlowID fid) {

	typedef typename graph_traits<myTypes::MyGraph>::out_edge_iterator OutEdgeIter;
	OutEdgeIter out, out_end;
	bool fold = false, fnew = false;
	tie(out, out_end) = out_edges(v, g);

	for (; out != out_end; ++out) {
		Usage_E usage = g[*out].flows[fid].usage;
		fold = fold || (usage == flow_old);
		fnew = fnew || (usage == flow_new);
	}
	return fold && fnew;
}

bool isJoinNode(const Vertex<myTypes::MyGraph> &v, const myTypes::MyGraph &g,
		FlowID fid) {

	typedef typename graph_traits<myTypes::MyGraph>::in_edge_iterator InEdgeIter;
	InEdgeIter in, in_end;
	bool fold = false, fnew = false;
	tie(in, in_end) = in_edges(v, g);

	for (; in != in_end; ++in) {
		Usage_E usage = g[*in].flows[fid].usage;
		fold = fold || (usage == flow_old);
		fnew = fnew || (usage == flow_new);
	}
	return fold && fnew;
}

void computeBlocks(myTypes::MyGraph &g, FlowID fid,
		vector<myTypes::MyGraph*> &blocks) {
	myTypes::VertexList vl;
	topological_sort(g, std::back_inserter(vl));

	cout << "A topological ordering: ";
	for (auto ii = vl.rbegin(); ii != vl.rend(); ++ii)
		cout << get(vertex_name, g)[*ii] << "_" << isForkNode(*ii, g, fid)
				<< "_" << isJoinNode(*ii, g, fid) << " ";

	myTypes::MyGraph *block;
	// scan the nodes in the topological order, create a block on each fork node
	// close the current block on the first join node
	// add any node in between to the current block
	bool added;
	for (auto vi = vl.rbegin(); vi != vl.rend(); ++vi) {
		string name = get(vertex_name, g)[*vi];

		if (isJoinNode(*vi, g, fid)) {
			blocks.push_back(block); // block is done
			cout << "\nadding join node " << name << *vi << endl;
			add_vertex(*vi, *block);
			added = true;
			print_network(*block);
		}

		if (isForkNode(*vi, g, fid)) {
			block = &g.create_subgraph(); // next block
			get_property(*block, graph_name) = fid;
			cout << "\nadding fork node " << name << *vi;
			add_vertex(*vi, *block);
			added = true;

		} else if (!added) {
			cout << "\nadding " << name << *vi;
			add_vertex(*vi, *block);
		}
		added = false;
	}
}

bool edgeExists(const Vertex<myTypes::MyGraph> &u,
		const Vertex<myTypes::MyGraph> &v, const myTypes::MyGraph &g) {
	return !(g.global_to_local(u) == g.null_vertex()
			|| g.global_to_local(v) == g.null_vertex());
}

void computeDependencyGraph(const vector<myTypes::MyGraph*> &blocks,
		myTypes::MyGraph root, myTypes::DAG &dependency) {
	// here only 2 flow-pairs is assumed
	typename graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;

//	cout << "printing b0\n";
//	print_graph(*blocks[0]);
//	cout << "printing b1\n";
//	print_graph(*blocks[1]);
//	cout << "printing b2\n";
//	print_graph(*blocks[2]);

	for (tie(ei, ei_end) = edges(root); ei != ei_end; ++ei) {

		auto enditr = end(root[*ei].flows), beginitr = begin(root[*ei].flows);

		auto oldFlow = find_if(beginitr, enditr,
				[](Flow f) {return f.usage==flow_old;});
		auto newFlow = find_if(beginitr, enditr,
				[](Flow f) {return f.usage==flow_new;});

		int fid_old = oldFlow - beginitr;
		int fid_new = newFlow - beginitr;

		if (oldFlow == enditr || newFlow == enditr || fid_old == fid_new) {
			continue; // no conflict
		}
		// unit demands are assumed
		if (root[*ei].capacity > 1) {
			continue; // sufficient capacity
		}
		cout << "\nhandling target edge ("
				<< get(vertex_name, root)[source(*ei, root)] << ','
				<< get(vertex_name, root)[target(*ei, root)] << ')';

		int b1_idx = 0, b2_idx = 0;
		for (auto b1 = blocks.begin(); b1 != blocks.end(); ++b1, ++b1_idx) {
			int b1_fid = get_property(**b1, graph_name);
			b2_idx = b1_idx + 1;

			for (auto b2 = b1 + 1; b2 != blocks.end(); ++b2, ++b2_idx) {
				int b2_fid = get_property(**b2, graph_name);
				if (b1_fid == b2_fid)
					continue;

				bool flag1, flag2;
				if (!(flag1 = edgeExists(source(*ei, root), target(*ei, root),
						**b1))
						|| !(flag2 = edgeExists(source(*ei, root),
								target(*ei, root), **b2))) {
//					cout << "\nflag1=" << flag1 << " b1_idx=" << b1_idx
//							<< " flag2=" << flag2 << " b2_idx=" << b2_idx
//							<< endl;
					continue;//  b1 and b2 do not share the edge
				}
				// now add the dependency edge
				cout << "\nadding dependency between blocks:" << b1_idx << ','
						<< b2_idx;
				// the edge must point to the block vertex whose old flow is assigned to the link *ei
				if (fid_old == b2_fid && fid_new == b1_fid)
					add_edge(b1_idx, b2_idx, dependency);
				else if (fid_old == b1_fid && fid_new == b2_fid)
					add_edge(b2_idx, b1_idx, dependency);
			}
		}
	}
	//add_edge(0, 1, dependency);
}

struct cycle_detector: public default_dfs_visitor {
	cycle_detector(myTypes::DAG &myDAG, bool &cycle, int &diameter) :
			myDAG(myDAG), has_cycle(cycle), diameter(diameter) {
	}
	void back_edge(const Edge<myTypes::DAG> &e, const myTypes::DAG &g) {
		has_cycle = true;
	}
	void discover_vertex(const Vertex<myTypes::DAG> &v, const myTypes::DAG &g) {
//		cout << "\ndiscover_vertex " << v;
		put(vertex_distance, myDAG, v, 0); // init
	}
	void start_vertex(const Vertex<myTypes::DAG> &v, const myTypes::DAG &g) {
//		cout << "\nstart_vertex " << v;
	}
	void finish_vertex(const Vertex<myTypes::DAG> &v, const myTypes::DAG &g) {
//		cout << "\nfinish_vertex " << v;
	}
	void finish_edge(const Edge<myTypes::DAG> &e, const myTypes::DAG &g) {
//		cout << "\nfinish_edge " << source(e, g);
		int h_src = get(vertex_distance, g)[source(e, g)];
		int h_tar = get(vertex_distance, g)[target(e, g)];
		h_src = max(h_tar + 1, h_src); // update the source vertex height
		put(vertex_distance, myDAG, source(e, g), h_src);
//		cout << "\nheight of v=" << source(e, g) << " is " << h_src;
		diameter = max(diameter, h_src); // the max height so far
	}
	void forward_or_cross_edge(const Edge<myTypes::DAG> &e,
			const myTypes::DAG &g) {
//		cout << "\nforward_or_cross_edge " << target(e, g);
//		int h_ = get(vertex_distance, g)[target(e, g)]; // height value from a previously finished vertex
//		cout << "\nheight of v=" << target(e, g) << " is " << h_;
		//height = h_; // new start value
		// n finish_edge will be called on source(e, g)
	}
	myTypes::DAG &myDAG;
	bool &has_cycle;
	int &diameter;
};

myTypes::Result evaluate(myTypes::DAG &g) {
	bool has_cycle = false;
	int diameter = 0;
	cycle_detector sd(g, has_cycle, diameter);
	depth_first_search(g, visitor(sd));
	return myTypes::Result(has_cycle, diameter);
}

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
//	cout << "\ndistance from " << s << " to " << t << " is " << distMap[t]
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

//	std::cout << "distances and parents:" << std::endl;
//	for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
//		std::cout << "parent(" << get(vertex_name, g)[*vi] << ") = "
//				<< get(vertex_name, g)[parent[*vi]] << endl;
//	}
	do {
		if (count_if(vi, vend,
				[v,target](Vertex<Graph> i) {return target.local_to_global(i)==v;})
				== 0) {  // check for duplicates
			cout << "\nadding vertex " << v;
			add_vertex(v, target);
			Edge<Graph> e = edge(parent[v], v, g).first;
			g[e].flows[flowID] = Flow(usage);
		}
		v = parent[v];
	} while (v != s);
}

/**
 * remove an edge from g such that the shortest path in g'=degenerate(g...) is the second shortest path in g
 * returns false if no such edge exist. That is, no second shortest path exists.
 */
template<class Graph>
bool degenerate(Graph &g, Vertex<Graph> s, Vertex<Graph> t) {
	Edge<Graph> dummy;
	Edge<Graph> &bestEdge = dummy;

	int bestSP = INF;
	vector<Vertex<Graph>> parent(num_vertices(g));
	size_t spLen = shortestPath(g, s, t, parent);

	for (Vertex<Graph> v = t; v != s; v = parent[v]) { // for each path edge from t to s
		cout << "\ncandidate edge " << edgeToStr(parent[v], v, g);

		auto pair = edge(parent[v], v, g);
		auto e = pair.first;
		auto ew = get(edge_weight, g, e);
		put(edge_weight, g, e, INF);
		spLen = shortestPath(g, s, t); // shortest path not using edge e
		put(edge_weight, g, e, ew);

		cout << "\nspLen=" << spLen;
		if (spLen <= bestSP) { // next best SP could have the same length as the current SP
			bestSP = spLen;
			bestEdge = e;
		}
	}

	cout << "\nremoving best edge " << edgeToStr(bestEdge, g);
	put(edge_weight, g, bestEdge, INF);
	print_network(g);

	return bestSP < INF;
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
	bool nextSPExist;

	add_vertex(s, pair1);
	add_vertex(s, pair2);

	cout << "\nold flow for pair1:";
	size_t spLen = shortestPath(g, s, t, f1old);
	if (spLen < INF) {
		addSP(g, pair1, BLUE, flow_old, f1old, s, t);
	} else {
		cout << "\ns,t is disconnected";
		return false;
	}

	cout << "\nnew flow for pair1:";
	nextSPExist = degenerate(g, s, t);
	if (nextSPExist) {
		shortestPath(g, s, t, f1new);
		addSP(g, pair1, BLUE, flow_new, f1new, s, t);

	} else { // no more s,t path is left
		cout << "\n SP2 does not exist";
		return false;
	}

	cout << "\nold flow for pair2:";
	nextSPExist = degenerate(g, s, t);
	if (nextSPExist) {
		shortestPath(g, s, t, f2old);
		addSP(g, pair2, RED, flow_old, f2old, s, t);

	} else { // no more s,t path is left
		cout << "\n SP3 does not exist";
		return false;
	}

	cout << "\nold flow for pair2:";
	nextSPExist = degenerate(g, s, t);
	if (nextSPExist) {
		shortestPath(g, s, t, f2new);
		addSP(g, pair2, RED, flow_old, f2new, s, t);

	} else { // no more s,t path is left
		cout << "\nSP4 does not exist";
		return false;
	}

	return true;
}

int main(int, char*[]) {

	// the network graph
	myTypes::MyGraph g(0);
	// the first flow pair
	myTypes::MyGraph pair1 = g.create_subgraph();
	// the second flow pair
	myTypes::MyGraph pair2 = g.create_subgraph();

	exampleNetwork(g, pair1, pair2);

	vector<myTypes::MyGraph*> blocks;
	computeBlocks(pair1, BLUE, blocks);
	computeBlocks(pair2, RED, blocks);
	cout << "blocks.size()=" << blocks.size() << "num_vertices="
			<< num_vertices(*blocks[0]) << endl;

	myTypes::DAG dep(blocks.size());
	computeDependencyGraph(blocks, g, dep);

//	cout << "\nprinting dependency graph:\n";

//	mt19937 rng;
//	rng.seed(uint32_t(time(0)));
//	generate_random_graph(dep, 10, 10, rng, false, false);

//	print_graph(g);
	myTypes::Result res = evaluate(dep);
	cout << "\nhas_cycle? " << res.first << " diameter=" << res.second << endl;

	myTypes::MyGraph pair = g.create_subgraph();
	generateFlowPairs(g, pair, pair, vertex(example.S, g),
			vertex(example.T, g));
	print_network(pair);

	//print_graph(*blocks[0], get(vertex_index, *blocks[0]));
	return 0;
}
