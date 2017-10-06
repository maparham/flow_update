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

#include "helpers.hpp"

using namespace boost;
using namespace std;


// Make convenient labels for the vertices
struct Example {
	int S;
	int T;
	int U;
	int V;
	int W;
	int numNodes = 5;
} example;

bool isForkNode(const myTypes::Vertex &v, const myTypes::MyGraph &g,
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

bool isJoinNode(const myTypes::Vertex &v, const myTypes::MyGraph &g,
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

bool edgeExists(const myTypes::Vertex &u, const myTypes::Vertex &v,
		const myTypes::MyGraph &g) {
	return !(g.global_to_local(u) == g.null_vertex()
			|| g.global_to_local(v) == g.null_vertex());
}

void computeDependencyGraph(const vector<myTypes::MyGraph*> &blocks,
		myTypes::MyGraph root, myTypes::MyGraph &dependency) {
	// here only 2 flow-pairs is assumed
	typename graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;

	cout << "printing b0\n";
	print_network(*blocks[0]);
	cout << "printing b1\n";
	print_network(*blocks[1]);
	cout << "printing b2\n";
	print_network(*blocks[2]);

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
	cycle_detector(bool & cycle) :
			has_cycle(cycle) {
	}
	void back_edge(const myTypes::Edge &e, const myTypes::MyGraph &g) {
		has_cycle = true;
	}
	bool & has_cycle;
};

bool has_cycle(const myTypes::MyGraph & g) {
	bool has_cycle = false;
	cycle_detector sd(has_cycle);
	depth_first_search(g, visitor(sd));
	return has_cycle;
}

int main(int, char*[]) {

	// declare a graph object
	myTypes::MyGraph g(0);

	put(vertex_name, g, example.S = add_vertex(g), "S");
	put(vertex_name, g, example.T = add_vertex(g), "T");
	put(vertex_name, g, example.W = add_vertex(g), "W");
	put(vertex_name, g, example.U = add_vertex(g), "U");
	put(vertex_name, g, example.V = add_vertex(g), "V");

	myTypes::Edge e;
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
//	// the first flow pair
	myTypes::MyGraph pair1 = g.create_subgraph();
	add_vertex(example.S, pair1);
	add_vertex(example.T, pair1);
	add_vertex(example.W, pair1);
	add_vertex(example.U, pair1);
	add_vertex(example.V, pair1);

	// the second flow pair
	myTypes::MyGraph pair2 = g.create_subgraph();
	add_vertex(example.S, pair2);
	add_vertex(example.T, pair2);
	add_vertex(example.W, pair2);
	add_vertex(example.U, pair2);
	add_vertex(example.V, pair2);

	vector<myTypes::MyGraph*> blocks;
	computeBlocks(pair1, BLUE, blocks);
	computeBlocks(pair2, RED, blocks);
	cout << "blocks.size()=" << blocks.size() << "num_vertices="
			<< num_vertices(*blocks[0]) << endl;

	myTypes::MyGraph dep(blocks.size());
	computeDependencyGraph(blocks, g, dep);

	cout << "\nprinting dependency graph:";
	print_graph(dep);
	cout << "has_cycle? " << has_cycle(dep);

	//print_graph(*blocks[0], get(vertex_index, *blocks[0]));

	return 0;
}
