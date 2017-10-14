#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "helpers.hpp"
#include "flowpairgenerator.hpp"

using namespace boost;
using namespace std;

template<class Graph>
bool isForkNode(const Vertex<Graph> &v, const Graph &g, FlowID fid) {

	typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIter;
	OutEdgeIter out, out_end;
	bool fold = false, fnew = false;
	tie(out, out_end) = out_edges(v, g);

	for (; out != out_end; ++out) {
		Usage_E usage = g[*out].flows[fid].usage;
		mylog << "{" << edgeToStr(*out, g) << fnames[usage] << "}";

		fold = fold || (usage == flow_old);
		fnew = fnew || (usage == flow_new);
	}
	return fold && fnew;
}
template<class Graph>
bool isJoinNode(const Vertex<Graph> &v, const Graph &g, FlowID fid) {

	typedef typename graph_traits<Graph>::in_edge_iterator InEdgeIter;
	InEdgeIter in, in_end;
	bool fold = false, fnew = false;
	tie(in, in_end) = in_edges(v, g);

	for (; in != in_end; ++in) {
		Usage_E usage = g[*in].flows[fid].usage;
		mylog << "{" << edgeToStr(*in, g) << fnames[usage] << "}";
		fold = fold || (usage == flow_old);
		fnew = fnew || (usage == flow_new);
	}
	return fold && fnew;
}

template<class FPair, class Graph, class Block>
void computeBlocks(FPair &p, Graph &g, FlowID fid, vector<Block*> &blocks) {
	myTypes::VertexList sorted;
	topological_sort(p, std::back_inserter(sorted));

	mylog << "A topomylogical ordering: ";
	for (auto ii = sorted.rbegin(); ii != sorted.rend(); ++ii)
		mylog << *ii << "_" << isForkNode(*ii, p, fid) << "_"
				<< isJoinNode(*ii, p, fid) << ",";

	Block *block = NULL;
	// scan the nodes in the topomylogical order, create a block on each fork node
	// close the current block on the first join node
	// add any node in between to the current block
	bool added;
	for (auto vi = sorted.rbegin(); vi != sorted.rend(); ++vi) {
		//string name = get(vertex_name, p)[*vi];

		if (isJoinNode(*vi, p, fid)) {
			blocks.push_back(block); // block is done
			mylog << "\nadding join node " << *vi << "\n";
			add_vertex(*vi, *block);
			added = true;
			//print_network(*block);
			block = NULL;
		}

		if (isForkNode(*vi, p, fid)) {
			block = &g.create_subgraph(); // next block
			get_property(*block, graph_name) = fid;
			mylog << "\nadding fork node " << *vi;
			add_vertex(*vi, *block);
			added = true;

		} else if (!added) {
			mylog << "\nadding " << *vi << " degrees=" << in_degree(*vi, p)
					<< "," << out_degree(*vi, p) << " ";
			if (block == NULL) {
				mylog << "\nblock graph was not created...skipping this node\n";
				continue;
			}
			add_vertex(*vi, *block);
		}
		added = false;
	}
}

bool edgeExists(const Edge<myTypes::MyGraph> &e, const myTypes::MyGraph &g) {
	auto ee = g.global_to_local(e);
	return edge(source(ee, g), target(ee, g), g).second;
}

void computeDependencyGraph(const vector<myTypes::MyGraph*> &blocks,
		myTypes::MyGraph root, myTypes::DAG &dependency) {
	// here only 2 flow-pairs is assumed
	typename graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;

//	mylog << "printing b0\n";
//	print_graph(*blocks[0]);
//	mylog << "printing b1\n";
//	print_graph(*blocks[1]);
//	mylog << "printing b2\n";
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
			continue; // sufficient capacity for 2 flows
		}

		mylog << "\nhandling target edge " << edgeToStr(*ei, root);

		int b1_idx = 0, b2_idx = 0;
		for (auto b1 = blocks.begin(); b1 != blocks.end(); ++b1, ++b1_idx) {
			int b1_fid = get_property(**b1, graph_name);
			b2_idx = b1_idx + 1;

			for (auto b2 = b1 + 1; b2 != blocks.end(); ++b2, ++b2_idx) {
				int b2_fid = get_property(**b2, graph_name);
				if (b1_fid == b2_fid)
					continue;

				bool flag1, flag2;
				if (!(flag1 = edgeExists(*ei, **b1))
						|| !(flag2 = edgeExists(*ei, **b2))) {
//					mylog << "\nflag1=" << flag1 << " b1_idx=" << b1_idx
//							<< " flag2=" << flag2 << " b2_idx=" << b2_idx
//							<< "\n";
					continue;//  b1 and b2 do not share the edge
				}
				// now add the dependency edge
				// the edge must point to the block vertex whose old flow is assigned to the link *ei
				if (fid_old == b2_fid && fid_new == b1_fid) {
					mylog << "\nadding dependency between blocks:" << b1_idx
							<< "->" << b2_idx;
					add_edge(b1_idx, b2_idx, dependency);

				} else if (fid_old == b1_fid && fid_new == b2_fid) {
					mylog << "\nadding   dependency between blocks:" << b2_idx
							<< "->" << b1_idx;
					add_edge(b2_idx, b1_idx, dependency);
				}
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
//		mylog << "\ndiscover_vertex " << v;
		put(vertex_distance, myDAG, v, 0); // init
	}
	void start_vertex(const Vertex<myTypes::DAG> &v, const myTypes::DAG &g) {
//		mylog << "\nstart_vertex " << v;
	}
	void finish_vertex(const Vertex<myTypes::DAG> &v, const myTypes::DAG &g) {
//		mylog << "\nfinish_vertex " << v;
	}
	void finish_edge(const Edge<myTypes::DAG> &e, const myTypes::DAG &g) {
//		mylog << "\nfinish_edge " << source(e, g);
		int h_src = get(vertex_distance, g)[source(e, g)];
		int h_tar = get(vertex_distance, g)[target(e, g)];
		h_src = max(h_tar + 1, h_src); // update the source vertex height
		put(vertex_distance, myDAG, source(e, g), h_src);
//		mylog << "\nheight of v=" << source(e, g) << " is " << h_src;
		diameter = max(diameter, h_src); // the max height so far
	}
	void forward_or_cross_edge(const Edge<myTypes::DAG> &e,
			const myTypes::DAG &g) {
//		mylog << "\nforward_or_cross_edge " << target(e, g);
//		int h_ = get(vertex_distance, g)[target(e, g)]; // height value from a previously finished vertex
//		mylog << "\nheight of v=" << target(e, g) << " is " << h_;
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

int main(int, char*[]) {

	// the network graph
	myTypes::MyGraph g(0);
	vector<myTypes::MyGraph*> blocks;

	do {
		g = myTypes::MyGraph(0);

		// the first flow pair
		myTypes::MyGraph pair1 = g.create_subgraph();
		// the second flow pair
		myTypes::MyGraph pair2 = g.create_subgraph();

		//exampleNetwork(g, pair1, pair2);

		mt19937 rng;
		rng.seed(uint32_t(time(0)));
		generate_random_graph(g, 5, 15, rng, false, false);
		graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
		for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
			put(edge_weight, g, *e_it, 1);
			g[*e_it].flows[BLUE].usage = flow_none;
			g[*e_it].flows[RED].usage = flow_none;
		}

		bool success = generateFlowPairs1(g, pair1, pair2, vertex(0, g),
				vertex(1, g));

		if (!success) {
			mylog << "\npairs could not be generated";
			//return 0;
			continue;
		}

		postGenerate(g);
		mylog << "\ngenerated flow pairs:\n";
		print_network(g);

		FlowEdgeFilter<myTypes::MyGraph> edgeFilter1(BLUE, g), edgeFilter2(RED,
				g);
		FlowVertexFilter<myTypes::MyGraph> vertexFilter1(BLUE, g),
				vertexFilter2(RED, g);

		FlowPair p_blue(g, edgeFilter1, vertexFilter1);
		FlowPair p_red(g, edgeFilter2, vertexFilter2);
		//print_network(p_blue);

		blocks.clear();
		computeBlocks(p_blue, g, BLUE, blocks);
		computeBlocks(p_red, g, RED, blocks);
		mylog << "blocks.size()=" << blocks.size() << "num_vertices="
				<< num_vertices(*blocks[0]) << "\n";

	} while (blocks.size() < 3);

	print_network_forced(g);

	for (auto *b : blocks) {
		cout << "\nprinting block for flow: " << get_property(*b, graph_name) << endl;
		print_graph(*b, get(vertex_name, *b));
	}

	myTypes::DAG dep(blocks.size());
	computeDependencyGraph(blocks, g, dep);
	cout << "\nprinting dependency graph:\n";
	print_graph(dep);
	myTypes::Result res = evaluate(dep);
	cout << "\nhas_cycle? " << res.first << " diameter=" << res.second << "\n";

	return 0;
}
