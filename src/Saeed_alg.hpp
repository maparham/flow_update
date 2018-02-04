#ifndef TWU_HPP_
#define TWU_HPP_

#include <utility>                   // for std::pair
#include <tuple>
#include <algorithm>                 // for std::for_each
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "flowutil.hpp"
#include "flowpairgenerator.hpp"
#include "helpers.hpp"

using namespace boost;
using namespace std;

template<class FPair, class Graph, class Block>
bool computeBlocks(FPair &p, Graph &g, FlowID fid, vector<Block*> &blocks) {
	myTypes::VertexList sorted;

	if (!topological_sort(p, std::back_inserter(sorted))) { // if not a DAG
		mylog << "\nnot a DAG, flow=" << fid;
		return false;
	}

	mylog << "A topological ordering: ";
	for (auto ii = sorted.rbegin(); ii != sorted.rend(); ++ii)
		mylog << *ii << "_" << isForkNode(*ii, p, fid) << "_"
				<< isJoinNode(*ii, p, fid) << ",";

	Block *block = NULL;
	// scan the nodes in the topological order, create a block on each fork node
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
			block = &(g.create_subgraph()); // next block
//			add_vertex(0, *block);
//			add_vertex(3, *block);
//			mylog << "\n*block:\n";
//			print_graph(*block);
//			edgeExists(*edges(g).first, *block);

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
	return blocks.size() > 0;
}

void computeDependencyGraph(const vector<myTypes::MyGraph*> &blocks,
		myTypes::MyGraph &root, myTypes::Directed &dependency) {
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

		mylog << "\nhandling target edge " << edgeToStr((*ei), root);

		mylog << "\n# of blocks=" << blocks.size();

		for (int b1 = 0; b1 < blocks.size(); ++b1) {

			int b1_fid = get_property(*blocks[b1], graph_name);

			for (int b2 = b1 + 1; b2 < blocks.size(); ++b2) {

				int b2_fid = get_property(*blocks[b2], graph_name);
				if (b1_fid == b2_fid) {
					mylog << "same fids=>continue";
					continue;
				}

				if (!edgeExists(*ei, *blocks[b1])
						|| !edgeExists(*ei, *blocks[b2])) {
					mylog << "blocks do not overlap on the edge "
							<< edgeToStr(*ei, root);
					continue; //  b1 and b2 do not share the edge
				}
				mylog << " edgeExists in both";
				// now add the dependency edge
				// the edge must point to the block vertex whose old flow is assigned to the link *ei
				if (fid_old == b2_fid && fid_new == b1_fid) {
					mylog << "\nadding dependency between blocks:" << b1 << "->"
							<< b2;
					add_edge(b1, b2, dependency);

				} else if (fid_old == b1_fid && fid_new == b2_fid) {
					mylog << "\nadding dependency between blocks1:" << b2
							<< "->" << b1;
					add_edge(b2, b1, dependency);

				} else {
					mylog << "\nsomething is wrong! **b2_idx=" << b2
							<< " b2_fid=" << b2_fid;
				}
			}
		}
	}
//	mylog << "dummy dependency";
//	add_edge(1, 0, dependency);
}

template<class Graph>
struct cycle_detector: public default_dfs_visitor {
	cycle_detector(Graph &myDAG, bool &cycle, int &diameter) :
			myDAG(myDAG), has_cycle(cycle), diameter(diameter) {
	}
	void back_edge(const Edge<Graph> &e, const Graph &g) {
		has_cycle = true;
	}
	void discover_vertex(const Vertex<Graph> &v, const Graph &g) {
//		mylog << "\ndiscover_vertex " << v;
		put(vertex_distance, myDAG, v, 0); // init
	}
	void start_vertex(const Vertex<Graph> &v, const Graph &g) {
//		mylog << "\nstart_vertex " << v;
	}
	void finish_vertex(const Vertex<Graph> &v, const Graph &g) {
//		mylog << "\nfinish_vertex " << v;
	}
	void finish_edge(const Edge<Graph> &e, const Graph &g) {
//		mylog << "\nfinish_edge " << source(e, g);
		int h_src = get(vertex_distance, g)[source(e, g)];
		int h_tar = get(vertex_distance, g)[target(e, g)];
		h_src = max(h_tar + 1, h_src); // update the source vertex height
		put(vertex_distance, myDAG, source(e, g), h_src);
//		mylog << "\nheight of v=" << source(e, g) << " is " << h_src;
		diameter = max(diameter, h_src); // the max height so far
	}
	void forward_or_cross_edge(const Edge<Graph> &e, const Graph &g) {
//		mylog << "\nforward_or_cross_edge " << target(e, g);
//		int h_ = get(vertex_distance, g)[target(e, g)]; // height value from a previously finished vertex
//		mylog << "\nheight of v=" << target(e, g) << " is " << h_;
		//height = h_; // new start value
		// n finish_edge will be called on source(e, g)
	}
	Graph &myDAG;
	bool &has_cycle;
	int &diameter;
};

template<class Graph>
myTypes::Result evaluate(Graph &g) {
	bool has_cycle = false;
	int diameter = 0;
	cycle_detector<Graph> sd(g, has_cycle, diameter);
	depth_first_search(g, visitor(sd));
	return myTypes::Result(has_cycle, diameter);
}

std::tuple<bool, int, int> two_flows_update(myTypes::MyGraph& g) {
	vector<FlowPair> fpairs;
	flowPairs(g, fpairs);
	FlowPair p_blue = fpairs[0];
	FlowPair p_red = fpairs[1];

	//		print_network1(p_blue, "BLUE pair:");
	//		print_network1(p_red, "RED pair:");

	vector<myTypes::MyGraph*> blocks;

	computeBlocks(p_blue, g, BLUE, blocks);
	computeBlocks(p_red, g, RED, blocks);
	if (blocks.size() == 0) { // old and new paths are the same in both pairs
		mylog << "\n bogus instance!\n";
		return {true, -1,-1};
	}
	mylog << "\nblocks.size()=" << blocks.size() << "num_vertices="
			<< num_vertices(*blocks[0]) << "\n";

#ifdef DEBUG1
	for (auto *b : blocks) {
		mylog << "\nprinting block for flow: "
				<< get_property(*b, graph_name) << "\n";
		print_network(*b);
	}
#endif
	myTypes::Directed dep(blocks.size());

	computeDependencyGraph(blocks, g, dep);

	mylog << "\nprinting dependency graph:\n";
	print_graph(dep);

	myTypes::Result res = evaluate(dep);
	mylog << "\nhascycle=" << res.first << '\n';
	return {res.first, res.second, blocks.size()};
}

#endif
