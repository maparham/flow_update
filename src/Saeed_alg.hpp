#ifndef TWU_HPP_
#define TWU_HPP_

#define DEBUG 0
#define Log2File 0

#if DEBUG
#if Log2File
#define PRINTF(...) fprintf(f,__VA_ARGS__)
#else
#define PRINTF printf
#endif
#else
#define PRINTF(format, args...) ((void)0)
#endif

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
bool computeBlocks(FPair &p, Graph &g, FlowID fid, vector<Block>& blocks) {
	myTypes::VertexList sorted;

	if (!topological_sort(p, std::back_inserter(sorted))) { // if not a DAG
		mylog << "\nnot a DAG, flow=" << fid;
		return false;
	}

//	mylog << "A topological ordering: ";
//	for (auto ii = sorted.rbegin(); ii != sorted.rend(); ++ii) {
//		mylog << *ii << "_" << isForkNode(*ii, p, fid) << "_"
//				<< isJoinNode(*ii, p, fid) << ",";
//	}
//	PRINTF("\n");
	// scan the nodes in the topological order, create a block on each fork node
	// close the current block on the first join node
	// add any node in between to the current block
	bool added = false; // vertex added?
	bool open = false; // block open?
	for (auto vi = sorted.rbegin(); vi != sorted.rend(); ++vi) {
		//string name = get(vertex_name, p)[*vi];

		if (isJoinNode(*vi, p, fid)) { // block is done
//			assert(open);
//			assert(!added);
			Block &block = blocks.back();
			PRINTF("adding join node %d, fid=%d\n", *vi, fid);
			add_vertex(*vi, block);
			added = true;
//			printf("block parent=%dV,%dE\n", num_vertices(*block.m_parent), num_edges(*block.m_parent));
			open = false;
		}

		if (isForkNode(*vi, p, fid)) {
			Block &b = g.create_subgraph();
			blocks.push_back(b); // pushed a copy
			Block &block = blocks.back(); // next block
			get_property(block, graph_name) = fid;
			PRINTF("adding fork node %d, fid=%d\n", *vi, fid);
			add_vertex(*vi, block);
//			added = true;
			open = true;

		} else if (!added) {
			if (!open) {
				PRINTF("no open block...skipping node %d\n", *vi);
				continue;
			}
			PRINTF("adding %d degrees=%d,%d, fid=%d\n", *vi, in_degree(*vi, p), out_degree(*vi, p), fid);
			Block &block = blocks.back();
			add_vertex(*vi, block);
		}
		added = false;
	}
	return blocks.size() > 0;
}

void computeDependencyGraph(const vector<myTypes::MyGraph> &blocks,
		myTypes::MyGraph &root, myTypes::Directed &dependency) {
// here only 2 flow-pairs is assumed
	typename graph_traits<myTypes::MyGraph>::edge_iterator ei, ei_end;

	for (tie(ei, ei_end) = edges(root); ei != ei_end; ++ei) {
		PRINTF("\nhandling target edge %s, cap=%d\n", edgeToStr((*ei), root).c_str(), root[*ei].capacity);

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

		//mylog << "\n# of blocks=" << blocks.size();

		for (int b1 = 0; b1 < blocks.size(); ++b1) {

			int b1_fid = get_property(blocks[b1], graph_name);

			for (int b2 = b1 + 1; b2 < blocks.size(); ++b2) {

				int b2_fid = get_property(blocks[b2], graph_name);
//				PRINTF("b1=%d fid1=%d, b2=%d fid2=%d\n", b1, b1_fid, b2, b2_fid);
				if (b1_fid == b2_fid) {
					PRINTF("same fids=>continue\n");
					continue;
				}

				if (!edgeExists(*ei, blocks[b1])
						|| !edgeExists(*ei, blocks[b2])) {
					PRINTF("blocks %d and %d do not overlap on the edge %s\n",
							b1, b2, edgeToStr(*ei, root).c_str());
					continue; //  b1 and b2 do not share the edge
				}
				PRINTF(" edgeExists in both\n");
				// now add the dependency edge
				// the edge must point to the block vertex whose old flow is assigned to the link *ei
				if (fid_old == b2_fid && fid_new == b1_fid) {
					PRINTF("adding dependency between blocks: %d->%d\n", b1, b2);
					add_edge(b1, b2, dependency);

				} else if (fid_old == b1_fid && fid_new == b2_fid) {
					PRINTF("adding dependency between blocks: %d->%d\n", b2, b1);
					add_edge(b2, b1, dependency);

				} else {
					cout << "something is wrong! **b2_idx=" << b2
							<< " b2_fid=" << b2_fid << "\n";
					assert(0);
				}
			}
		}
	}
//	EI<myTypes::Directed> itr, itr_end;
//	for (tie(itr, itr_end) = edges(dependency); itr != itr_end; ++itr) {
//		const Vertex<myTypes::Directed> &tar = target(*itr, dependency);
//		if (out_degree(tar, dependency) == 0 && newLinkCount[tar] > 1) {
//			put(dependency, *itr, edge_weight, -2); // account for 1 preparation phase
//		}
//	}
}

template<class G, class B>
class BlockGraphVisitor: public default_dfs_visitor {
	const vector<B> &blocks;
	map<int, int> newLinkCount, oldLinkCount; // counters
public:
	struct Result {
		bool cyclic = false;
		int diameter = 0;
		int headBlockIdx = -1;
		int tailBlockIdx = -1;
	};
	BlockGraphVisitor(G &myDAG, const vector<B> &blocks, BlockGraphVisitor<G, B>::Result &result) :
			myDAG(myDAG), blocks(blocks), result(result) {
		for (int b = 0; b < blocks.size(); ++b) {
			int fid = get_property(blocks[b], graph_name);
			int newLinks = countLinks(blocks[b], fid, flow_new);
			int oldLinks = countLinks(blocks[b], fid, flow_old);
			newLinkCount[b] = newLinks;
			oldLinkCount[b] = oldLinks;
		}
	}
	void back_edge(const Edge<G> &e, const G &g) {
		PRINTF("back_edge (%d,%d)", source(e, g), target(e, g));
		result.cyclic = true;
	}
	void discover_vertex(const Vertex<G> &v, const G &g) {
		PRINTF("discover_vertex %d\n", v);
		if (out_degree(v, g) == 0) {
			put(vertex_distance, myDAG, v, 0);
		} else {
			put(vertex_distance, myDAG, v, -INF);
		}
		headMap[v] = v;
		parentMap[v] = v;
	}
	void finish_vertex(const Vertex<G> &v, const G &g) {
		PRINTF("finish_vertex %d\n", v);
		int h = get(vertex_distance, g)[v];
		if (h == -INF) { // v is the first finished vertex in a cycle
			put(vertex_distance, myDAG, v, 0);
			headMap[v] = v;
		}
	}
	void finish_edge(const Edge<G> &e, const G &g) {
		PRINTF("finish_edge (%d,%d)\n", source(e, g), target(e, g));
		int src = source(e, g);
		int tar = target(e, g);
		if(inPath(tar,src)) {
			return;
		}

		int src_dist = get(vertex_distance, g)[src];
		int tar_dist = get(vertex_distance, g)[tar];

		int w = 1;
		w += out_degree(tar, g)==0 && newLinkCount[tar]>1? 1: 0;
		w += in_degree(src, g)==0 && oldLinkCount[src]>1? 1: 0;

		// update the source vertex distance
		if ( tar_dist > -INF && tar_dist + w > src_dist) {
			src_dist = tar_dist + w;
			headMap[src] = headMap[tar];
			parentMap[src] = tar;
			put(vertex_distance, myDAG, source(e, g), src_dist);
			// the max height so far
			if (result.diameter < src_dist) {
				result.diameter = src_dist;
				result.headBlockIdx = headMap[src];
				result.tailBlockIdx = src;
				PRINTF("on longest path, result.diameter=%d\n",result.diameter);
//				assert(inPath(result.tailBlock,result.headBlockIdx));
			}
		}
	}
private:
	Result &result;
	G &myDAG;
	map<int,int> headMap, parentMap;
	bool inPath(const int &s, const int &v) {
		int x = s;
		while (parentMap[x] != x) {
			x = parentMap[x];
			if(x==v) {
				return true;
			}
		}
		return false;
	}
};

template<class Graph, class B>
typename BlockGraphVisitor<Graph, B>::Result evaluate(Graph &g, const vector<B> &blocks) {
	typename BlockGraphVisitor<Graph, B>::Result res;
	BlockGraphVisitor<Graph, B> dv(g, blocks, res);
	depth_first_search(g, visitor(dv));
	return res;
}

template<class G>
std::tuple<bool, int, int> two_flows_update(G& g) {
	vector<FlowPair> fpairs;
	flowPairs(g, fpairs);
	FlowPair p_blue = fpairs[0];
	FlowPair p_red = fpairs[1];
	vector<G> blocks;

	computeBlocks(p_blue, g, BLUE, blocks);
	computeBlocks(p_red, g, RED, blocks);
	if (blocks.size() == 0) { // old and new paths are the same in both pairs
		mylog << "\n bogus instance!\n";
		return {true, -1, -1};
	}
	PRINTF("\nblocks.size()=%d, num_vertices=%d\n", blocks.size(), num_vertices(blocks[0]));
#if DEBUG
	for (int i=0;i<blocks.size();++i) {
		PRINTF("\nblock%d: %luV,%luE, fid=%d\n",i,
				num_vertices(blocks[i]), num_edges(blocks[i]), get_property(blocks[i], graph_name));
		print_network_forced(blocks[i]);
	}
#endif
	myTypes::Directed dep(blocks.size());
	computeDependencyGraph(blocks, g, dep);

#if DEBUG
	printf("\nprinting dependency graph\n");
	print_graph(dep);
#endif
	auto res = evaluate(dep, blocks);
	PRINTF("res.diameter=%d, res.tailBlockIdx=%d, res.headBlockIdx=%d\n", res.diameter, res.tailBlockIdx,
			res.headBlockIdx);

	int rounds = res.diameter + 1; // since we have at least one block

	// when rounds<3, it could be the case that some block needs 3 rounds anyway
	for (int i = 0; i < blocks.size(); ++i) {
		auto &b = blocks[i];
		const auto [itr_begin, itr_end] = edges(b);
		for (int fid : { BLUE, RED }) {
			int c_newlinks = count_if(itr_begin, itr_end, [&](Edge<G> e) {
				return b[e].flows[fid].usage==flow_new;
			});
			int c_oldlinks = count_if(itr_begin, itr_end, [&](Edge<G> e) {
				return b[e].flows[fid].usage==flow_old;
			});
			rounds = max(rounds, (c_newlinks > 1) + (c_oldlinks > 1) + 1); // ensure the minimum required rounds
		}
	}
	return {res.cyclic, rounds, blocks.size()};
}

#undef DEBUG
#endif
