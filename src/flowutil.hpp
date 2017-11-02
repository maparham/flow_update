#ifndef FLOWUTIL_HPP_
#define FLOWUTIL_HPP_

#include "helpers.hpp"

template<class Graph>
bool isForkNode(const Vertex<Graph> &v, const Graph &g, int fid) {

	typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIter;
	OutEdgeIter out, out_end;
	bool fold = false, fnew = false;
	tie(out, out_end) = out_edges(v, g);

	for (; out != out_end; ++out) {
		Usage_E usage = g[*out].flows[fid].usage;
		//mylog << "{" << edgeToStr(*out, g) << fnames[usage] << "}";

		fold = fold || (usage == flow_old);
		fnew = fnew || (usage == flow_new);
	}
	return fold && fnew;
}
template<class Graph>
bool isJoinNode(const Vertex<Graph> &v, const Graph &g, int fid) {

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

// check if edge e from the root graph is contained in the subgraph g
bool edgeExists(const Edge<myTypes::MyGraph> &e, const myTypes::MyGraph &g) {
	return g.find_edge(e).second;
//	auto u_global = source(e, g.root());
//	auto v_global = target(e, g.root());
//	auto res1 = g.find_vertex(u_global);
//	auto res2 = g.find_vertex(v_global);
//	return res1.second && res2.second;
}

// check the edge is present in the filtered graph g
bool edgeExists(const Vertex<FlowPair> source, const Vertex<FlowPair> target,
		const FlowPair &g) {
	return edge(source, target, g).second;
}
#endif /* FLOWUTIL_HPP_ */
