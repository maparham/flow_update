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
//		mylog << "{" << edgeToStr(*in, g) << fnames[usage] << "}";
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

void setMinimalCapacity(Edge_label& el) {
	el.capacity = 1;
	if (el.flows[BLUE]() == flow_none || el.flows[RED]() == flow_none) { // then 1 is enough
		return;
	}
	if (el.flows[BLUE]() == flow_both || el.flows[RED]() == flow_both || el.flows[BLUE]() == el.flows[RED]()) {
		el.capacity = 2;
	}
}

template<class G>
void setMinimalCapacities(G& g) {
	graph_traits<myTypes::MyGraph>::edge_iterator e, e_end;
	for (tie(e, e_end) = edges(g); e != e_end; ++e) {
		setMinimalCapacity(g[*e]);
	}
}

void setEdgeUsage(Edge_label& el, const Usage_E usage, const int fid) {
	if (el.flows[fid].usage != flow_none) {
		el.flows[fid] = flow_both;
	} else {
		el.flows[fid] = usage;
	}
	setMinimalCapacity(el);
}

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

template<class G = myTypes::MyGraph>
bool equal(const Path<G>& p1, const Path<G>& p2) {
	if (p1.size() != p2.size()) {
		return false;
	}
	for (int i = 0; i < p1.size(); ++i) {
		if (p1[i] != p2[i]) {
			return false;
		}
	}
	return true;
}

template<class G = myTypes::MyGraph>
bool prefixPath(const Path<G>& p1, const Path<G>& p2) {
	if (p1.size() > p2.size()) {
		return false;
	}
	for (int i = 0; i < p1.size(); ++i) {
		if (p1[i] != p2[i]) {
			return false;
		}
	}
	return true;
}

// set flow_none to all links
template<class G = myTypes::MyGraph>
void clearFlows(G &g) {
	typename graph_traits<G>::edge_iterator e, e_end;
	for (tie(e, e_end) = edges(g); e != e_end; ++e) {
		g[*e].flows[BLUE].usage = flow_none;
		g[*e].flows[RED].usage = flow_none;
	}
}

#endif /* FLOWUTIL_HPP_ */
