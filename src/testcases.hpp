#ifndef TESTCASES_HPP_
#define TESTCASES_HPP_

#include "helpers.hpp"

// Make convenient labels for the vertices
struct Example {
	int S;
	int T;
	int U;
	int V;
	int W;
	int numNodes = 5;
} example;

void setWeights(myTypes::MyGraph &g) {
	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

void paperExample(myTypes::MyGraph &g) {
	put(vertex_name, g, example.S = add_vertex(g), "0");
	put(vertex_name, g, example.T = add_vertex(g), "1");
	put(vertex_name, g, example.W = add_vertex(g), "2");
	put(vertex_name, g, example.U = add_vertex(g), "3");
	put(vertex_name, g, example.V = add_vertex(g), "4");

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
}

void minimalExample(myTypes::MyGraph &g) {
	example.S = add_vertex(g);
	example.T = add_vertex(g);
	example.U = add_vertex(g);
	Edge<myTypes::MyGraph> e;
	e = add_edge(example.S, example.T, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(example.S, example.U, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(example.U, example.T, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	setWeights(g);
}

void singleEdge(myTypes::MyGraph &g) {
	example.S = add_vertex(g);
	example.T = add_vertex(g);
	Edge<myTypes::MyGraph> e;
	e = add_edge(example.S, example.T, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_both);
	g[e].flows[BLUE] = Flow(flow_both);

	setWeights(g);
}

void example_cyclic(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(4);
	Edge<myTypes::MyGraph> e;

	e = add_edge(0, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(0, 3, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

//	e = add_edge(2, 3, g).first;
//	g[e].capacity = 1;
//	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(2, 1, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_both);
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(3, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(3, 1, g).first;
	g[e].capacity = 1;
//	g[e].flows[RED] = Flow(flow_new);
	g[e].flows[BLUE] = Flow(flow_new);

	setWeights(g);
}

void example1(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(6);
	Edge<myTypes::MyGraph> e;
	const int S = 0;
	const int T = 1;
	e = add_edge(S, 2, g).first;
	g[e].capacity = 2;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_both);

	e = add_edge(S, 3, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(2, 4, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(2, 5, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(3, 1, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(4, T, g).first;
	g[e].capacity = 2;
	g[e].flows[BLUE] = Flow(flow_both);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(5, T, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

void example2(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(4);
	Edge<myTypes::MyGraph> e;
	const int S = 0;
	const int T = 1;
	e = add_edge(S, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

//	e = add_edge(S, 3, g).first;
//	g[e].capacity = 2;
//	g[e].flows[BLUE] = Flow(flow_old);
//	g[e].flows[RED] = Flow(flow_both);

	e = add_edge(2, T, g).first;
	g[e].capacity = 2;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(S, T, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

// only feasible with batch update support
void example3(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(4);
	Edge<myTypes::MyGraph> e;
	const int S = 0;
	const int T = 1;
	e = add_edge(S, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(S, 3, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(2, T, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(3, T, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

// for auto path allocation
void example4(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(4);
	Edge<myTypes::MyGraph> e;
	const int S = 0;
	const int T = 1;
	e = add_edge(S, 2, g).first;
	e = add_edge(S, 3, g).first;
	e = add_edge(2, T, g).first;
	e = add_edge(3, T, g).first;

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

void longDependency(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(0);
	Edge<myTypes::MyGraph> e;
	int s = 0, t = 1;

	e = add_edge(s, 2, g).first;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(2, 3, g).first;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(3, 5, g).first;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(5, 7, g).first;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(7, 9, g).first;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(9, 11, g).first;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(11, 13, g).first;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(13, t, g).first;
	g[e].flows[BLUE] = Flow(flow_old);

	// level1
	e = add_edge(s, 3, g).first;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(s, 4, g).first;
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(4, 5, g).first;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(3, 4, g).first;
	g[e].flows[BLUE] = Flow(flow_both);

	// level2
	e = add_edge(4, 6, g).first;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(6, 7, g).first;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(5, 6, g).first;
	g[e].flows[RED] = Flow(flow_both);

	// level3
	e = add_edge(6, 8, g).first;
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(8, 9, g).first;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(7, 8, g).first;
	g[e].flows[BLUE] = Flow(flow_both);

	// level4
	e = add_edge(8, 10, g).first;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(10, 11, g).first;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(9, 10, g).first;
	g[e].flows[RED] = Flow(flow_both);

	// level5
	e = add_edge(10, 12, g).first;
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(12, 13, g).first;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(11, 12, g).first;
	g[e].flows[BLUE] = Flow(flow_both);

	// level6
	e = add_edge(12, 14, g).first;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(14, t, g).first;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_both);

	e = add_edge(13, 14, g).first;
	g[e].flows[RED] = Flow(flow_both);

	setMinimalCapacities(g);
}

class StefanGraph {
	using G = myTypes::MyGraph;

	void addRandomPath(const int fid, const Usage_E usage) {
		int r, jump;
		Vertex<G> next, currNode = s;
		Edge<G> e;
		enum {
			UPPER, LOWER
		} current = UPPER;
		for (int i = -1; i < k - 1; ++i) {
			r = rand() % 2;
			jump = (currNode == s) ? i + 1 : i + 1 + rand() % (k - i - 1);

			if (r == UPPER) {
				next = (current == LOWER) ? upper[jump] : upper[i + 1];
				//mylog << "next on upper:" << get(vertex_name, g, next);
				current = UPPER;

			} else { // r==LOWER
				next = (current == UPPER) ? lower[jump] : lower[i + 1];
				//mylog << "next on lower:" << get(vertex_name, g, next);
				current = LOWER;
			}
			auto p = edge(currNode, next, g);
			e = (p.second) ? p.first : add_edge(currNode, next, g).first; // stupid BGL doesn't provide a better way

			currNode = next;
			setEdgeUsage(g[e], usage, fid);
		}
		// finish off at t
		auto p = edge(currNode, t, g);
		e = (p.second) ? p.first : add_edge(currNode, t, g).first;
		setEdgeUsage(g[e], usage, fid);
	}

public:
	StefanGraph(G &g, int length) :
			g(g), k(length) {

		upper.resize(k);
		lower.resize(k);
		put(vertex_name, g, example.S = s = add_vertex(g), "s");
		put(vertex_name, g, example.S = t = add_vertex(g), "t");

		for (int i = 0; i < k; ++i) {
			upper[i] = add_vertex(g);
			lower[i] = add_vertex(g);
			put(vertex_name, g, upper[i], "U" + to_string(i));
			put(vertex_name, g, lower[i], "L" + to_string(i));
		}

		addRandomPath(BLUE, flow_old);
		addRandomPath(BLUE, flow_new);
		addRandomPath(RED, flow_old);
		addRandomPath(RED, flow_new);
	}
	G &g;
	Vertex<G> s, t;
	int k;
	vector<Vertex<G>> upper, lower;
};

#endif /* TESTCASES_HPP_ */
