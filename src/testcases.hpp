/*
 * testcases.hpp
 *
 *  Created on: 16 Oct 2017
 *      Author: mahmoud
 */

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

void exampleNetwork(myTypes::MyGraph &g) {
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

	graph_traits<myTypes::MyGraph>::edge_iterator e_it, e_end;
	for (tie(e_it, e_end) = edges(g); e_it != e_end; ++e_it) {
		put(edge_weight, g, *e_it, 1);
	}
}

void example1(myTypes::MyGraph &g) {
	g = myTypes::MyGraph(4);
	Edge<myTypes::MyGraph> e;

	e = add_edge(0, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(0, 4, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(2, 1, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_both);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(2, 3, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(3, 1, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(4, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

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
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(s, 3, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(s, 4, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(2, 4, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(2, 5, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(2, 6, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);

	e = add_edge(3, 2, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(4, 5, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(5, 6, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(5, 7, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_old);

	e = add_edge(5, t, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_new);

	e = add_edge(6, 7, g).first;
	g[e].capacity = 1;
	g[e].flows[BLUE] = Flow(flow_new);
	g[e].flows[RED] = Flow(flow_old);

	e = add_edge(7, t, g).first;
	g[e].capacity = 1;
	g[e].flows[RED] = Flow(flow_old);
}

#endif /* TESTCASES_HPP_ */
