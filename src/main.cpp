#include <stdio.h>
#include <execinfo.h>
#include <gurobi_MIP.hpp>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>

#include "testcases.hpp"
#include "Saeed_alg.hpp"
#include "gurobi_MIP.hpp"

using namespace std;

void handler(int sig) {
	void* array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}

void two_flows_algorithm_sim() {
	// the network graph
	myTypes::MyGraph g(0);

	randomNetwork randnet;
	int rounds = -1, nofBlocks = -1;
	double count = 0, acceptance = 0;
	bool feasible;
	srand(time(NULL));

	do {
		++count;
		g = myTypes::MyGraph(0);

//		exampleNetwork(g);
//		example_cyclic(g);
//		example3(g);
//		longDependency(g);
		StefanGraph sg(g, 4);
//		setVertexNames(g);

// random graph
//		int nof_vert = rand() % 20 + 6;
//		int nof_edges = nof_vert + rand() % (1 * nof_vert);
//		if (!randnet.generate(g, nof_vert, nof_edges)) {
//			continue;
//		}

		if (trivial(g)) {
			continue;
		}

		auto [cyclic, diameter_, nofBlocks_] = two_flows_update(g);
//		 feasible = !cyclic;
//		 rounds = diameter_;
		nofBlocks = nofBlocks_;

		tie(feasible, rounds) = runILP(g, 0, 1);

		if (feasible == cyclic) {
			printf("trouble!!\n");
			print_network(g);
			return;
		}

		mylog << "\nfeasible? " << feasible << " rounds=" << rounds
				<< "\n";
		if (feasible) {
			++acceptance;
		}

		if ((size_t) count % 1000 == 0)
			printf("%f\n", acceptance / count);

	} while ((rounds < 5 || !feasible));

	cout << "\nnofBlocks=" << nofBlocks << " diameter=" << rounds << " isDAG="
			<< feasible;
	print_network_forced(g);
	save_dot_file(
			"rounds" + to_string(rounds) + "DAG" + to_string(feasible)
					+ ".dot", g);
}

int main(int, char*[]) {
	signal(SIGSEGV, handler); // install our handler
	srand(time(NULL));

	myTypes::MyGraph g(0);

	//generate/load graph
	//paperExample(g);
	//singleEdge(g);
	//minimalExample(g);
	example4(g);
	//StefanGraph(g, 10);

	print_network(g);

//	runILP(g, 0, 1);
//	two_flows_algorithm_sim();
	OverlapFlows<> alg(g);
	alg.allocate();

	print_network(g);
	/*
	 cout << "\npair1:\n";
	 print_graph (flowpairs[0]);
	 cout << "\npair2:\n";
	 print_graph(flowpairs[1]);
	 */

	return 0;
}
