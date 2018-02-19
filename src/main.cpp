#include <stdio.h>
#include <execinfo.h>
#include <gurobi_MIP.hpp>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <time.h>

#include "testcases.hpp"
#include "Saeed_alg.hpp"
#include "gurobi_MIP.hpp"
#include<all_flowpairs.hpp>

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

template<class G>
int maxRounds(G& g) {
	int maxd = 0, x = 0;
	ForEach_allocation<G>(g, [&]() {
		++x;
//		printf("Allocation %d\n",x);
//			print_network_forced(g);
			/*
			 * 		setMinimalCapacities(g);
			 auto [feasible, rounds1] =
			 runILP(g, 0, 1);
			 if(feasible) {
			 maxd = max(maxd,rounds1);
			 printf("rounds=%d\n", rounds1);
			 }
			 */

			auto [cyclic, rounds2, nofBlocks] = two_flows_update(g);
			PRINTF("rounds=%d, cyclic=%d, nofBlocks=%d\n",
					rounds2, cyclic, nofBlocks);
			maxd = max(maxd,rounds2);
//			assert(feasible==!cyclic);
		});
	printf("#Allocations=%d\n", x);
	return maxd;
}

int main(int, char*[]) {
	signal(SIGSEGV, handler); // install our handler
	srand(time(NULL));
	using G = myTypes::MyGraph;
	G g(0);

	const char* f1 = "/Users/mahmoud/eclipse-workspace/2flowupdate/data/Garr201112.graphml_directed.gv";
	//loadFromFile(g, f1);
//	paperExample(g);
	print_network(g);
	PRINTF("%d links, %d nodes\n", num_edges(g), num_vertices(g));

	clock_t t0 = clock();

	const int max_rounds = maxRounds(g);

	printf("max diameter=%d\n", max_rounds);

//	runILP(g, 0, 1);
//	two_flows_algorithm_sim();

	/*
	 cout << "\npair1:\n";
	 print_graph (flowpairs[0]);
	 cout << "\npair2:\n";
	 print_graph(flowpairs[1]);
	 */

	clock_t elapsed_secs = double(clock() - t0) / CLOCKS_PER_SEC;
	printf("total time=%d sec\n", elapsed_secs);
	return 0;
}
