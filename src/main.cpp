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

#define TEST 0

using namespace std;

template<class G>
int maxRounds(G& g) {
	int max_all = 0;

	vector<vector<int>> bestPaths;
	VI<> v1, vend;

	for (tie(v1, vend) = vertices(g); v1 != vend; ++v1) {
		for (VI<> v2 = v1 + 1; v2 != vend; ++v2) {
			int max_st = 0;
			size_t x = 0;
			clock_t t = clock();
			const size_t X = 200000;
			ForEach_allocation<G> forEach(g, *v1, *v2, [&](size_t& left, const vector<int>& idx) {
				++x;
//				setMinimalCapacities(g);
					auto [cyclic, rounds, nofBlocks] = two_flows_update(g);
//				if(cyclic) {
//					return -1;
//				}
					PRINTF("rounds=%d, cyclic=%d, nofBlocks=%d, ST=%d,%d\n",
							rounds, cyclic, nofBlocks, *v1,*v2);
					if(x%X == 0) {
						float tpa = (double(clock() - t) / CLOCKS_PER_SEC)/X;
						printf("%f sec/case, left=%lu, wait=%.0fS, ST=(%d,%d), max rounds=%d\n",
								tpa, left, tpa*left,v1,v2, max_all);
						t = clock();
						fflush(stdout);
					}
					max_st = max(max_st, rounds);
					if(max_all < rounds) {
						printf("\nwinner! cyclic=%d, rounds=%d, paths=%d, ST=%d,%d\n",cyclic,rounds,forEach.paths.size(),*v1, *v2);
						max_all = rounds;
// report the result

//					for(int i : idx) {
//						Path<G> p = forEach.paths[i];
//						copy(p.begin(), p.end(), ostream_iterator<int>(cout, ","));
//						printf("\n");
//					}
						/*
						 clearFlows(g);
						 addSP(g, BLUE, flow_old, forEach.paths[idx[0]], *v1,*v2);
						 addSP(g, BLUE, flow_new, forEach.paths[idx[1]], *v1,*v2);
						 addSP(g, RED, flow_old, forEach.paths[idx[2]], *v1,*v2);
						 addSP(g, RED, flow_new, forEach.paths[idx[3]], *v1,*v2);
						 print_network(g);
						 auto [feasible,rounds1]=runILP(g, *v1,*v2);
						 assert(rounds==rounds1);
						 */
					}
				});
			if (max_st > 0) {
//				printf("ST=(%d,%d)=>%d ", *v1, *v2, max_st);
			}
		}
	}

	return max_all;
}
#if TEST
int main(int, char*[]) {
	init();
	using G = myTypes::MyGraph;
	std::ifstream filelist("data/input_test.txt");
	string path;
	while (filelist >> path) {
		if (path[0] == '#') {
			continue;
		}
		printf("\nLoading %s\n", path.c_str());
		G g(0);
		loadFromFile(g, path.c_str());
		print_network_forced(g);
		printf("%d links, %d nodes\n", num_edges(g), num_vertices(g));
		{
			auto [cyclic, rounds, nofBlocks] = two_flows_update(g);
			printf("Alg: rounds=%d, cyclic=%d, nofBlocks=%d\n", rounds, cyclic, nofBlocks);
		}
		{
			auto [feasible,rounds] = runILP(g, 0, 1);
			printf("ILP: rounds=%d, feasible=%d\n", rounds, feasible);
		}
	}
	return 0;
}
#else
int main(int, char*[]) {
	init();
	using G = myTypes::MyGraph;
	clock_t t0 = clock();

	//ofstream resultFile(Reporter::resultDir() + "overall.txt");
//	resultFile << "Name\t" << "Fail\t" << "Success\t"
//			<< "Ratio\t" << "MaxStackSize" << '\n';

	std::ifstream filelist("data/input.txt");
	string path;
	int maxr = 0;
	while (filelist >> path) {
		if (path[0] == '#') {
			continue;
		}
		printf("\nLoading %s\n", path.c_str());

		G g(0);
//		const char* f1 = "/Users/mahmoud/eclipse-workspace/2flowupdate/data/Garr201112.graphml_directed.gv";
		loadFromFile(g, path.c_str());
//		print_network(g);
		printf("%d links, %d nodes\n", num_edges(g), num_vertices(g));

		const int rounds = maxRounds(g);
		maxr = max(maxr, rounds);
		printf("max rounds=%d\n", rounds);
		clock_t t1 = clock();
		printf("elapsed: %d sec\n", double(clock() - t1) / CLOCKS_PER_SEC);
	}

	clock_t elapsed_secs = double(clock() - t0) / CLOCKS_PER_SEC;
	printf("total time=%d sec\n", elapsed_secs);
	return 0;
}
#endif
