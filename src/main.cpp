#include <stdio.h>
#include <execinfo.h>
#include <gurobi_MIP.hpp>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <time.h>
#include <omp.h>

#include "testcases.hpp"
#include "Saeed_alg.hpp"
#include "gurobi_MIP.hpp"
#include<all_flowpairs.hpp>
#include<reporter.hpp>

#define TEST 0

using namespace std;

template<class G, typename V>
bool verify(G &g, const int rounds, const bool cyclic, vector<int> idx, ForEach_allocation<G> &forEach, const V &S,
		const V &T) {
	//					for(int i : idx) {
	//						Path<G> p = forEach.paths[i];
	//						copy(p.begin(), p.end(), ostream_iterator<int>(cout, ","));
	//						printf("\n");
	//					}
	// verify the result using the MIP solver
	clearFlows(g);
	addSP(g, BLUE, flow_old, forEach.paths[idx[0]], S, T);
	addSP(g, BLUE, flow_new, forEach.paths[idx[1]], S, T);
	addSP(g, RED, flow_old, forEach.paths[idx[2]], S, T);
	addSP(g, RED, flow_new, forEach.paths[idx[3]], S, T);
	print_network(g);
auto [feasible,rounds1]=runILP(g, S, T);
																												assert(feasible == !cyclic);
	if (feasible) {
		assert(rounds == rounds1);
	}

}

template<class G>
int maxRounds(G& g, Reporter &report) {
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
						//verify(g, rounds, cyclic, idx, forEach, *v1, *v2);
						max_all = rounds;
					}
					// report the result
					report << cyclic << '\t' << rounds << '\t' << nofBlocks << '\t'
					<< *v1 << ',' << *v2 << '\t' << forEach.paths.size()
					<< '\n';
				});
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

	std::ofstream all("data/overall.txt", std::ios_base::app);
	int maxr = 0, i = 0, tid;

#pragma omp parallel default(shared) private(tid,i)
	{
		std::ifstream filelist("data/input.txt");
		int tid = omp_get_thread_num();
		string path;
		while (filelist >> path) {
			if (path[0] == '#') {
				continue;
			}
			if (i++ % omp_get_num_threads() != omp_get_thread_num()) { // then not this thread's job
//			printf("skipping %s; i=%d, thread=%d\n",path.c_str(),i, tid);
				continue;
			}
			printf("Hello from thread %d, nthreads %d\n", tid, omp_get_num_threads());
			clock_t t1 = clock();
			Reporter rep(path);
			printf("\nLoading %s; thread=%d, i=%d\n", path.c_str(), tid, i);
//			exit(0);
			G g(0);
			auto [nodes, links]=loadFromFile(g, path.c_str());
			printf("%d links, %d nodes\n", links, nodes);

			const int rounds = maxRounds(g, rep);

			maxr = max(maxr, rounds);
			printf("max rounds=%d\n", rounds);
			clock_t time = (long double) (clock() - t1) / CLOCKS_PER_SEC;
			printf("elapsed: %u sec\n", time);
			all << rep.name << '\t' << nodes << '\t' << links << '\t' << maxr << '\t' << time << '\n';
			all.flush();
		}
	}
	all.close();
	clock_t elapsed_secs = double(clock() - t0) / CLOCKS_PER_SEC;
	printf("total time=%Lf sec\n", elapsed_secs);
	return 0;
}
#endif
