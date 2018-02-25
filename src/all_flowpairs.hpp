/*
 * all_flowpairs.hpp
 *
 *  Created on: 16 Feb 2018
 *      Author: mahmoud
 */

#ifndef ALL_FLOWPAIRS_HPP_
#define ALL_FLOWPAIRS_HPP_

#define DEBUG 0

#if DEBUG
#define PRINTF printf
#else
#define PRINTF(format, args...) ((void)0)
#endif

template<class G = myTypes::MyGraph>
class ForEach_allocation {
	const int S, T;
	G &g;

	struct Path_Enum: public default_dfs_visitor {
		const Vertex<G> S, T;
		map<int, vector<Path<G>>> &subpathMap;
		Path_Enum(G &g, map<int, vector<Path<G>>> &subpathMap, const int S, const int T) :
				g(g), subpathMap(subpathMap), S(S), T(T) {
		}
		void discover_vertex(const Vertex<G> &v, const G &g) {
//			PRINTF("discover_vertex %d\n", v);
			//		put(vertex_distance, g, v, 0); // init
			if (v == T) {
//				PRINTF("T visited, paths.size=%d\n", paths.size());
				subpathMap[T].push_back(Path<G>( { T }));
			}
		}
		void finish_edge(const Edge<G> &e, const G &g) {
			const Vertex<G> parent = source(e, g);
			const Vertex<G> child = target(e, g);
//			PRINTF("finish_edge (%d,%d)\n", parent, child);
			vector<Path<G>> &parentPaths = subpathMap[parent];
			const vector<Path<G>> &childPaths = subpathMap[child];

			for (Path<G> p : childPaths) {
				if (inPath(parent, p)) {
//					printf("%d already included\n",parent);
					continue;
				}
				assert(p.back() == T);
				Path<G> pp( { parent });
				pp.insert(pp.begin() + 1, p.begin(), p.end());
				parentPaths.push_back(pp);
			}
//			PRINTF("parentPaths=%d\n", parentPaths.size());
		}
		void finish_vertex(const Vertex<G> &v, const G &g) {
//			PRINTF("finish_vertex %d\n", v);
			path.pop_back();
		}
		G &g;
		Path<G> path;
		vector<Path<G>> paths;
	};
	public:
	vector<Path<G>> paths;
	ForEach_allocation(G &g, const int S, const int T, std::function<void(size_t&, const vector<int>)> fn) :
			g(g), S(S), T(T) {
		map<int, vector<Path<G>> > pmap;
		Path_Enum pe(g, pmap, S, T);
		depth_first_search(g, visitor(pe)); // compute all S-T paths
		paths = pmap[S]; // all S-T paths

		for (int i = 0; i < paths.size(); ++i) {
			printPath(paths[i]);
		}
//		printf("paths.size=%d ", paths.size());
		if (paths.size() < 3) {
			return;
		}
		Combination comb1(paths.size(), 2);
		Combination comb2(paths.size(), 2);

		size_t num = pow(comb1(), 2);
//		printf(" #allocations=%lu\n", num);
		do {
			Combination comb2 = comb1;
			while (comb2.next()) { // initial +1 offset
				clearFlows(g);
				vector<int> idx = comb1.toVector(comb2); // join the vectors

#if DEBUG
				PRINTF("idx=");
				copy(idx.begin(), idx.end(), ostream_iterator<int>(cout, ","));
				PRINTF("\n");
#endif
//						assert(idx[0]!=idx[2] || idx[1]!=idx[3]);
				addSP(g, BLUE, flow_old, paths[idx[0]], S, T);
				addSP(g, BLUE, flow_new, paths[idx[1]], S, T);
				addSP(g, RED, flow_old, paths[idx[2]], S, T);
				addSP(g, RED, flow_new, paths[idx[3]], S, T);
				fn(--num, idx);
//				if (idx == vector<int>( { 0, 4, 2, 0 })) {
//					exit(0);
//				}
			}
		} while (comb1.next());
	}
};

#undef DEBUG
#endif /* ALL_FLOWPAIRS_HPP_ */
