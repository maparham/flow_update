#ifndef GUROBI_FLOWUPDATE_CPP_
#define GUROBI_FLOWUPDATE_CPP_

#include <functional>
#include "gurobi_c++.h"
#include "flowutil.hpp"

using namespace std;

class FlowUpdateLP {
private:
	const vector<FlowPair>& flowpairs;
	const int *fidmap;
	const int S, T;
	const int rounds, nofPairs;

	template<class Graph>
	void foreach_v(Graph& g, std::function<void(int)> f) {
		assert(num_vertices(g) > 0);
		typename graph_traits<Graph>::vertex_iterator vi, vend;
		for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
//			printf("foreach_v=%d\n", *vi);
			f(*vi);
		}
	}

	void foreach_r(std::function<void(int)> f) {
		for (int r = 0; r <= rounds; ++r) {
			f(r);
		}
	}

	void foreach_i_v(std::function<void(int, int)> f) {
		for (int i = 0; i < nofPairs; ++i) {
			foreach_v(flowpairs[i], [&](int v) {
				f(fidmap[i], v);
			});
		}
	}

	void foreach_flowpair(std::function<void(const FlowPair&, int)> f) {
		for (int i = 0; i < nofPairs; ++i) {
//			printf("foreach_flowpair: %d\n", num_vertices(flowpairs[i]));
			f(flowpairs[i], fidmap[i]);
		}
	}

	template<class Graph>
	void foreach_edge(Graph& g, std::function<void(int, int, Edge_label)> f) {
		typename graph_traits<Graph>::edge_iterator e, e_end;
		for (tie(e, e_end) = edges(g); e != e_end; ++e) {
			f(source(*e, g), target(*e, g), g[*e]);
		}
	}

	template<class Graph>
	void foreach_outedge(Graph& g, int v,
			std::function<void(int, int, Edge_label)> f) {
		typename graph_traits<Graph>::out_edge_iterator e, e_end;
		for (tie(e, e_end) = out_edges(v, g); e != e_end; ++e) {
			f(source(*e, g), target(*e, g), g[*e]);
		}
	}

	string vi(string varName, int v, int i) {
		ostringstream name;
		name << varName << "_" << v << "_" << i;
		return name.str();
	}

	string rvi(string varName, int r, int v, int i) {
		ostringstream name;
		name << varName << "_" << r << "_" << v << "_" << i;
		return name.str();
	}

	string ri(string varName, int r, int i) {
		ostringstream name;
		name << varName << "_" << r << "_" << i;
		return name.str();
	}

	string ruv(string varName, int r, int u, int v) {
		ostringstream name;
		name << varName << "_" << r << "_" << u << "_" << v;
		return name.str();
	}

	string ruvi(string varName, int r, int u, int v, int i) {
		ostringstream name;
		name << varName << "_" << r << "_" << u << "_" << v << "_" << i;
		return name.str();
	}

public:
	FlowUpdateLP(const vector<FlowPair>& flowpairs, const int fid[], const int n, const int s, const int t) :
			flowpairs(flowpairs), fidmap(fid), rounds(
					(n - 1) * flowpairs.size()), S(s), T(t), nofPairs(
					flowpairs.size()) {

	}

	void add_constraints(GRBModel& lp) {
		foreach_i_v([&](int i, int v) {
			// schedule constraints
				GRBLinExpr sum_r = 0;
				foreach_r([&](int r) {
							if (r > 0) {
								const GRBVar& x_r_vi = lp.addVar(0, 1, 0, GRB_BINARY, rvi("x", r, v, i));
								sum_r += x_r_vi;
							}
						});
				if(v!=T) {
					lp.addConstr(sum_r == 1, "schedule");
				}
			});

		// objective variable
		GRBVar R = lp.addVar(0, rounds, 1, GRB_INTEGER, "R");

		lp.update();

		foreach_flowpair([&](const FlowPair &flowpair, int i) {

			foreach_r([&](int r) {
//						printf("for r=%d\n",r);
				// R >= r.x^r_{v,fid}
				foreach_v(flowpair, [&](int v) {
//							printf("for v=%d\n",v);
							if (r > 0) {
								const GRBVar& x_r_vi = lp.getVarByName(rvi("x", r, v, i));
								lp.addConstr(R - r * x_r_vi >= 0, rvi("rounds", r, v, i));
							}
						});

				// add all y^r_{u,v,i} for all (u,v) \in P_i, then initialize
				foreach_edge(flowpair,
						[&](int u, int v, Edge_label lbl) {
//							printf("for e=(%d,%d)\n",u,v);
//							mylog<<"\nadding "<<ruvi("y",r,u,v,i);
							const GRBVar& y_r_uvi = lp.addVar(0, 1, 0, GRB_CONTINUOUS, ruvi("y", r, u, v, i));

							GRBLinExpr sum_r = 0;
// sum_{r'<=r}{x^r_vi}
							foreach_r([&](int r_) {
										if (r_ <= r && r_ > 0) {
											const GRBVar& x_rprime_vi = lp.getVarByName(rvi("x", r_, u, i));
											sum_r += x_rprime_vi;
										}
									});

							//  y for old path links
							if (lbl.flows[i].usage == flow_old ) {
								if (r == 0) {
									lp.addConstr(y_r_uvi == 1, ruvi("old",0,u,v,i));
								}
								else {
									lp.addConstr(y_r_uvi + sum_r == 1, ruvi("old",r,u,v,i));
								}
							}
							//  y for new path links
							if (lbl.flows[i].usage == flow_new ) {
								if (r == 0) {
									lp.addConstr(y_r_uvi == 0, ruvi("new",0,u,v,i));
								} else {
									lp.addConstr(y_r_uvi - sum_r == 0, ruvi("new",r,u,v,i));
								}
							}
							// unchanged links
							if(lbl.flows[i].usage==flow_both) {
								lp.addConstr(y_r_uvi == 1, ruvi("both",r,u,v,i));
							}
						});

				if (r == 0) {
					return; // skip the rest
				}

				lp.update();

				// fork^r_{vi}
//				GRBLinExpr m_r_i = 0;
				foreach_v(flowpair,
						[&](int v) {
							const GRBVar& fork_r_vi = lp.addVar(0, 1, 0, GRB_BINARY, rvi("fork", r, v, i));
							const GRBVar& x_r_vi = lp.getVarByName(rvi("x", r, v, i));
							lp.addConstr(fork_r_vi - isForkNode(v, flowpair, i) * x_r_vi == 0, rvi("forknode", r, v, i));
//							m_r_i += fork_r_vi;
						});

				lp.update();

				// f^r_{uvi}
				foreach_edge(flowpair,
						[&](int u, int v, Edge_label lbl) {

							const GRBVar& f_r_uvi = lp.addVar(0, 1, 0, GRB_CONTINUOUS, ruvi("f", r, u, v, i));
							GRBLinExpr sum1_uz = 0, sum2_uz = 0;
//							foreach_outedge(flowpair, u, [&](int u_, int z, Edge_label lbl) {
//										if (u != u_) {
//											exit(-2);
//										}
//										const GRBVar& y_r_uzi = lp.getVarByName(ruvi("y", r, u, z, i));
//										const GRBVar& y_rminus1_uzi = lp.getVarByName(ruvi("y", r - 1, u, z, i));
//										sum1_uz += y_r_uzi;
//										sum2_uz += y_rminus1_uzi;
//									});
							const GRBVar& y_r_uvi = lp.getVarByName(ruvi("y", r, u, v, i));
							const GRBVar& y_rminus1_uvi = lp.getVarByName(ruvi("y", r - 1, u, v, i));
							const GRBVar& fork_r_ui = lp.getVarByName(rvi("fork", r, u, i));

							// TODO
							lp.addConstr(f_r_uvi - y_r_uvi - fork_r_ui<= 0, ruvi("transientflow", r, u, v, i));
							lp.addConstr(f_r_uvi - y_rminus1_uvi - fork_r_ui<= 0, ruvi("transientflow", r - 1, u, v, i));

							const GRBVar& alpha_r_uvi = lp.addVar(0, 1, 0, GRB_BINARY, ruvi("alpha", r, u, v, i));

							// TODO
							lp.addConstr(alpha_r_uvi - 1 + y_r_uvi <= 0, ruvi("transientdemand", r, u, v, i));
							lp.addConstr(alpha_r_uvi - fork_r_ui <= 0, ruvi("transientdemand", r, u, v, i));
						});

				lp.update();

				// transient flow check
				foreach_v(flowpair,
						[&](int v) {
							GRBLinExpr rhs, inflows = 0, outflows = 0;
							GRBVar outflow1, outflow2;
							GRBVar inflow1, inflow2;
							int outcount = 0,incount=0;
							foreach_edge(flowpair, [&](int x, int y, Edge_label lbl) {
										const GRBVar& f_r_xyi = lp.getVarByName(ruvi("f", r, x, y, i));
										if (x == v) { // out-edge
											outflows += f_r_xyi;
											outflow2 = outflow1;
											outflow1 = f_r_xyi;
											++outcount;
										}
										if (y == v) { // in-edge
											inflows += f_r_xyi;
											inflow2 = inflow1;
											inflow1 = f_r_xyi;
											++incount;
										}
									});
//							if(r>1) {
//								return;
//							}
							const GRBVar& fork_r_vi = lp.getVarByName(rvi("fork", r, v, i));
							const GRBVar& join_r_vi = lp.addVar(0, 1, 0, GRB_BINARY, rvi("join", r, v, i));

							if (outcount == 2) { // a fork node ==> force equal out-flows
								lp.addConstr(outflow1 - fork_r_vi >= 0, rvi("outflow1", r, v, i));
								lp.addConstr(outflow2 - fork_r_vi >= 0, rvi("outflow2", r, v, i));
							}
							if (incount == 2) { // a join node ==> merge the in-flows
								lp.addConstr(inflow1 - join_r_vi >= 0, rvi("inflow1", r, v, i));
								lp.addConstr(inflow2 - join_r_vi >= 0, rvi("inflow2", r, v, i));

							} else { // prevent ordinary nodes from turning into a sink
								lp.addConstr( join_r_vi == 0, rvi("ordinarynode", r, v, i));
							}

							if (v == S) {
								rhs = fork_r_vi + 1;
							}
							else if (v == T) {
								rhs = -1-join_r_vi;
							}
							else {
								rhs = fork_r_vi - join_r_vi;
							}
							lp.addConstr(outflows - inflows - rhs == 0, rvi("transientcheck", r, v, i));
						});
			}); // r
}); // i

	// capacity constraints
			foreach_edge(
		 flowpairs[0].m_g, // main graph edges
		 [&](int u,int v, Edge_label lbl) {
		 foreach_r([&](int r) {
		 if(r<1) {
		 return;
		 }
		 GRBLinExpr demand_r_uv=0;
		 foreach_flowpair([&](FlowPair flowpair, int i) {
		 if(edgeExists(u,v,flowpair)) {
		 const GRBVar &f_r_uvi = lp.getVarByName(ruvi("f",r,u,v,i));
		 const GRBVar &alpha_r_uvi = lp.getVarByName(ruvi("alpha",r,u,v,i));

		 demand_r_uv += f_r_uvi;
//		 demand_r_uv -= alpha_r_uvi;
		 }
		 });
		 lp.addConstr(demand_r_uv <= lbl.capacity, ruv("capacity",r,u,v));
		 });
		 });
	}
};

pair<bool, int> runILP(myTypes::MyGraph& g, const int S, const int T) {

	vector<FlowPair> flowpairs;
	flowPairs(g, flowpairs);
	// TODO
	//flowpairs.pop_back();
	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.set(GRB_IntParam_OutputFlag, false);

		int n = num_vertices(g);
		FlowUpdateLP myLP(flowpairs, fid, n, S, T);
		myLP.add_constraints(model);
		model.optimize();

		int numvars = model.get(GRB_IntAttr_NumVars);
		mylog << "\nnumvars=" << numvars;
		auto vars = model.getVars();

		// report
		int optimstatus = model.get(GRB_IntAttr_Status);
		if (optimstatus == GRB_OPTIMAL) {
			for (int j = 0; j < numvars; j++) {
				GRBVar v = vars[j];
				if (v.get(GRB_DoubleAttr_X) != 0.0) {
					cout << endl << v.get(GRB_StringAttr_VarName) << " "
							<< v.get(GRB_DoubleAttr_X) << endl;
				}
			}
			double objval = model.get(GRB_DoubleAttr_ObjVal);
			mylog << "\nOptimal objective: " << objval << '\n';

//			model.write("debug.lp");
			return {true, objval};

		} else if (optimstatus == GRB_INFEASIBLE) {
			mylog << "\nModel is infeasible" << '\n';
//			model.computeIIS();
//			model.write("IISmodel.lp");
			return {false, -1};
		}

	} catch (GRBException e) {
		mylog << "\nError code = " << e.getErrorCode() << '\n';
		mylog << e.getMessage() << '\n';
	} catch (...) {
		mylog << "\nException during optimization" << '\n';
	}
}

#endif
