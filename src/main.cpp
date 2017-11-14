#include <stdio.h>
#include <execinfo.h>
#include <gurobi_MIP.hpp>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>

#include "testcases.hpp"
//#include "Saeed_alg.hpp"
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

int main(int, char*[]) {
	signal(SIGSEGV, handler); // install our handler

	//two_flows_algorithm();

	myTypes::MyGraph g(0);
	vector<myTypes::MyGraph*> blocks;

	//generate/load graph
	paperExample(g);
	//singleEdge(g);
	//minimalExample(g);

	// extract pairs
	FlowEdgeFilter<myTypes::MyGraph> edgeFilter1(BLUE, g), edgeFilter2(RED, g);
	FlowVertexFilter<myTypes::MyGraph> vertexFilter1(BLUE, g), vertexFilter2(
			RED, g);

	FlowPair p_blue(g, edgeFilter1, vertexFilter1);
	FlowPair p_red(g, edgeFilter2, vertexFilter2);
	vector<FlowPair> flowpairs = { p_blue, p_red };
	int fid[] = { BLUE, RED };

	print_network(g);
	cout << "\npair1:\n";
	print_graph(flowpairs[0]);
	cout << "\npair2:\n";
	print_graph(flowpairs[1]);

	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

		int n = num_vertices(g);
		FlowUpdateLP myLP(flowpairs, fid, n, example.S, example.T);
		myLP.add_constraints(model);
		model.optimize();

		int numvars = model.get(GRB_IntAttr_NumVars);
		cout << "\nnumvars=" << numvars;
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
			cout << "\nOptimal objective: " << objval << endl;
		} else if (optimstatus == GRB_INFEASIBLE) {
			cout << "\nModel is infeasible" << endl;
			model.computeIIS();
			model.write("IISmodel.lp");
		}
		model.write("debug.lp");

	} catch (GRBException e) {
		cout << "\nError code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "\nException during optimization" << endl;
	}
	return 0;
}
