\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  R
Subject To
 schedule: x_1_0_1 + x_2_0_1 = 1
 schedule: x_1_2_1 + x_2_2_1 = 1
 old_0_0_1_1: y_0_0_1_1 = 1
 new_0_0_2_1: y_0_0_2_1 = 0
 new_0_2_1_1: y_0_2_1_1 = 0
 rounds_1_0_1: - x_1_0_1 + R >= 0
 rounds_1_1_1: - x_1_1_1 + R >= 0
 rounds_1_2_1: - x_1_2_1 + R >= 0
 old_1_0_1_1: x_1_0_1 + y_1_0_1_1 = 1
 new_1_0_2_1: - x_1_0_1 + y_1_0_2_1 = 0
 new_1_2_1_1: - x_1_2_1 + y_1_2_1_1 = 0
 transientflow_1_0_1_1: - 3 y_1_0_1_1 - 3 y_1_0_2_1 + f_1_0_1_1 <= 0
 transientflow_0_0_1_1: - 3 y_0_0_1_1 - 3 y_0_0_2_1 + f_1_0_1_1 <= 0
 transientflow_1_0_2_1: - 3 y_1_0_1_1 - 3 y_1_0_2_1 + f_1_0_2_1 <= 0
 transientflow_0_0_2_1: - 3 y_0_0_1_1 - 3 y_0_0_2_1 + f_1_0_2_1 <= 0
 transientflow_1_2_1_1: - 3 y_1_2_1_1 + f_1_2_1_1 <= 0
 transientflow_0_2_1_1: - 3 y_0_2_1_1 + f_1_2_1_1 <= 0
 forknode_1_0_1: - x_1_0_1 + fork_1_0_1 = 0
 forknode_1_1_1: fork_1_1_1 = 0
 forknode_1_2_1: fork_1_2_1 = 0
 outflow1_1_0_1: f_1_0_2_1 - fork_1_0_1 >= 0
 outflow2_1_0_1: f_1_0_1_1 - fork_1_0_1 >= 0
 ordinarynode_1_0_1: join_1_0_1 = 0
 transientcheck_1_0_1: f_1_0_1_1 + f_1_0_2_1 - fork_1_0_1 = 1
 inflow1_1_1_1: f_1_2_1_1 - join_1_1_1 >= 0
 inflow2_1_1_1: f_1_0_1_1 - join_1_1_1 >= 0
 transientcheck_1_1_1: - f_1_0_1_1 - f_1_2_1_1 = -1
 ordinarynode_1_2_1: join_1_2_1 = 0
 transientcheck_1_2_1: - f_1_0_2_1 + f_1_2_1_1 - fork_1_2_1 + join_1_2_1
   = 0
 rounds_2_0_1: - 2 x_2_0_1 + R >= 0
 rounds_2_1_1: - 2 x_2_1_1 + R >= 0
 rounds_2_2_1: - 2 x_2_2_1 + R >= 0
 old_2_0_1_1: x_1_0_1 + x_2_0_1 + y_2_0_1_1 = 1
 new_2_0_2_1: - x_1_0_1 - x_2_0_1 + y_2_0_2_1 = 0
 new_2_2_1_1: - x_1_2_1 - x_2_2_1 + y_2_2_1_1 = 0
 transientflow_2_0_1_1: - 3 y_2_0_1_1 - 3 y_2_0_2_1 + f_2_0_1_1 <= 0
 transientflow_1_0_1_1: - 3 y_1_0_1_1 - 3 y_1_0_2_1 + f_2_0_1_1 <= 0
 transientflow_2_0_2_1: - 3 y_2_0_1_1 - 3 y_2_0_2_1 + f_2_0_2_1 <= 0
 transientflow_1_0_2_1: - 3 y_1_0_1_1 - 3 y_1_0_2_1 + f_2_0_2_1 <= 0
 transientflow_2_2_1_1: - 3 y_2_2_1_1 + f_2_2_1_1 <= 0
 transientflow_1_2_1_1: - 3 y_1_2_1_1 + f_2_2_1_1 <= 0
 forknode_2_0_1: - x_2_0_1 + fork_2_0_1 = 0
 forknode_2_1_1: fork_2_1_1 = 0
 forknode_2_2_1: fork_2_2_1 = 0
 outflow1_2_0_1: f_2_0_2_1 - fork_2_0_1 >= 0
 outflow2_2_0_1: f_2_0_1_1 - fork_2_0_1 >= 0
 ordinarynode_2_0_1: join_2_0_1 = 0
 transientcheck_2_0_1: f_2_0_1_1 + f_2_0_2_1 - fork_2_0_1 = 1
 inflow1_2_1_1: f_2_2_1_1 - join_2_1_1 >= 0
 inflow2_2_1_1: f_2_0_1_1 - join_2_1_1 >= 0
 transientcheck_2_1_1: - f_2_0_1_1 - f_2_2_1_1 = -1
 ordinarynode_2_2_1: join_2_2_1 = 0
 transientcheck_2_2_1: - f_2_0_2_1 + f_2_2_1_1 - fork_2_2_1 + join_2_2_1
   = 0
Bounds
 R <= 2
 y_0_0_1_1 <= 1
 y_0_0_2_1 <= 1
 y_0_2_1_1 <= 1
 y_1_0_1_1 <= 1
 y_1_0_2_1 <= 1
 y_1_2_1_1 <= 1
 f_1_0_1_1 <= 3
 f_1_0_2_1 <= 3
 f_1_2_1_1 <= 3
 y_2_0_1_1 <= 1
 y_2_0_2_1 <= 1
 y_2_2_1_1 <= 1
 f_2_0_1_1 <= 3
 f_2_0_2_1 <= 3
 f_2_2_1_1 <= 3
Binaries
 x_1_0_1 x_2_0_1 x_1_1_1 x_2_1_1 x_1_2_1 x_2_2_1 fork_1_0_1 fork_1_1_1
 fork_1_2_1 join_1_0_1 join_1_1_1 join_1_2_1 fork_2_0_1 fork_2_1_1
 fork_2_2_1 join_2_0_1 join_2_1_1 join_2_2_1
Generals
 R
End
