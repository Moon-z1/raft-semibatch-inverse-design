# RAFT半间歇共聚进料优化

本代码用于复现报告中的进料优化、代理模型比较和绘图结果。运行环境为 MATLAB R2024a，需要 Optimization Toolbox 和 Statistics and Machine Learning Toolbox。

## 使用方法

将 MATLAB 当前文件夹设为本目录，运行：

```matlab
reproduce_report                 % 检查已有结果记录
reproduce_report('figures')      % 将 report_figures 中的 FIG 导出为 PNG
reproduce_report('diagnostics')  % 检查第4.3节六条解析目标的补充结果
```

这三个入口读取已有结果，不重新搜索进料方案。需要重新计算时，使用下面各节的入口。

各目录中存在同名函数，不要对整个代码包使用 `addpath(genpath(...))`。进入对应的 `matlab` 目录后，先运行 `restoredefaultpath` 和 `addpath(pwd)`；若有 `core` 子目录，再运行 `addpath(fullfile(pwd,'core'))`。

## 各节计算入口

- **第2.1节：算法比较。** `raft_course_report/matlab` 中的 `run_ai_evidence`、`run_constraint_ai` 计算各方法结果；`raft_final/matlab` 中的 `run_final_evidence` 汇总配对比较，`repeat_timing` 测量运行时间。
- **第2.2节：已知目标恢复。** 在 `raft_revision/matlab` 中运行 `prepare_reachability`、`run_reachability`。
- **第2.3节：七条目标曲线。** 在 `raft_target_scan/matlab` 中依次运行 `run_scan_constrained`、`analyze_scan`、`verify_scan`，`make_figures` 用于绘图。
- **第3.1节：后段进料。** 在 `raft_late_feed/matlab` 中运行 `run_bridge('search')` 计算，`run_bridge('analyze')` 汇总，`verify_bridge` 核验结果。
- **第3.2节：初始组成。** 在 `raft_initial_composition/matlab` 中运行 `run_mayo`，用 `verify_mayo` 核验结果。
- **第3.3节：450次评价下的操作条件比较。** 在 `raft_feed_comparison/matlab` 中运行 `run_initial` 和 `run_segments`，再运行 `analyze_fair` 汇总。
- **第4.1节：误差代理预测。** 在 `raft_final/matlab` 中运行 `run_final_evidence`，`make_final_figures` 用于绘图。
- **第4.2—4.3节：轨迹代理。** 在 `raft_trajectory_pilot/matlab` 中运行 `restore_compact` 展开训练记录，`analyze_pilot` 分析，`audit_pilot` 核验；`run_pilot` 用于重新搜索。
- **第4.3节补充测试：六条解析目标。** 进入 `raft_target_feasibility`，分别运行 `run_original(101)`、`run_original(202)` 和 `run_original(303)` 搜索原操作条件下的方案，运行 `run_expanded(11)` 搜索扩展操作条件下的方案；`reproduce_diagnostics('replay')` 对最终方案重新积分，`reproduce_diagnostics('tables')` 导出误差汇总。
- **第4.4节：八条新增目标。** 在 `raft_next/matlab` 中运行 `run_extension`，再运行 `analyze_extension` 分析。
- **附录A：预算与泵速对照。** 入口为 `raft_step45/matlab/run_diagnostics`。
- **附录C：1500次评价对照。** 在 `raft_high_budget/matlab` 中运行 `run_resolution('search')`，用 `verify_resolution` 核验，`make_resolution_report_figure` 绘图。

## 数据说明

`outputs` 保存计算结果，`historical` 保存部分对照共用的初始方案。请保留目录结构，部分计算需要读取其他目录中的初始记录。搜索入口通常会跳过已有结果；若要从头重算，建议另复制一份代码，备份并移出对应结果文件，保留训练数据和共同初始方案。

`analysis.mat` 等分析文件由相应分析入口生成。MAT文件中的RMSE和最大偏差通常为无量纲原值，报告中乘100后以%表示。重新计算的耗时可能与原记录不同。

`report_figures` 保存MATLAB图文件，`apply_report_style` 设置坐标框、字体和图例格式。
