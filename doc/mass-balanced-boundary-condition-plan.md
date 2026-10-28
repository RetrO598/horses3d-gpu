# Mass-balanced boundary condition：迭代算法与开发计划

日期：2026-10-28。状态：设计提案，尚未实现或通过算例验证。

## 1. 范围与论文依据

参考 Gray et al. (2018), *Modeling Boundary Layer Ingestion Using a Coupled Aeropropulsive Analysis*, Journal of Aircraft 55, 1191–1199, DOI: 10.2514/1.C034601。本地 PDF 第 3–4 页（期刊页 1193–1194），III.B、III.D、式 (3) 及图 5 是本方案依据。

论文中 fan-face 是流出 CFD 域的静压出口，fan-exit 是流入 CFD 域的总压/总温入口。论文通过调整 fan-face 静压，使两截面质量流量一致：

$$
R_{p_s}^{FF}=\dot m^{FF}-\dot m^{FE}=0.
$$

论文完整算法的外层是 nonlinear Gauss–Seidel，内层 Broyden 同时处理静压及短舱直径两个残差。本次只保留质量匹配约束与静压未知量，固定几何、入口总压、总温和流向，采用带保护的标量割线求根。后文的稳定窗口、限幅、退回策略、软件结构及数值初值均是面向本项目的设计建议，不是论文给出的参数或原算法细节。

当前开发：稳态可压缩 Euler / laminar NS，随后验证 SA；一组可包括多个入口 patch 和多个出口 patch，但只有一个独立压力调节量。风扇图谱、效率、轴功率、FPR、pyCycle、几何调整、总压/总温反馈不进入当前迭代。

| 本文术语 | 定义与约定 |
|---|---|
| 入口组 I | CFD 域的 totalinflow，可包含多个 zone |
| 出口组 O | CFD 域的 outflow，可包含多个 zone |
| 质量流量 | 有符号净质量通量的面积积分，单位 kg/s |
| 外迭代 k | 每次固定压力、求取稳定 CFD 响应后的压力迭代 |
| 内迭代 n | 固定所有受控压力的一次完整 CFD 迭代 |

## 2. 先确定约束是否独立

目标应用是论文中的部件替代：两个截面之间不解析流动，靠质量约束保证流出的质量等于重新注入的质量。它们可处于同一个外流计算域中，该域还存在远场等其他边界。

如果 I 和 O 恰好是一个无质量源、无其他通流边界的封闭流道的全部入口和出口，则任意可行出口压力下的收敛稳态解本来就满足质量守恒。此时 R=0 不能额外决定出口压力，求根问题退化。控制器可能暂时消除储质量造成的不平衡，却不能据此选定物理工作点。这样的算例只适合验证积分和连续性，不适合验证压力求根；若需要选定工作点，应另给目标流量或部件压升关系。

一个压力自由度只能约束一个总量残差。多个 patch 的总流量相等，不意味着 patch 之间逐一相等，也不规定截面流量分布。首版不支持多个独立耦合组；后续若支持，必须检查残差 Jacobian 的独立性。

## 3. 流量定义与离散实现

法向统一使用 CFD 域外法向。对每个 zone z 定义：

$$
\Phi_z=\int_{\Gamma_z}\rho\mathbf u\cdot\mathbf n\,dA,
\qquad
\dot m_I=-\sum_{z\in I}\Phi_z,
\qquad
\dot m_O=\sum_{z\in O}\Phi_z,
\qquad R=\dot m_O-\dot m_I.
$$

因此正常运行时两个流量均为正。禁止对积分结果或积分点通量直接取绝对值后相减；局部回流必须保留其负贡献。另记录逆向流量占比，供异常判定使用。

### 3.1 两种积分各司其职

1. **控制器使用数值质量通量**：对组内边界上的连续方程数值通量积分，与残差装配所用通量一致。
2. **诊断使用状态流量**：保留现有 `mass-flow` 对内部状态 rho*u·n 的积分，与数值流量并列输出。

推荐离散式：

$$
\Phi_z^h=\sum_{f\in z}\sum_{i,j}w_iw_j\widehat F_{\rho,n}(Q^-,Q^+)J_f.
$$

当前 NS `computeBoundaryFlux` 中，`storage(1)%FStar` 在减去黏性通量后已经乘了面 Jacobian。若在这个位置读取质量分量，积分只乘两个求积权重，不能再乘一次 Jacobian。需核实法向与投影的符号约定，用均匀流和全域守恒检查确认。SA/RANS 路径单独核对，不能假设其存储位置相同。

数值通量的采样必须对应同一最终 Q、同一组压力：首版在需要评估控制器时，对完整步结束状态执行或复用一次新鲜的残差评估，再立即积分/保存边界质量通量；禁止读取 RK 最后一个中间 stage 或 Jacobian 扰动评估留下的旧缓存。可使用 state/pressure epoch 标记。现有 `compute time derivative after timestep` 机制可复用，但各积分器需确认是否真正对应最终状态。

### 3.2 无量纲与并行

内部压力保持当前约定 p*=p/p_ref；无量纲流量转 kg/s 的比例为 rho_ref*V_ref*L_ref^2，与项目参考量定义保持一致。输入和日志显式使用 Pa、kg/s，内部计算只用一致的无量纲量。

GPU 在设备端完成每个 zone 的积分和归约，只回传少量标量；所有 rank 按相同组顺序参与 MPI_Allreduce，无局部面的 rank 贡献 0。MPI 后 root 计算新压力并广播；新压力通过专用 setter 更新到已有 OpenACC 对象，等待其完成后才启动下一次边界 kernel。不能只修改 host 的 pExt，也不能每次重新 copyin 整个多态 BC 对象。

## 4. 控制变量与残差

首版支持一组内所有出口共享一个静压 p_s。将来若确需保留出口间压差，可扩展为 p_z=p_z,base+Delta p，共享一个 Delta p；不要为每个出口创建自由压力又只设置一个总量约束。

压力给定后，设固定边界参数下的收敛流场为 Q(p_s)，外层目标为：

$$
R(p_s)=\dot m_O(Q(p_s),p_s)-\dot m_I(Q(p_s),p_s)=0.
$$

**两个流量都必须重算。** 即使总压入口的总压/总温固定，其质量流量仍随流场和出口压力变化。不得将首次测得的入口流量永久当作目标值。

固定缩放 m_ref 使用配置的典型质量流量，或第一次稳定评估的 max((|m_I|+|m_O|)/2,m_floor)，整个一次求根过程不再改变。求根可用 r=R/m_ref；收敛判据采用：

$$
|R|\le m_{abs}+\varepsilon_m\max\left(\frac{|\dot m_I|+|\dot m_O|}{2},m_{floor}\right).
$$

m_abs 和 m_floor 均为有物理单位的配置值，按工况尺度设定。对于预计应有正向通流的组，两个流量接近零不能算正常成功；应输出 no-throughflow 状态。小压力增量也不能单独作为成功判据。

## 5. 推荐首版迭代：固定压力 CFD + 带保护的割线法

### 5.1 初始化

1. 解析组名、zone 名称、出口初压 p0 和正压范围 [p_min,p_max]。p0 是现有 outflow 的初始值。
2. 校验 I/O 非空、不重叠、没有重复 zone；入口类型为 totalinflow、出口类型为 outflow；压力范围合法、全局面积非零。单个 MPI rank 面数为零是正常情况。
3. 保持入口 p0/T0/direction 固定，初始化 Q；若 restart，先恢复已接受压力，再恢复或重置求根历史。
4. 在 p0 下完成固定压力内迭代，获得第一组有效响应 (p0,R0)。若已经满足全部收敛要求，结束，不进行不必要扰动。

### 5.2 固定压力内迭代与有效采样

一轮压力试探期间，所有 RK stage、隐式 Newton/Jacobian 评估、线性求解、FAS/coarse-grid 调用均使用同一压力。控制器只能在完整 CFD 步的安全位置工作，不能在 `FlowState` 回调内积分、通信或改变压力。

首版追求可诊断性：每个压力点达到配置的 CFD 残差阈值，并通过流量稳定性检查后才进入割线历史。稳定性在同一压力下比较两个连续窗口 A、B 的各组流量：

$$
e_{stat}=\frac{\max(|\bar m_{I,A}-\bar m_{I,B}|,|\bar m_{O,A}-\bar m_{O,B}|)}{m_{ref}}.
$$

还需检查窗口内波动幅度或标准差，避免周期波动恰好具有相同均值而误判。接受条件为 CFD 残差达标、e_stat 达标、窗口波动达标，且压力变化后至少完成最小等待步数。窗口不能跨压力点复用；噪声容限应显著小于目标质量误差，例如不超过其 10%。若当前工况持续非定常，首版返回 inner-not-steady，不默认把均值当稳态解。

内层收敛误差需要小于外层所需的流量精度。先用严格嵌套建立基线；通过测试后才能加入早期宽松、后期收紧的内层容差或每若干步的在线更新。

### 5.3 第二个压力点与斜率

在压力范围内选取一个小扰动 p1=p0+delta p_probe；若向上空间不足则向下扰动。固定 p1 重复内迭代，取得有效 R1。扰动要超过舍入和流量噪声水平，同时受最大压力步长约束。

$$
S_k=\frac{R_k-R_{k-1}}{p_k-p_{k-1}},\qquad
\Delta p_N=-\frac{R_k}{S_k}.
$$

斜率必须由**两个充分稳定的压力点**估计，而非固定 Q 下仅重算边界状态得到的瞬时导数。入口变化与外场反馈都应包含在 S 中。

若入口流量近似不变，通常提高静压会降低出口流量，因此 R>0 时倾向于升压。但本方案不硬编码 S<0：在耦合外流域中，入口也会响应，必须以实际稳定响应判定。

只有 |delta p| 超过压力分辨率且 |delta R| 超过流量噪声阈值时斜率才有效。无效时扩大探测（仍受限幅）或利用已有异号区间；不做接近零分母的除法。

### 5.4 压力提案、限幅与异号区间

无有效异号区间时，使用松弛割线步：

$$
p_{trial}=\Pi_{[p_{min},p_{max}]}
\left[p_k+\operatorname{clip}(\omega\Delta p_N,-\Delta p_{max},\Delta p_{max})\right].
$$

保存有效稳定评估点；若存在 R_a*R_b<0，建立 [p_a,p_b]。在异号区间内提议割线点；候选点落在区间外、过度靠近端点或不能提供充分区间收缩时，改用二分点。已有区间时最终试探点保持在该区间内；最大步长约束也需满足。二分兜底的收敛性质以 R 连续且内层落在同一解支为前提。

没有异号区间又没有有效斜率时，按有界探测扩展寻找响应和区间；不得假设上下限一定跨根。在达到配置范围、探测次数或外迭代预算后仍无根，明确报告 no-bracket / no-root-in-range / insensitive，保留最好可行解，不能宣告质量平衡完成。

### 5.5 试探点结果处理

- CFD 内迭代不稳定或出现负密度/负压：拒绝该点，不用于割线；恢复最近接受的 Q 和压力，缩小步长后重试。回滚状态必须包含积分器所需历史，不能只还原 Q 而留下不一致 BDF 历史。
- CFD 稳定且跨越 R=0：更新异号区间，用区间法继续。
- CFD 稳定但 |R| 明显恶化：保留该点作为有价值的区间/斜率信息，主解回退到较优可行点，减小松弛或改用区间法，避免直接沿同方向放大。
- 持续压力饱和、持续小斜率、大面积回流、压力边界失去调节能力：停止压力更新并输出明确失败状态。

允许累计有限次拒绝后退出；错误信息必须给出当前压力、范围、流量、CFD 残差及拒绝原因。

### 5.6 成功与退出

成功要求同时满足：CFD 残差达标、有效稳定采样、质量误差达标、流动方向/通流量有效。候选根保持同一压力，再完成一个独立确认窗口，通过后标记 balanced。确认期间不需要再次估计斜率或改变压力。

达到最大内迭代/外迭代数、用户停止或压力范围耗尽属于未完成，输出相应状态。现有“CFD 残差达标即退出”在初始化及主循环中的两处逻辑都必须修改：启用控制器时只有 CFD 与全部质量约束共同满足才能正常结束。压力改变后，旧残差和旧质量误差立即失效。

### 5.7 伪代码

```text
validate_group_and_configuration()
restore_or_initialize(Q, pressure, controller_state)
while outer_evaluations < max_outer:
    apply_pressure_to_host_and_device(pressure)
    clear_windows_and_mark_response_invalid()
    advance_CFD_with_fixed_pressure_until_valid_response_or_failure()
    if inner_failure:
        restore_last_accepted_full_solver_state()
        shrink_trial_or_stop()
        continue
    refresh_final_state_boundary_fluxes()
    allreduce_signed_zone_fluxes()
    R = sum_outward_flux(outlets) + sum_outward_flux(inlets)
    record_valid_pressure_response_and_update_bracket()
    if mass_tolerance_met and flow_valid:
        confirm_at_same_pressure_with_new_window()
        if confirmation_passes: return BALANCED
    pressure = safeguarded_secant_or_bracket_or_probe()
return EXPLICIT_FAILURE_STATUS
```

## 6. 初始调参建议

以下仅为开发起点，必须通过算例调整，不是论文验证值：

| 参数 | 建议起点 | 用途 |
|---|---|---|
| relative mass tolerance | 1e-4，之后验证 1e-5 | 外层质量误差 |
| relaxation omega | 0.3–0.5，建立可靠区间后可用 1 | 早期抑制过大割线步 |
| pressure probe | 特征绝对压力的 0.1%–0.5% | 估计稳定响应斜率 |
| max pressure step | 特征绝对压力的 1%–2% | 抑制压力突跳 |
| minimum settle steps | 20–100 个完整内迭代 | 只是下限，不能代替稳定判据 |
| window size | 每窗 20 个完整内迭代 | 在同一固定压力下检测稳定性 |
| max outer evaluations | 30–50，包含探测/失败试探 | 限制成本并保留明确状态 |

压力上下限应来自目标工况，首版要求显式输入，不能假设出口静压必须小于某个固定总压，也不能用 0 Pa 作为有效下限。绝对质量容差、m_floor、流量噪声阈值、CFD 内层阈值、最大内层步数都需配置或按明确尺度生成并打印。

## 7. 项目落点与模块划分

以下相对路径均以仓库根目录为基准，所列现有实现已阅读。

| 代码位置 | 当前能力 | 建议修改 |
|---|---|---|
| `Solver/src/libs/physics/common/bcs/TotalInflowBC.f90` | 总压/总温/流向与特征入口状态 | 当前参数固定；保留现有实现 |
| `Solver/src/libs/physics/common/bcs/OutflowBC.f90` | pExt；亚声速施加静压，超声速外推；OpenACC 对象驻留 | 添加受控压力 setter 和设备同步；增加逆流/可控面积诊断 |
| `Solver/src/libs/physics/common/bcs/BoundaryConditions.f90` | BC 构造、全局 BCs、设备数据初始化 | 提供类型安全的可调静压接口；控制器初始化后解析组关联 |
| `Solver/src/libs/monitors/SurfaceIntegrals.f90` | 内部状态质量流量积分；GPU reduction；MPI_Allreduce | 提取共享的轻量积分功能；数值通量积分独立命名，避免改变已有 monitor 语义 |
| `Solver/src/NavierStokesSolver/SpatialDiscretization.f90` | `computeBoundaryFlux` 组装数值通量 | 在最终状态采样路径暴露/保存质量通量，控制其有效期 |
| `Solver/src/NavierStokesSolverRANS/SpatialDiscretization.f90` | 独立 RANS 通量路径 | SA 阶段适配并验证等价积分，包含实际连续方程通量中的附加项 |
| `Solver/src/libs/timeintegrator/TimeIntegrator.f90` | 完整步结束、监测、两处稳态退出、SaveRestart | 安全采样与压力更新钩子；联合收敛；controller checkpoint |
| `Solver/src/libs/io/ReadInputFile.f90` 及相关配置代码 | 全局控制参数读取与模板 | 加入组配置、容差、压力范围及配置错误提示 |

建议新增 `MassBalanceController.f90`，放入 timeintegrator 库：负责组状态、内外迭代状态机、标量求根、日志、restart；通过低层接口访问 BC 和流量。不要让 BC 类反向依赖 monitors 或 timeintegrator，避免模块循环。若抽取积分公共模块，应放在能被上层使用、但不依赖控制器的位置。

首版复用 `totalinflow` 和 `outflow`，不增加第三种独立局部状态边界。mass balance 是若干既有边界的非局部约束，注册组即可表达。局部 FlowState 保持只计算外侧状态。

控制器数据建议包括：组名、I/O zone IDs、压力及界限、固定 m_ref、两个有效响应点、异号区间、窗口统计、等待步数、迭代计数、当前状态、最后接受 solver checkpoint 标识。配置接口先支持一个组，内部以可扩展数据类型表达。

配置语义草案（不是当前可运行 control 语法，实际语法需遵循现有 parser）：

```text
mass balance enabled = true
mass balance group = bypass
mass balance inlet zones = reinjection
mass balance outlet zones = extraction
mass balance pressure mode = shared
mass balance initial pressure = <Pa>
mass balance minimum pressure = <Pa>
mass balance maximum pressure = <Pa>
mass balance relative tolerance = 1.0e-4
mass balance absolute tolerance = <kg/s>
mass balance reference flow = <kg/s>
mass balance flow floor = <kg/s>
mass balance pressure probe = <Pa>
mass balance maximum pressure step = <Pa>
mass balance relaxation = 0.5
mass balance minimum settle steps = 50
mass balance window size = 20
mass balance maximum outer evaluations = 40
```

初始组压力与原 outflow pressure 字段若同时提供，应要求一致或明确定义唯一来源，禁止静默覆盖。组内入口总温和方向仍为必需输入，总压本身不能完整定义现有入口。

## 8. 并行、积分器与失败处理细节

- 初版只启用 steady-state 模式。真实非定常中一般有 dM/dt=m_in-m_out，不能默认把瞬时差压成零。统计均值控制作为独立后续功能设计。
- 当前 outflow 使用 |u_n/a| 判定亚声速，因此超声速出口静压不起作用；若可控面积太小或实测 |dR/dp| 很小，返回 insensitive。入口实现中也有 Mach 限制，需记录其触发情况，不能只靠出口诊断阻塞。
- 现有 outflow 对强回流不能直接视作可靠的总条件入口；持续回流比例超过配置阈值时返回 unsupported-flow-regime，并输出定位信息。
- 有限差分 Jacobian 和隐式内迭代必须冻结控制器。压力改变后，依赖边界状态的 Jacobian、预条件器或边界缓存应按现有策略失效并重建。
- FAS、多重网格、p-adaptation 和 FMG 不在首版支持范围；启用时显式拒绝此组合，随后逐项适配，不能让粗网格内迭代调用控制器。
- 全域质量守恒与本组平衡分别监测。全域守恒不能证明某一组已平衡；本组平衡也不能证明远场或其他开口守恒。
- 保证 MPI collective 顺序一致，失败决策也需全局一致，避免某个 rank 先退出而其他 rank 阻塞。
- restart 可先使用与 hsol 同名的带版本 sidecar。必须绑定迭代号、时间及组配置指纹，同步保存已接受 Q 与对应压力。缺少 sidecar 时不能假装精确续算：要求显式 reset-controller 模式，从指定初压重新稳定，清空割线/窗口历史。

建议日志每次评估输出：group、CFD iteration、outer evaluation、pressure(Pa)、m_I、m_O、R(kg/s)、相对误差、状态流量与数值流量之差、CFD residual、流量稳定指标、逆流比例、可控亚声速面积比例、斜率、松弛、区间端点、接受/拒绝原因。每次压力更新打印旧值和新值，最终结果区分 balanced 与各类 failure。

## 9. 分阶段开发及验收

| 阶段 | 开发内容 | 通过条件 |
|---|---|---|
| P0：流量观测 | 最终状态数值通量积分、组解析、符号/单位、独立日志；暂不调压 | 平面及曲面均匀流积分正确；多 patch 求和不重复；积分与离散全域连续性相符 |
| P1：压力通路 | outflow setter、device 同步、完整步钩子、联合退出条件、显式模式限制 | 单次压力改变确实改变下一步边界状态；CPU/GPU 均有效；所有子步保持固定压力 |
| P2：单组求根 | 内层稳定窗口、压力探测、保护割线、区间兜底、失败状态、回退 | 对已知可达工作点，从偏高/偏低初压均收敛；误差达到配置目标；无 NaN 或静默饱和 |
| P3：并行与续算 | MPI 归约/广播、全局失败决策、checkpoint/sidecar、恢复 | 串行/多 rank、CPU/GPU 在归约误差允许范围内一致；恢复后工作点一致；缺失 sidecar 行为明确 |
| P4：SA 与稳健性 | RANS 质量通量路径、SA 通量/状态兼容、阻塞/逆流、回归文档 | Euler、laminar NS、SA 各有通过算例；不支持的组合有可诊断错误；关闭控制器保持原行为 |

### 有物理意义的测试设计

1. **积分验证**：给定 rho、速度和任意法向面积，校验正负号、曲面 Jacobian、多个 patch 汇总、Pa/kg/s 转换；含一个 rank 没有组内面的 MPI 用例。
2. **连续性验证**：普通封闭管道，固定出口压力，用于校验数值流量与体积质量变化/离散残差一致；明确不拿它证明压力唯一性。
3. **控制器主算例**：两个无解析连接的流道/支路，共享同一受控组。支路 A 从固定外部入口流向受控静压出口；支路 B 从固定总条件入口流向独立固定外部出口。两支路流量没有被同一封闭流道连续性自动绑定。在可实现范围内匹配流量，建立非退化 R(p)。可先用准一维无黏关系给出参考，再用固定压力扫描核实 CFD 根的存在、连续性和唯一性。
4. **初压与容差扫描**：至少三个初压（根两侧及接近根），比较最终压力和流量；收紧目标误差后，内层精度也随之收紧，检查误差平台。
5. **失败用例**：压力范围不跨根、双支路均阻塞导致斜率接近零、逆流、内层不收敛、空组、重复 zone、非法类型及零流量。检查有界退出及错误状态，不只检查进程是否结束。
6. **并行/restart**：同一主算例使用不同 MPI 分区，CPU/GPU（在可用硬件上）及中途 restart，检查质量目标和最终压力在事先定义的数值误差内一致。
7. **关闭功能回归**：现有 totalinflow/outflow 算例禁用控制器时，边界结果与原路径一致。

主算例验收应同时报告最终 CFD 残差、组质量误差、全域守恒误差、压力历史、内外迭代成本。不能只凭两条 monitor 曲线看似重合判断成功。

## 10. 后续扩展边界

通过单组稳态基线后，可先评估降低内层求解成本、多个独立组的 Jacobian 耦合以及多重网格兼容。多个组若相互影响，独立标量控制器可能振荡，应考虑向量 Broyden 或有全局再检查的块迭代。

风扇耦合只保留未来接口：获取截面流量和质量加权总状态、更新入口总压/总温。未来需要另行定义部件方程、能量/功率约束及耦合未知量。本次固定入口参数，不计算或迭代这些输入；质量相等本身不构成风扇能量模型。
