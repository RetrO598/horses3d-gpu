# NTNU 风洞空盘零压差算例

算例直接使用 `MESH/Wind_turbine_disk_empty_cgns_g2_mesh.h5` 的 HOPR HDF5 网格，
含 38,558 个二阶六面体单元。风洞尺寸为
11.15 × 2.71 × 1.80 m，盘径为 0.894 m，
空盘厚度为 0.05 m。`rotor_cylinder` 在网格的区域表中有标签，控制文件为此声明了类型。
该区域的 2,560 个双侧面仍按流体内部面连接，求解器报告该区域有 0 个边界面。

| 网格物理面 | HORSES3D 边界 | 设置 |
| --- | --- | --- |
| `inlet` | `Inflow` | x = 0 m，均匀来流沿 +x，速度 10 m/s |
| `outlet` | `Outflow` | x = 11.15 m，静压 101325 Pa |
| `disk_upstream` | `Outflow` | x = 3.635 m，流体离开上游计算域，静压 101325 Pa |
| `disk_downstream` | `Inflow` | x = 3.685 m，沿 +x 补入流体，给定状态与风洞入口相同 |
| `disk_side` | `FreeSlipWall` | 空盘圆柱侧面，自由滑移和绝热 |
| `sideA`、`sideB`、`top`、`bottom` | `FreeSlipWall` | 风洞四壁，沿用仓库风机教程设置 |

两盘面给定的静压均为 101325 Pa。盘后入口采用均匀来流参考态，
因此两盘面的瞬时质量流量分别由各自的边界条件决定。算例关闭 actuator line，
省略塔架和机舱的浸没边界模型，并用求解器默认的均匀 +x 来流初始化。

论文第 2.1 和第 3 节提供了 NTNU 风洞尺寸、盘径和来流速度，以及可压缩 NS、
Vreman LES、Lax–Friedrichs、BR1、Gauss–Legendre 和显式 RK3 的设置。
本初测采用 p = 2、CFL 0.5 和 20 步。以 1 m 为参考长度，
由论文的叶尖弦长和弦长雷诺数得到 ν ≈ 1.5015 × 10⁻⁵ m²/s
和 Re ≈ 6.66 × 10⁵。参考静压和温度取 HORSES3D 默认大气值，
101325 Pa、288.889 K，对应 Ma = 0.0293438。
前后盘面的独立进出口是本次空盘测试的设置。

## 运行

先在 `Solver/` 运行 `./configure`，生成 `SETUP/Makefile`、`RESULTS/` 和求解器链接。
构建 NS 求解器时设置 `WITH_HDF5=YES` 和 `HDF5_DIR`。
算例使用求解器的默认 ProblemFile 库，均匀来流由默认初始化例程设置。
从本算例目录运行：

```sh
./horses3d.ns wind_turbine_empty_disk.control > RESULTS/run.log 2>&1
```

`RESULTS/` 中的 `.surface` 文件记录风洞和空盘的质量流量以及两个盘面的平均静压。
质量流量按流体域外法向为正，盘前预期为正，盘后预期为负，盘侧预期接近零。
监测值取自各面的流体侧状态，使用求解器无量纲量。
本算例的质量流量换算因子为 `rho_ref * U_ref * L_ref^2`，约为 12.215 kg/s，
压力换算因子为 `rho_ref * U_ref^2`，约为 122.15 Pa。
