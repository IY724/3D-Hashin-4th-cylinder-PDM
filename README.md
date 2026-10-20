# 3D-Hashin- PDM of the 4th type cylinder

**整案例说明文档（中文）：[readme.pdf](./readme.pdf)** —— 三个子程序各自对应的案例、建模与收敛注意事项、3D Hashin 四种判据式与损伤演化矩阵 C(d) 的完整展开、INP 用法。

**Full case explanation (Chinese): [readme.pdf](./readme.pdf)**


> **Status: research in progress.** This repository holds the original release only:
> the single-layer 3D Hashin UMAT in `hashin_self/`, the two full-vessel input files,
> and `readme.pdf`. The version under active revision is not published here, so this
> snapshot stays as it is. Please do not use it for engineering decisions.
>
> **状态：仍在研究中。** 本仓库只收录原始版本：`hashin_self/` 的单层 3D Hashin 子程序、
> 两份整瓶输入卡，以及 `readme.pdf`。持续修改的版本不在这里，此快照保持原样，
> 请勿直接用于工程定案。

> **实现范围说明：** 本库是个人摸索阶段整理的版本，子程序未考虑铺层角度转换，
> 退化方式按层合板处理。结果仅供学习与对比，不宜直接用于工程判断。
>
> **Scope note:** This is a personal exploratory implementation. The UMAT does not
> apply ply-angle transformation, and degradation follows a laminate-level model.
> Results are for learning and comparison only, not for engineering decisions.

## Contents

| Path | What it is |
|---|---|
| `hashin_self/hashin.for` | 437 行。单层 3D Hashin 判据 + 基于断裂能的线性软化（CDM）；剪切刚度 `DDSDDE(4,4)/(5,5)/(6,6)` 不乘 `(1−df)` |
| `hashin_self/hashin2.for` | 437 行。与 `hashin.for` 仅差三行：剪切刚度保留 `(1−df)` 折减 |
| `hashin_self/onlyhahsin.for` | 122 行。只按文献判据做刚度折减，不含断裂能演化 |
| `hashin_self/hashin_constant.for` | 184 行。固定刚度折减系数版本，系数由 `PROPS` 传入 |
| `Job-hashinnewcdm.inp` | 18.5 MB。整瓶输入卡，配 `hashin.for` 的能量法路线 |
| `Job-hashinnewconstant.inp` | 18.5 MB。整瓶输入卡，配 `hashin_constant.for` 的固定折减路线 |
| `readme.pdf` | 3.4 MB。中文说明文档：三个子程序各自对应的案例、建模与收敛注意事项、3D Hashin 四种判据式与损伤演化矩阵 C(d) 的完整展开、INP 用法 |

UMAT 子程序与 INP 的详细说明见上文 [readme.pdf](./readme.pdf)。

Because I am also a graduate student working on this topic, it’s hard to avoid having some questions.

if any question about this project,send related questions to the email Z25150070@s.upc.edu.cn
