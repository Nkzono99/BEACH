title: 表面電荷更新の数値仕様

Lang: [日本語](SurfaceChargeNumerics.md) | [English](SurfaceChargeNumerics.en.md)

# 表面電荷更新の数値仕様

このページは、表面電荷の保存量、batch 末尾の更新順、浮遊導体の連立方程式、並列集約を調べるための
数値・実装リファレンスです。モデルを選ぶだけなら[表面はどう帯電するか](SurfaceModels.html)を参照してください。

## 保存する量と符号

`q_elem(i)` は要素 $i$ の総電荷 [C] です。P0 panel の場評価時だけ、要素面積 $A_i$ で割って
表面電荷密度へ変換します。

$$
\sigma_i=\frac{q_i}{A_i}
$$

マクロ粒子 $p$ が要素 $i$ へ吸収された場合の差分は

$$
\Delta q_i\mathrel{+}=q_p w_p
$$

です。電子は負、正イオンは正の電荷を堆積し、吸収後の粒子は追跡から除かれます。collision に使う
ordered triangle と片側 field trace の符号は同じ mesh geometry から作り、surface model は triangle の
winding を変更しません。

## batch 末尾の更新順

1. batch 開始時の `q_elem` から、粒子追跡中に固定する場を作る。
2. 各粒子の最初の mesh hit による電荷差分を thread ごとに集める。
3. 表面放出 source の反作用電荷を加える。
4. 有効なら `neutral_return` の global な再重み付けを適用する。
5. thread 間を加算し、MPI rank 間で global な電荷差分を集約する。
6. `q_elem` へ差分を一度だけ加える。
7. conductor object があれば、総電荷を保って等電位化する。
8. 更新前後の正味差分、統計、履歴を計算する。

同じ batch 中の粒子は手順 1 の場を共有します。不完全な collision query または photo-ray query で batch が
成立しない場合は、部分的な粒子配列や放出差分を確定反映しません。

## 浮遊導体の連立方程式

`surface_model="conductor"` の要素は `mesh_id` ごとに group を作ります。要素 $i$ が group $g(i)$ に
属するとき、未知の要素電荷 $q_j$ と group 電位 $V_g$ に対して

$$
\sum_j A_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed}
$$

を課します。単位総電荷を持つ source triangle $T_j$ の P0 potential coefficient は

$$
A_{ij}=\frac{1}{A_j}\int_{T_j}
\frac{1}{|\mathbf c_i-\mathbf y|}\,dA_{\mathbf y}
$$

です。$\mathbf c_i$ は target 要素重心、$A_j$ は source 要素面積です。自己項を含む解析 panel potential を
principal-value 側規約で評価し、$\phi_i^\mathrm{fixed}$ は non-conductor 電荷と一様外場が作る電位を
`k_coulomb` で割った量です。

各 group には総電荷保存も課します。

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}
$$

この dense square system を部分 pivot 付き Gauss 消去で解き、conductor 要素だけを置換します。object 間の
電荷と non-conductor 要素は変更しません。現行 model は要素重心 collocation と P0 triangle influence を使い、
`field_boundary.mode="free"` だけを受理します。

## OpenMP、MPI、再開

particle loop は thread ごとの配列へ吸収電荷を集め、loop 後に加算します。光電子 closure に必要な放出・帰還量と
local な電荷差分は MPI 全体で集約するため、全 rank が同じ global `q_elem` を保持します。conductor relaxation は
その同一 state へ各 rank で決定論的に適用します。

checkpoint は確定済みの `q_elem` と関連 state を保存します。再開互換性と必須 file は
[実行・再開する](Execution.html)および[入力パラメータ](Parameters.html)を参照してください。

## `tol_rel` の定義

conductor relaxation まで含む更新前後の差を

$$
\Delta\mathbf q=\mathbf q^\mathrm{after}-\mathbf q^\mathrm{before}
$$

とすると、監視値は

$$
\mathrm{tol\_rel\ metric}
=\frac{\|\Delta\mathbf q\|_2}{\max(\|\mathbf q^\mathrm{after}\|_2,q_\mathrm{floor})}
$$

です。`tol_rel` は output metric であり、現行実装の early-stop 条件ではありません。

## 実装とテスト

- particle absorption と batch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- surface model facade: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- floating-conductor solver: [`bem_surface_models_conductor.f90`](../src/physics/bem_surface_models_conductor.f90)
- batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)
- model regression: [`test_surface_models.f90`](../tests/fortran/test_surface_models.f90)
