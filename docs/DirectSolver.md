title: Direct solver

Lang: [日本語](DirectSolver.md) | [English](DirectSolver.en.md)

# Direct solver

Direct solverは、各評価点で全source要素の寄与を加算します。`triangle_p0`のsource離散化を近似せずに
評価するため、小規模計算に使えるだけでなく、TreecodeやFMMの近似誤差を測る基準にもなります。

| 特性 | 内容 |
| --- | --- |
| 計算量 | source数$N$、評価点数$M$に対して$O(MN)$ |
| 要素kernel | triangle P0（固定） |
| 通常の場境界 | free |
| 幾何plan | 不要。評価時に現在の`q_elem`を読む |

```toml
[sim]
field_solver = "direct"
field_bc_mode = "free"
```

<a id="triangle-p0"></a>

## 面電荷を三角形上で積分する

BEACHは要素面積$A_i$上に一定面密度$\sigma_i=q_i/A_i$を置きます。`triangle_p0`は暗黙の
唯一の要素kernelであり、設定で選択しません。

$$
\phi_i(\mathbf{r})=
k_c\int_{T_i}\frac{\sigma_i}{\lVert\mathbf{r}-\mathbf{r}'\rVert}\,dS',
\qquad
\mathbf{E}_i(\mathbf{r})=
k_c\int_{T_i}\sigma_i
\frac{\mathbf{r}-\mathbf{r}'}{\lVert\mathbf{r}-\mathbf{r}'\rVert^3}\,dS'.
$$

BEACHは三角形の辺に沿う対数項と立体角を使う解析kernelでこの積分を評価します。
遠方と近傍の寄与を同じpanel kernelで評価します。

面上の電位は連続で、自己電位にはprincipal valueを使います。電場の法線成分には面電荷によるjumpがあり、
真空側の極限は

$$
\mathbf{E}^{\pm}=\mathbf{E}^{\mathrm{PV}}
\pm\frac{\sigma_i}{2\epsilon_0}\mathbf{n}_i
$$

となります。この符号を決めるため、各要素に`normal_plus`または`normal_minus`のvacuum sideが必要です。
OBJ入力では`[mesh].surface_side`、templateでは各`[[mesh.templates]].surface_side`で指定します。
`outward_closed`は、向きが整合した閉じたtwo-manifoldにだけ使用できます。

現行要件は次のとおりです。

- 有限で面積が正の三角形
- 全要素でvacuum sideが解決済み
- solverは`direct`、`treecode`、`fmm`、または`direct`/`fmm`を選ぶ`auto`

旧`[field]` tableと`sim.softening`は削除済みです。残した設定はunknown table / keyとして停止します。

設定例は[`examples/panel_direct.toml`](../examples/panel_direct.toml)にあります。

<a id="要素中心の電位出力"></a>

## 要素中心電位では自己項を加える

`potential_history.csv`などの要素中心電位は、自要素を含む解析panel積分を使います。面上で有限な
triangle P0自己電位をそのまま含み、別の経験的自己係数は加えません。

要素中心を全要素について評価するため、Directの電位履歴は$O(N^2)$です。履歴を細かく保存すると、
粒子位置の場評価とは別に大きなコストになることがあります。

## commit済み電荷を直接評価する

内部座標は`field_normalization`で正規化され、電場と電位はSIへ戻されます。[<sup>1</sup>](FieldSolvers.html#長さの正規化)

Directは固定geometryのtreeや展開係数を持たず、各評価で現在の`mesh%q_elem`を直接読みます。batch中は
参照する要素電荷を固定し、batch末尾にcommitした電荷を次batchの場へ反映する更新順序は他のsolverと同じです。

## Directを基準に誤差を切り分ける

Directにはsolverによる近似誤差がありませんが、次の誤差は残ります。

- surface meshの解像度と形状誤差
- 粒子の時間刻みと衝突判定誤差
- periodic/outer modelを含む場合のmodel化誤差

Directとの差が小さければ、FMMやTreecodeの近似誤差が小さいことは確認できます。ただし、それだけでは
物理解全体の収束を保証できません。source meshと粒子時間刻みも独立に細分化してください。

## Code reference

- Directの電場・電位と自己項: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- triangle P0解析kernel: [`bem_panel_kernel.f90`](../src/physics/panel/bem_panel_kernel.f90)
- panel geometry: [`bem_panel_geometry.f90`](../src/physics/panel/bem_panel_geometry.f90)
- 設定の検証: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)
- 主な回帰テスト: [`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90)、[`test_panel_kernel.f90`](../tests/fortran/test_panel_kernel.f90)
