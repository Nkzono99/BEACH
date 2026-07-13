title: Direct solver

Lang: [日本語](DirectSolver.md) | [English](DirectSolver.en.md)

# Direct solver

Direct solverは、各評価点で全source要素の寄与を加算します。選択したsource離散化を近似せずに評価するため、
小規模計算に使えるだけでなく、TreecodeやFMMの近似誤差を測る基準にもなります。結果を大きく左右するのは、
要素電荷を重心の点電荷として扱うか、三角形上の面電荷として積分するかという選択です。

| 特性 | 内容 |
| --- | --- |
| 計算量 | source数$N$、評価点数$M$に対して$O(MN)$ |
| 対応kernel | point、triangle P0 |
| 通常の場境界 | free |
| 幾何plan | 不要。評価時に現在の`q_elem`を読む |

```toml
[sim]
field_solver = "direct"
field_bc_mode = "free"
```

## 要素電荷を重心の点電荷として扱う

`field.element_kernel="point"`では、要素$i$の総電荷$q_i$を三角形重心$\mathbf{c}_i$に置きます。
評価点$\mathbf{r}$の電場と電位は

$$
\mathbf{E}(\mathbf{r})=
k_c\sum_{i=1}^{N}q_i
\frac{\mathbf{r}-\mathbf{c}_i}
{\left(\lVert\mathbf{r}-\mathbf{c}_i\rVert^2+\epsilon^2\right)^{3/2}},
$$

$$
\phi(\mathbf{r})=
k_c\sum_{i=1}^{N}
\frac{q_i}{\sqrt{\lVert\mathbf{r}-\mathbf{c}_i\rVert^2+\epsilon^2}}
$$

です。$\epsilon$は`sim.softening`です。softeningは点電荷の数値特異性を緩和します。三角形上の
面電荷を解像する場合はtriangle P0を使い、point kernelを使う場合はsoftening依存性を確認します。

`softening=0`で評価点がsource重心と一致した寄与は、ゼロ除算を避けるため通常の一点評価ではスキップされます。
したがって、この設定で重心上の物理的な自己場が定義されたことにはなりません。

<a id="triangle-p0"></a>

## 面電荷を三角形上で積分する

`field.element_kernel="triangle_p0"`では、要素面積$A_i$上に一定面密度
$\sigma_i=q_i/A_i$を置きます。

$$
\phi_i(\mathbf{r})=
k_c\int_{T_i}\frac{\sigma_i}{\lVert\mathbf{r}-\mathbf{r}'\rVert}\,dS',
\qquad
\mathbf{E}_i(\mathbf{r})=
k_c\int_{T_i}\sigma_i
\frac{\mathbf{r}-\mathbf{r}'}{\lVert\mathbf{r}-\mathbf{r}'\rVert^3}\,dS'.
$$

BEACHは三角形の辺に沿う対数項と立体角を使う解析kernelでこの積分を評価します。
遠方と近傍の寄与を同じpanel kernelで評価し、重心点電荷へのfallbackは行いません。

面上の電位は連続で、自己電位にはprincipal valueを使います。電場の法線成分には面電荷によるjumpがあり、
真空側の極限は

$$
\mathbf{E}^{\pm}=\mathbf{E}^{\mathrm{PV}}
\pm\frac{\sigma_i}{2\epsilon_0}\mathbf{n}_i
$$

となります。この符号を決めるため、各要素に`normal_plus`または`normal_minus`のvacuum sideが必要です。
OBJ入力では`[mesh].surface_side`、templateでは各`[[mesh.templates]].surface_side`で指定します。
`outward_closed`は、向きが整合した閉じたtwo-manifoldにだけ使用できます。

triangle P0の現行要件は次のとおりです。

- `sim.softening = 0`
- 有限で面積が正の三角形
- 全要素でvacuum sideが解決済み
- Phase 1では全表面が`insulator`
- solverは`direct`、`fmm`、またはその二つを選ぶ`auto`

設定例は[`examples/panel_direct.toml`](../examples/panel_direct.toml)にあります。

<a id="要素中心の電位出力"></a>

## 要素中心電位では自己項を加える

`potential_history.csv`などの要素中心電位は、任意点の電位評価とは自己項の扱いが異なります。

triangle P0では、自要素を含む解析panel積分をそのまま使います。point kernelでは、自要素の有限な代表値を
別に加えます。$h_i=\sqrt{A_i}$として、自己係数は

$$
C_{ii}=\begin{cases}
1/\epsilon, & \epsilon>0,\\
2\sqrt{\pi}/h_i, & \epsilon=0
\end{cases}
$$

で、自己電位は$k_c C_{ii}q_i$です。後者はpoint kernelの出力用自己項であり、triangle P0の厳密な
三角形自己積分とは異なります。要素中心電位を物理的に解釈するときは、kernelの違いを保ったまま比較してください。

要素中心を全要素について評価するため、Directの電位履歴は$O(N^2)$です。履歴を細かく保存すると、
粒子位置の場評価とは別に大きなコストになることがあります。

## commit済み電荷を直接評価する

内部座標は`field_normalization`で正規化され、電場と電位はSIへ戻されます。式と選択肢は
[場の評価](FieldSolvers.html#長さの正規化)を参照してください。

Directは固定geometryのtreeや展開係数を持たず、各評価で現在の`mesh%q_elem`を直接読みます。batch中は
snapshotを固定し、batch末尾にcommitした電荷を次batchの場へ反映する更新順序は他のsolverと同じです。

## Directを基準に誤差を切り分ける

Directにはsolverによる近似誤差がありませんが、次の誤差は残ります。

- pointかtriangle P0かというsource離散化
- surface meshの解像度と形状誤差
- 粒子の時間刻みと衝突判定誤差
- point kernelのsoftening依存性
- periodic/outer modelを含む場合のmodel化誤差

Directとの差が小さければ、FMMやTreecodeの近似誤差が小さいことは確認できます。ただし、それだけでは
物理解全体の収束を保証できません。source meshと粒子時間刻みも独立に細分化してください。

## Code reference

- Directの電場・電位と自己項: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- triangle P0解析kernel: [`bem_panel_kernel.f90`](../src/physics/panel/bem_panel_kernel.f90)
- panel geometry: [`bem_panel_geometry.f90`](../src/physics/panel/bem_panel_geometry.f90)
- 設定の検証: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)
- 主な回帰テスト: [`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90)、[`test_panel_kernel.f90`](../tests/fortran/test_panel_kernel.f90)
