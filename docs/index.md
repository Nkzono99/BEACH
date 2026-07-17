title: BEACH ドキュメント

Lang: [日本語](index.md) | [English](index.en.md)

# BEACH

BEACH は、絶縁体表面への電荷蓄積と、その電荷が作る電場中の荷電粒子軌道を計算する
シミュレータです。

粒子照射による表面帯電、電位変化、粒子の吸収・脱出を、三角形境界要素メッシュ上で評価できます。

**[10 分チュートリアルを始める](Tutorial.html)** · [インストール](Installation.html)

## 最短で動かす

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

設定を作成し、検査してからシミュレーションを実行します。結果は `outputs/latest` に保存されます。

## 入力から出力まで

BEACH は、同じ電場を共有する粒子の追跡と、末尾の表面電荷更新をひとまとまりとして繰り返します。
この計算単位をバッチ（batch）と呼びます。

```text
粒子条件・表面メッシュ・境界条件
                 ↓
┌──────────── バッチ n ────────────┐
│ 表面電荷から電場を計算             │
│             ↓                     │
│ 粒子追跡 → 衝突・吸収 → 表面電荷更新 │
└───────────────────┬───────────────┘
                    ├──→ 次のバッチ n + 1 の電場へ（繰り返す）
                    └──→ 表面電荷・表面電位・粒子統計・バッチ履歴
```

各バッチでは、電場を計算して粒子を追跡し、表面に吸収された粒子の電荷を三角形要素へ
蓄積します。更新した表面電荷は、次のバッチの電場計算に反映されます。

<div align="center">
  <img src="images/potential_history.gif" alt="電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## 目的から選ぶ

- **初めて実行する：** [インストール](Installation.html) → [10 分チュートリアル](Tutorial.html)
- **研究ケースを作る：** [シミュレーションケースを設計する](ConfigurationRecipes.html) → [`beach.toml`を作成・検証する](Configuration.html) → [実行する](Execution.html) → [結果を検証する](ValidationGuide.html)
- **外部シースを構成する：** [境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)。標準は`kinetic_1d`、rough surface線形screeningは高度な`unified_linear_response`です。
- **数値モデルを理解する：** [計算モデルの全体像](Algorithms.html) → [表面電荷更新](SurfaceModels.html) → [場の評価](FieldSolvers.html)
- **開発する：** [開発・運用ワークフロー](Workflow.html) → [物理リリースの検証](PhysicsReleaseVerification.html) → [Fortran 依存関係マップ](FortranDependencyMap.html)
