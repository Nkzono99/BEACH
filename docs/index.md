title: BEACH ドキュメント

Lang: [日本語](index.md) | [English](index.en.md)

# BEACH

BEACHは、絶縁体表面への電荷蓄積と、その電荷が作る電場中の荷電粒子軌道を計算する
シミュレータです。

粒子照射による表面帯電、電位変化、粒子の吸収・脱出を、三角形境界要素メッシュ上で評価できます。

**[10分チュートリアルを始める](Tutorial.html)** · [インストール](Installation.html)

## 最短で動かす

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

設定を作成し、検査してからシミュレーションを実行します。結果は`outputs/latest`に保存されます。

## 入力から出力まで

```text
粒子条件・表面メッシュ・境界条件
                ↓
              BEACH
                ↓
表面電荷・表面電位・粒子統計・batch履歴
```

BEACHは各batchで電場を計算して粒子を追跡し、表面に吸収された粒子の電荷を三角形要素へ
蓄積します。蓄積した電荷は次のbatchの電場計算に反映されます。

<div align="center">
  <img src="images/potential_history.gif" alt="電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

詳しい使い方と数値モデルは、左側のナビゲーションから参照してください。
