# ADR 0007: 外部境界の実行時ownership契約

- Status: accepted
- Date: 2026-07-23

> 公開TOMLを分割したまま維持する判断は
> [ADR 0008](0008-external-boundary-authoring-facade.md)で置き換えた。ここで定義した実行時ownership契約と
> particle loop境界は引き続き有効である。[ADR 0009](0009-remove-linear-debye-and-legacy-sheath.md)は
> `linear_debye` / `legacy_sheath` 分岐を、[ADR 0010](0010-remove-unified-linear-response.md)は
> unified 3D 分岐を削除したため、以下のモデル名とselector列挙は履歴上の記録である。

## Context

外部条件は、`sim.reservoir_potential_model`、`sim.sheath_injection_model`、
`sim.open_boundary_model`、`outer_plasma.model`、`outer_plasma.return_model`、
`coupling.particle_transfer_mode`、queue設定へ段階的に追加されてきた。各層が同じ文字列の組合せを
独自に判定すると、流入補正の二重適用、z-highのowner競合、return後の残り時間の取りこぼしが起こり得る。

一方、外部場の構築、粒子の流入写像、通常open面、interface輸送は別の責務である。これらを一つの
汎用境界クラスへ統合すると、現時点で存在しない差し替え点まで抽象化することになる。

## Decision

公開TOMLは変更せず、実行境界で既存設定を次の最小な内部契約へ解決する。

- z-high流入写像のowner
- outerが所有しない通常open面の処理
- z-high interface輸送
- delayed queueの有無

粒子loopはrun開始時に解決した契約をread-onlyで共有する。batch injectionも同じresolverをhot loop外で使う。
外部場の構築は既存のelectrostatic snapshotとouter coupler、個々のreturn/escape計算は既存のlinear、
kinetic 1D、unified 3D modelが引き続き担う。

1D transferを使う`linear_debye`と`kinetic_1d`は、z-high returnと流入写像を同じrefresh済みprofileで
所有する。`infinity_barrier`またはlegacy sheath補正との併用は二重補正になるため設定時に拒否する。
`open_boundary_model`は独立に保ち、outerが所有しないopen面へだけ適用する。

instant return後は同じlocal stepの残り時間を再追跡し、z-highへ再到達すれば同じ契約へ再dispatchする。
external eventは元のlocal stepあたり最大8回、box eventは各local continuationあたり最大8回とし、
別のfailure statusでfail closedにする。z-highの二次軌道補正後も周期・反射面との同時性が保たれる場合は
横方向作用を先に合成する。補正によってz-highと横面の先後・同時関係が元のface maskからどちら向きに
変わっても作用順序を推測しない。別のopen面も含む場合はownerが一意でないため、これらをすべて拒否する。

kinetic 1D/unified modelの無限遠電位は0へ固定し、kinetic ambient electron/ionおよび必要な
photoelectron speciesはexact-oneとする。これは曖昧だった既存入力を拒否する意図的なvalidity tighteningである。

## Rejected alternatives

- 公開`[boundary_environment]` bundleを直ちに追加する。既存TOMLを移行させる必要があり、外部場と時間発展の
  ownershipがまだ一つの公開概念として安定していない。
- 抽象基底型、procedure pointer、任意面interfaceを導入する。現在の所有面はz-highだけで、二つ目の実装要求がない。
- JSON Schemaへ全model組合せ表を複製する。組合せ意味論はFortran resolverとPython semantic validatorで検査し、
  schemaはenum、型、zero-gaugeなど局所制約に留める。

## Consequences

BEACHのbatch loop、Boris更新、mesh hit、電荷commitの順序は変わらない。外部条件の選択と粒子追跡の
orchestrationだけが分離され、hot loop内の文字列dispatchと暗黙のowner競合がなくなる。

z-highの二次処理は、候補終点が外側にあるchord検出済みeventの法線時刻補正に限定する。終点が内側へ戻る
途中越境の探索、x/y軌道やmesh hitの二次化はtrusted coreの変更になるため行わない。内部Fortran呼出しでは
旧`defer_z_high_interface` booleanを`boundary_contract`へ置き換えるが、公開TOMLとPython APIは変更しない。

公開bundle、動的outer fieldの時間契約、任意面transportは、具体的な二つ目の利用例が出るまで延期する。

## Validation

1. resolverの全canonical組合せ、未知値、owner競合をunit testする。
2. linear/kinetic流入がrefresh済みouter stateを使い、legacy cutoffを弱めないことを確認する。
3. 同一step内の複数return、global crossing fraction、external event上限をdriver testで固定する。
4. z-highとperiodic/reflect cornerの合成、補正後に同時性または横面との順序が崩れるevent、
   別open cornerのfail-closedをstepper testで固定する。
5. kinetic ambient/photoelectronのexact-one条件とzero-gaugeをFortran/Python/schemaで確認する。
