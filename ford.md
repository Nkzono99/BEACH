project: BEACH
summary: BEM + Accumulated Charge simulator with a Fortran runtime and Python post-processing utilities.
author: beach contributors
license: Apache-2.0
doc_license:
version: 0.7.3
src_dir: src
src_dir: app
page_dir: docs
media_dir: docs/media
output_dir: build/ford-docs
project_github: https://github.com/Nkzono99/BEACH
project_download: https://github.com/Nkzono99/BEACH
graph: true
coloured_edges: true
source: true
search: true
preprocess: true
md_base_dir: .

# BEACH Fortran Documentation

このサイトは、BEACH の Fortran 実装を中心にまとめた参照ドキュメントです。

- `src/` と `app/` の module / submodule / program を FORD で可視化
- `!>` コメントを元に API ページを生成
- `docs/` 配下の運用ドキュメントをそのまま Pages へ掲載
- `Fortran 依存関係マップ` で `use` 依存と `submodule` 親子関係を俯瞰

現行仕様の source of truth は `SPEC.md` と Fortran 実装本体です。  
このドキュメントサイトは、その理解を助けるためのナビゲーションとして使う想定です。
