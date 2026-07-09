project: BEACH
summary: BEM + Accumulated Charge simulator with a Fortran runtime and Python post-processing utilities.
author: beach contributors
license: Apache-2.0
doc_license:
version: 1.4.0
src_dir: src
src_dir: app
output_dir: build/ford-docs
project_github: https://github.com/Nkzono99/BEACH
project_download: https://github.com/Nkzono99/BEACH
graph: true
coloured_edges: true
source: true
search: true
preprocess: true

# BEACH Fortran Documentation

このサブサイトは、BEACH の Fortran 実装を中心にまとめた API 参照ドキュメントです。

- `src/` と `app/` の module / submodule / program を FORD で可視化
- `!>` コメントを元に API ページを生成
- Starlight で生成する利用者向けドキュメント本体から `/fortran/` としてリンク

現行仕様の source of truth は `SPEC.md` と Fortran 実装本体です。  
このサブサイトは、実装 API とソースコメントを確認するための参照先として使う想定です。
