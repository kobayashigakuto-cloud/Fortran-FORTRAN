# Legacy Bulk Water

水バルク系の解析に使用しているF77形式のプログラムを保存しています。

## ディレクトリ構成

- `active/`
  - 現在も用途別に使用しているプログラム

- `integrated/`
  - 複数の解析処理を統合・更新したプログラム
  - Modern Fortran版の基礎となったコード

## Modern Fortran版

`integrated`を基に、共通処理の整理、動的配列の導入、モジュール化などを行った版は、`Modern-BulkWater`に保存しています。
