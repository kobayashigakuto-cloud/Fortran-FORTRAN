# Modern Fortran

旧形式のFortranプログラムを基に、Modern Fortranで再設計した解析プログラムを保存しています。

## Design Goals

- 共通処理の整理
- Moduleによる機能分離
- allocatable配列による柔軟なデータ管理
- 明示的な型宣言
- 可読性・保守性・拡張性の向上
- 重複処理の削減
- 複数の解析機能の統合

各プログラムの機能、入出力形式、実行方法については、それぞれのディレクトリ内のREADMEを参照してください。

## Migration Policy

Legacy版の処理結果と比較しながら、機能単位で段階的に移行しています。最終的にはModern Fortran版のみで解析処理を完結できる構成を目指しています。

## Legacy Version

移行元のプログラムは `../Legacy-Fortran/` に保存しています。
