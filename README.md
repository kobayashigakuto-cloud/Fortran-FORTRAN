# Fortran-FORTRAN

研究で使用する数値解析プログラムと計算支援ツールを公開しています。

研究初期に作成したレガシーなFORTRANコードを基に、可読性・保守性・拡張性を高めたModern Fortranへの段階的な移行を進めています。

また、Linux環境における解析処理や計算ジョブの自動化を目的として、Shell Scriptによる補助ツールも開発しています。

## Repository Structure

- `Modern-Fortran/`
  - Modern Fortranで再設計した解析プログラム
  - Moduleによる機能分離
  - 複数の解析処理の統合

- `Legacy-Fortran/`
  - 研究で作成・使用してきた旧形式の解析プログラム
  - Modern Fortran版との比較・移行元として保存

## Technologies

- FORTRAN
- Modern Fortran
- Shell Script
- Linux

## Objectives

- 分子動力学シミュレーションの解析手法開発
- レガシーコードのModern Fortranへの移行
- 解析コードの可読性・保守性・拡張性向上
- 計算・解析工程の自動化
- 研究コードの標準化と再利用性向上

## Current Status

現在は、バルク水を対象とした解析プログラムを掲載しています。既存コードについても、整理と動作確認が完了したものから順次追加します。
