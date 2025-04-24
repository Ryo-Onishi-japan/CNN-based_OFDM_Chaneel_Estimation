事前設定
- test.pyを開き、右下のpythn環境を3.11にする
- addpath("module","model","train_data")

| パラメータ名          | 最適化範囲               | デフォルト値             | 型       | 変換     |
|---------------------|-----------------------|---------------------|--------|--------|
| InitialLearnRate    | [1e-5, 1e-1]  | 1e-3 | real   | log    |
| L2Regularization    | [1e-6, 1e-2]  | 1e-4 | real   | log    |
| MiniBatchPower      | [3, 11]               | 7（2^7=128）          | integer| none   |
| GradientDecayFactor | [0.60, 0.99]      | 0.9                 | real   | none   |