# CSV 結果繪圖與分析指南

## 概述

`plot_sensitivity_results.py` 腳本會讀取敏感度分析生成的 CSV 檔案，自動生成圖表並提供詳細的現象解釋。

## 快速開始

### 1. 先執行敏感度分析

```bash
# 執行分析，生成 CSV 檔案
./pdn_sensitivity
```

這會生成：
- `sensitivity_rx.csv`
- `sensitivity_ry.csv`
- `sensitivity_rz.csv`
- `sensitivity_all.csv`

### 2. 執行繪圖腳本

```bash
# 基本使用（生成所有圖表和報告）
python3 plot_sensitivity_results.py

# 指定輸出目錄
python3 plot_sensitivity_results.py --output-dir my_plots

# 只生成報告，不生成圖表
python3 plot_sensitivity_results.py --no-plots
```

## 生成的圖表

### 1. 電阻值 vs 總電阻圖（plot_*_vs_total.png）

**檔案**：
- `plot_rx_vs_total.png`
- `plot_ry_vs_total.png`
- `plot_rz_vs_total.png`

**內容**：
- **上圖**：R 值 vs 總電阻的關係
  - 藍色實線：實際數據
  - 紅色虛線：線性擬合（顯示趨勢）
  - 文字框：敏感度值
- **下圖**：百分比變化關係
  - 顯示 R 值變化百分比 vs 總電阻變化百分比
  - 紅色虛線：1:1 參考線

**觀察重點**：
- **線性關係**：如果點大致在直線上 → 線性關係
- **斜率**：斜率越大 → 敏感度越高
- **敏感度值**：顯示在圖上，值越大影響越大

### 2. 分支電阻統計圖（plot_*_branch_stats.png）

**檔案**：
- `plot_rx_branch_stats.png`
- `plot_ry_branch_stats.png`
- `plot_rz_branch_stats.png`

**內容**：四個子圖顯示分支電阻的統計
- **平均分支電阻**：整體電阻水平
- **最小分支電阻**：網路中的最低電阻
- **最大分支電阻**：網路中的最高電阻（瓶頸）
- **電阻範圍**：最大 - 最小（變異程度）

**觀察重點**：
- **範圍大小**：範圍大 → 網路電阻分佈不均勻
- **最大電阻**：識別網路瓶頸
- **趨勢**：觀察這些統計值如何隨 R 值變化

### 3. 敏感度比較圖（plot_sensitivity_comparison.png）

**內容**：
- 柱狀圖比較 Rx、Ry、Rz 的敏感度
- 每個柱子上標註敏感度值
- 文字框顯示最高和最低敏感度

**觀察重點**：
- **哪個 R 值最重要**：敏感度最高的柱
- **相對重要性**：比較三個 R 值的敏感度差異
- **設計建議**：高敏感度的 R 值需要更精確的建模

## 生成的報告

### 報告檔案：`sensitivity_report.txt`

報告包含每個 R 值的詳細分析：

#### 1. 參數範圍
```
1. Parameter Range:
   Rx: 0.050000 Ω → 0.200000 Ω
   Change: 0.150000 Ω (+200.00%)
```
- 顯示 R 值的掃描範圍
- 變化量和變化百分比

#### 2. 總電阻響應
```
2. Total Resistance Response:
   TotalR: 0.900000 Ω → 1.800000 Ω
   Change: 0.900000 Ω (+100.00%)
```
- 顯示總電阻如何響應 R 值的變化
- 變化量和變化百分比

#### 3. 敏感度計算
```
3. Sensitivity:
   Sensitivity = ΔTotalR / ΔRx
              = 0.900000 / 0.150000
              = 6.000000
```
- 詳細的敏感度計算過程
- 最終敏感度值

#### 4. 解讀（Interpretation）

根據敏感度值自動分類：

**高敏感度（> 10）**：
```
⚠️  HIGH SENSITIVITY (12.50)
→ This R value has STRONG impact on total resistance
→ Accurate R value modeling is CRITICAL
→ Small errors in R will cause large errors in results
```

**中等敏感度（5-10）**：
```
⚡ MODERATE SENSITIVITY (6.00)
→ This R value has MODERATE impact on total resistance
→ Accurate R value modeling is IMPORTANT
→ Some tolerance for R value errors is acceptable
```

**低敏感度（< 5）**：
```
✓ LOW SENSITIVITY (2.00)
→ This R value has WEAK impact on total resistance
→ Simplified R value modeling may be acceptable
→ Larger tolerance for R value errors is acceptable
```

#### 5. 關係分析（Relationship Analysis）
```
5. Relationship Analysis:
   Correlation: 0.9998
   → Strong linear relationship
```
- 計算 R 值與總電阻的相關係數
- 判斷關係強度（強/中/弱）

#### 6. 影響評估（Impact Assessment）
```
6. Impact Assessment:
   Impact Ratio: 0.5000
   (TotalR changes 0.50x the rate of Rx)
   → TotalR is LESS sensitive than the R value itself
```
- 計算影響比率
- 判斷總電阻的變化速度相對於 R 值變化

#### 7. 分支電阻統計
```
7. Branch Resistance Statistics:
   Average: 0.100000 Ω
   Min: 0.050000 Ω
   Max: 0.200000 Ω
   Range: 0.150000 Ω
```
- 分支電阻的統計資訊
- 幫助了解網路電阻分佈

## 現象解釋

### 現象 1：線性關係

**觀察**：R 值 vs 總電阻的圖呈現直線

**解釋**：
- 在簡化模型中（總電阻 = 總和），關係是線性的
- 總電阻 = Σ(所有分支電阻)
- 如果所有分支電阻都按比例變化，總電阻也按比例變化

**意義**：
- 預測容易：知道 R 值變化，可以直接計算總電阻變化
- 敏感度是常數：不隨 R 值大小改變

### 現象 2：高敏感度

**觀察**：敏感度 > 10

**解釋**：
- R 值的小變化會導致總電阻的大變化
- 例如：敏感度 = 20 表示 R 值增加 1%，總電阻增加 20%

**意義**：
- **R 值建模很重要**：需要精確的 R 值
- **誤差放大**：R 值的誤差會被放大到總電阻
- **設計建議**：投入更多資源做精確建模

### 現象 3：低敏感度

**觀察**：敏感度 < 5

**解釋**：
- R 值的變化對總電阻影響較小
- 例如：敏感度 = 2 表示 R 值增加 1%，總電阻只增加 2%

**意義**：
- **R 值建模可以簡化**：不需要極度精確
- **誤差容忍**：可以容忍較大的 R 值誤差
- **設計建議**：可以使用簡化模型

### 現象 4：敏感度差異

**觀察**：Rx、Ry、Rz 的敏感度不同

**解釋**：
- 不同方向的 R 值對總電阻的影響不同
- 可能原因：
  - 網路結構不對稱（例如：更多水平分支）
  - 分支數量不同（例如：水平分支比垂直分支多）

**意義**：
- **優先級排序**：高敏感度的 R 值優先處理
- **資源分配**：把資源集中在最重要的 R 值上

### 現象 5：百分比變化非 1:1

**觀察**：總電阻變化百分比 ≠ R 值變化百分比

**解釋**：
- 如果敏感度 ≠ 1，變化比例不同
- 例如：R 值增加 100%，總電阻可能只增加 50%（敏感度 = 0.5）

**意義**：
- **非比例響應**：總電阻的變化速度與 R 值不同
- **影響評估**：需要考慮敏感度來評估影響

### 現象 6：分支電阻範圍

**觀察**：最大分支電阻 >> 最小分支電阻

**解釋**：
- 網路中存在電阻差異很大的分支
- 高電阻分支成為瓶頸

**意義**：
- **識別瓶頸**：找出限制網路性能的分支
- **優化目標**：優先降低高電阻分支的電阻

## 實際應用建議

### 根據敏感度值做決策

#### 高敏感度（> 10）
```
行動：
1. 使用精確的 R 值提取方法
2. 驗證 R 值的準確性
3. 考慮不確定性分析
4. 在設計中預留安全邊際

資源分配：高優先級
```

#### 中等敏感度（5-10）
```
行動：
1. 使用標準的 R 值建模方法
2. 定期驗證 R 值
3. 監控 R 值的變化

資源分配：中等優先級
```

#### 低敏感度（< 5）
```
行動：
1. 可以使用簡化的 R 值模型
2. 較大的誤差容忍度
3. 快速建模方法可接受

資源分配：低優先級
```

### 比較分析

**如果 Rx 敏感度 > Ry 敏感度 > Rz 敏感度**：
- 優先優化 Rx 的建模
- Ry 次之
- Rz 可以簡化

**如果所有敏感度都很高**：
- 所有 R 值都很重要
- 需要全面的精確建模
- 找到「好的 R 值」確實很重要

**如果所有敏感度都很低**：
- R 值對結果影響較小
- 可以使用簡化模型
- 找到「好的 R 值」不太重要

## 範例輸出

執行腳本後，你會得到：

```
Found 4 CSV file(s)
============================================================
Generating plots...
Loaded sensitivity_rx.csv: 20 data points
Saved: plots/plot_rx_vs_total.png
Saved: plots/plot_rx_branch_stats.png
...
============================================================
Analysis Complete!
============================================================
Report saved to: sensitivity_report.txt
Plots saved to: plots/
```

查看報告：
```bash
cat sensitivity_report.txt
```

查看圖表：
```bash
ls plots/*.png
```

## 進階使用

### 自訂分析

修改腳本以加入自訂分析：

```python
# 在 plot_sensitivity_results.py 中加入
def custom_analysis(df):
    # 你的自訂分析邏輯
    pass
```

### 批次處理

處理多個分析結果：

```bash
for dir in analysis_*/; do
    cd "$dir"
    python3 ../plot_sensitivity_results.py --output-dir "../plots_$(basename $dir)"
    cd ..
done
```

## 總結

這個腳本提供了：
1. **自動化分析**：讀取 CSV，自動生成圖表和報告
2. **詳細解釋**：每個現象都有詳細的解釋
3. **實用建議**：根據敏感度值提供設計建議
4. **視覺化**：多種圖表幫助理解結果

通過這個工具，你可以：
- 快速了解 R 值對結果的影響
- 量化「找到好的 R 值」的重要性
- 做出基於數據的設計決策


