# PDN Sensitivity Analysis Tool

這個工具用於分析 pixel model 中 R 值（Rx, Ry, Rz）對 PDN 網路特性的影響，幫助了解找到「好的 R 值」的重要性。

## 功能

- **Pixel Model**: 表示 PDN 網格中每個像素的電阻模型（Rx, Ry, Rz）
- **PDN Network**: 簡化的 PDN 網路結構，支援網格建模
- **Sensitivity Analysis**: 掃描 R 值並分析其對網路電阻的影響

## 詳細說明文件

- **[代碼詳細解釋](CODE_EXPLANATION.md)**: 完整解釋每個類別、方法和代碼邏輯
- **[研究方法說明](RESEARCH_METHODOLOGY.md)**: 詳細說明研究問題、方法和步驟
- **[敏感度分析詳細說明](SENSITIVITY_ANALYSIS_EXPLANATION.md)**: 完整解釋敏感度分析的概念、方法和應用
- **[實例演練](EXAMPLE_WALKTHROUGH.md)**: 通過具體例子逐步演示分析方法

## 編譯

```bash
mkdir build
cd build
cmake ..
make
```

## 使用方式

### 基本使用

```bash
./pdn_sensitivity
```

這會執行預設的敏感度分析：
- 掃描 Rx 從 0.5x 到 2.0x
- 掃描 Ry 從 0.5x 到 2.0x  
- 掃描 Rz 從 0.5x 到 2.0x
- 同時掃描所有 R 值

結果會輸出到 CSV 檔案：
- `sensitivity_rx.csv`: Rx 的敏感度分析結果
- `sensitivity_ry.csv`: Ry 的敏感度分析結果
- `sensitivity_rz.csv`: Rz 的敏感度分析結果
- `sensitivity_all.csv`: 所有 R 值同時變化的結果

### 自訂分析

修改 `src/sensitivity_analysis.cpp` 中的參數：

```cpp
// 基礎 R 值
float baseRx = 0.1f;  // 0.1 Ω
float baseRy = 0.1f;  // 0.1 Ω
float baseRz = 0.5f;  // 0.5 Ω

// 掃描範圍和步數
float startRatio = 0.5f;  // 起始倍數
float endRatio = 2.0f;    // 結束倍數
int steps = 20;           // 步數

// 網格大小
int rows = 10, cols = 10, layers = 1;
```

## 輸出格式

CSV 檔案包含以下欄位：
- `Rx, Ry, Rz`: 電阻值（Ω）
- `TotalResistance`: 網路總電阻
- `AvgBranchResistance`: 平均分支電阻
- `MaxBranchResistance`: 最大分支電阻
- `MinBranchResistance`: 最小分支電阻
- `PathResistance`: 特定路徑的電阻

## 分析結果解讀

1. **敏感度 (Sensitivity)**: `ΔTotalR / ΔR` 表示 R 值變化對總電阻的影響程度
   - 值越大，表示該 R 值對網路影響越大
   - 可以幫助識別哪些 R 值最關鍵

2. **電阻範圍**: 顯示 R 值變化時，網路電阻的變化範圍
   - 範圍越大，表示該 R 值對結果影響越大

3. **分支電阻統計**: 了解 R 值變化對網路中個別分支的影響

## 結果分析和繪圖

### CSV 結果繪圖和分析

**安裝依賴**（選擇一種方式）：

方式 1：使用完整版本（需要 pandas）
```bash
pip install pandas matplotlib numpy
python3 plot_sensitivity_results.py
```

方式 2：使用簡化版本（只需要 matplotlib，不需要 pandas）
```bash
pip install matplotlib numpy
python3 plot_sensitivity_results_simple.py
```

執行敏感度分析後，使用以下腳本生成圖表和詳細分析：

```bash
# 先執行分析生成 CSV
./pdn_sensitivity

# 生成圖表和報告（使用簡化版本，不需要 pandas）
python3 plot_sensitivity_results_simple.py

# 查看報告
cat sensitivity_report.txt

# 查看圖表
ls plots/*.png
```

腳本會自動：
- 讀取所有 CSV 檔案
- 生成多種圖表（電阻 vs 總電阻、分支統計、敏感度比較）
- 生成詳細的分析報告（包含現象解釋）
- 提供設計建議

詳細說明請參考 [PLOT_ANALYSIS_GUIDE.md](PLOT_ANALYSIS_GUIDE.md)

## 視覺化

工具提供了多種視覺化方式來查看 PDN 網路結構：

### Python 視覺化（推薦）

```bash
# 基本使用
python3 visualize_pdn.py --rows 5 --cols 5

# 不同類型的視覺化
python3 visualize_pdn.py --type grid --rows 10 --cols 10
python3 visualize_pdn.py --type heatmap --direction x --rows 10 --cols 10
python3 visualize_pdn.py --type graph --rows 5 --cols 5

# 保存為圖片
python3 visualize_pdn.py --rows 5 --cols 5 --output pdn_grid.png --no-show
```

詳細說明請參考 [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)

## 擴展功能

可以進一步擴展以分析：
- IR drop 的影響（需要加入電壓源和電流負載）
- 特定路徑的電阻變化
- 多層網路的影響
- 不同網格大小的影響

## 注意事項

- 目前是簡化版本，專注於 R 值敏感度分析
- 實際的 IR drop 計算需要完整的 SPICE 模擬
- 可以根據需要擴展功能

