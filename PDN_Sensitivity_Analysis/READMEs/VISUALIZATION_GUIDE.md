# PDN 視覺化指南

這個工具提供了多種方式來視覺化 PDN 網路結構。

## 方法 1: Python 視覺化（推薦）

### 安裝依賴

```bash
pip install matplotlib numpy
```

### 基本使用

```bash
# 顯示 5x5 網格的視覺化
python3 visualize_pdn.py --rows 5 --cols 5

# 生成不同類型的視覺化
python3 visualize_pdn.py --type grid --rows 10 --cols 10
python3 visualize_pdn.py --type heatmap --direction x --rows 10 --cols 10
python3 visualize_pdn.py --type graph --rows 5 --cols 5

# 保存為圖片（不顯示）
python3 visualize_pdn.py --rows 5 --cols 5 --output pdn_grid.png --no-show
```

### 視覺化類型

#### 1. Grid 視覺化（預設）
顯示網格結構，包含：
- 節點位置和標籤
- 水平分支（Rx）和垂直分支（Ry）
- 電阻值標註
- 顏色編碼（綠色=低電阻，橙色=中電阻，紅色=高電阻）

```bash
python3 visualize_pdn.py --type grid --rows 5 --cols 5
```

#### 2. Heatmap 視覺化
以熱圖形式顯示電阻值分佈：
- 顏色深度表示電阻大小
- 數值標註在每個網格點
- 適合觀察電阻分佈模式

```bash
# 顯示 Rx 的熱圖
python3 visualize_pdn.py --type heatmap --direction x --rows 10 --cols 10

# 顯示 Ry 的熱圖
python3 visualize_pdn.py --type heatmap --direction y --rows 10 --cols 10
```

#### 3. Graph 視覺化
以圖形形式顯示網路拓撲：
- 節點和連接的清晰表示
- 適合理解網路結構

```bash
python3 visualize_pdn.py --type graph --rows 5 --cols 5
```

## 方法 2: Graphviz 視覺化（DOT 格式）

### 安裝 Graphviz

```bash
# Ubuntu/Debian
sudo apt-get install graphviz

# macOS
brew install graphviz
```

### 生成 DOT 檔案

使用 C++ 程式生成 DOT 檔案（需要整合到主程式中）：

```cpp
// 在程式中呼叫
generateDOTFile(network, "pdn_network.dot");
```

### 轉換為圖片

```bash
# PNG 格式
dot -Tpng pdn_network.dot -o pdn_network.png

# SVG 格式（可縮放）
dot -Tsvg pdn_network.dot -o pdn_network.svg

# PDF 格式
dot -Tpdf pdn_network.dot -o pdn_network.pdf
```

## 方法 3: ASCII 文字視覺化

在程式中直接輸出 ASCII 文字表示：

```cpp
printNetworkASCII(network);
```

輸出範例：
```
=== PDN Network ASCII Visualization ===

(0,0) --Rx=0.100--> (0,1) --Rx=0.100--> (0,2)
  | Ry=0.100          | Ry=0.100          |
(1,0) --Rx=0.100--> (1,1) --Rx=0.100--> (1,2)
  | Ry=0.100          | Ry=0.100          |
(2,0) --Rx=0.100--> (2,1) --Rx=0.100--> (2,2)
```

## 整合到敏感度分析

### 在分析前後視覺化

```python
# 在 sensitivity_analysis.cpp 的 main 函數中加入
# 或者在 Python 腳本中：

from visualize_pdn import PDNVisualizer, create_sample_pdn

# 創建視覺化器
viz = PDNVisualizer(10, 10, 1)
pixel_models = create_sample_pdn(10, 10)

# 分析前的狀態
viz.visualize_grid_2d(pixel_models, save_file="before_analysis.png", show=False)

# 執行敏感度分析...

# 分析後的狀態（修改後的 R 值）
# 更新 pixel_models 並重新視覺化
viz.visualize_grid_2d(pixel_models, save_file="after_analysis.png", show=False)
```

## 自訂視覺化

### 修改顏色方案

在 `visualize_pdn.py` 中修改 `_get_resistance_color` 函數：

```python
def _get_resistance_color(self, resistance, min_r, max_r):
    # 自訂顏色邏輯
    if resistance < 0.05:
        return 'green'
    elif resistance < 0.15:
        return 'yellow'
    else:
        return 'red'
```

### 添加電壓/電流資訊

可以擴展視覺化以顯示：
- 節點電壓（用顏色或數值）
- 分支電流（用線條粗細）
- IR drop（用顏色深淺）

## 範例腳本

創建一個完整的視覺化腳本：

```python
#!/usr/bin/env python3
import sys
sys.path.append('.')

from visualize_pdn import PDNVisualizer, create_sample_pdn

def main():
    rows, cols = 8, 8
    
    # 創建視覺化器
    viz = PDNVisualizer(rows, cols)
    
    # 創建範例 PDN
    pixel_models = create_sample_pdn(rows, cols)
    
    # 生成所有類型的視覺化
    print("Generating grid visualization...")
    viz.visualize_grid_2d(pixel_models, save_file="pdn_grid.png", show=False)
    
    print("Generating Rx heatmap...")
    viz.visualize_resistance_heatmap(pixel_models, direction='x', 
                                    save_file="pdn_rx_heatmap.png", show=False)
    
    print("Generating Ry heatmap...")
    viz.visualize_resistance_heatmap(pixel_models, direction='y', 
                                    save_file="pdn_ry_heatmap.png", show=False)
    
    print("All visualizations generated!")

if __name__ == '__main__':
    main()
```

## 輸出檔案

視覺化工具會生成以下檔案：
- `pdn_grid.png`: 網格視覺化
- `pdn_rx_heatmap.png`: Rx 熱圖
- `pdn_ry_heatmap.png`: Ry 熱圖
- `pdn_network.dot`: Graphviz DOT 檔案
- `pdn_network.png`: 從 DOT 生成的圖片

## 提示

1. **大型網路**：對於大型網路（>20×20），建議使用熱圖視覺化，因為網格視覺化會太擁擠
2. **對比分析**：生成分析前後的視覺化，可以清楚看到 R 值變化的影響
3. **批次處理**：可以寫腳本批次生成多個視覺化，用於報告或演示

