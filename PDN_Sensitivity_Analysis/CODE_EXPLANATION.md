# 代碼和研究方法詳細解釋

## 一、項目整體目標

### 研究問題
**「Pixel Model 中的 R 值（Rx, Ry, Rz）對 PDN 分析結果的影響到底有多大？」**

這個問題的意義：
- 如果 R 值不準確，會導致 IR drop 分析結果不準確
- 我們想知道：R 值的誤差會造成多大的結果誤差？
- 這可以幫助判斷「找到好的 R 值」的重要性

### 研究方法
使用**敏感度分析（Sensitivity Analysis）**：
1. 系統性地改變 R 值（例如從 0.5 倍到 2.0 倍）
2. 觀察網路總電阻的變化
3. 計算敏感度：`敏感度 = Δ總電阻 / ΔR值`
4. 如果敏感度高 → R 值很重要；敏感度低 → R 值影響較小

---

## 二、代碼結構和文件說明

### 文件組織

```
PDN_Sensitivity_Analysis/
├── include/
│   ├── PixelModel.h          # Pixel Model 類別定義
│   └── PDNNetwork.h          # PDN 網路類別定義
├── src/
│   ├── PixelModel.cpp        # Pixel Model 實作
│   ├── PDNNetwork.cpp        # PDN 網路實作
│   └── sensitivity_analysis.cpp  # 敏感度分析主程式
└── ...
```

---

## 三、核心類別詳細解釋

### 1. PixelModel（像素模型）

#### 作用
代表 PDN 網格中**一個像素點**的電阻模型。

#### 為什麼需要這個？
在 PDN 建模中，我們把整個晶片劃分成網格，每個網格點（pixel）需要一個模型來描述：
- 水平方向的電阻（Rx）
- 垂直方向的電阻（Ry）
- 層間方向的電阻（Rz）

#### 代碼解析

```cpp
class PixelModel {
private:
    std::string _name;    // 模型名稱
    float _rx;           // 水平電阻（X 方向）
    float _ry;           // 垂直電阻（Y 方向）
    float _rz;           // 層間電阻（Z 方向）
    float _voltage;      // 電壓值
```

**關鍵方法**：

1. **`setResistances(rx, ry, rz)`**
   ```cpp
   void setResistances(float rx, float ry, float rz) {
       setRx(rx);  // 會檢查 rx >= 0
       setRy(ry);  // 會檢查 ry >= 0
       setRz(rz);  // 會檢查 rz >= 0
   }
   ```
   - **作用**：設定三個方向的電阻值
   - **為什麼要檢查**：電阻不能是負數（物理上不合理）

2. **`scaleResistances(ratio)`**
   ```cpp
   void scaleResistances(float ratio) {
       _rx *= ratio;
       _ry *= ratio;
       _rz *= ratio;
   }
   ```
   - **作用**：按比例縮放所有電阻值
   - **用途**：在敏感度分析中，我們需要改變 R 值，這個方法很方便

#### 實際例子
```cpp
// 創建一個 pixel model
PixelModel model("pixel_0_0", 0.1, 0.1, 0.5, 1.0);
// 表示：Rx=0.1Ω, Ry=0.1Ω, Rz=0.5Ω, 電壓=1.0V

// 在敏感度分析中，我們可能這樣改變它：
model.scaleResistances(1.5);  // 所有 R 值變成 1.5 倍
// 現在：Rx=0.15Ω, Ry=0.15Ω, Rz=0.75Ω
```

---

### 2. PDNNetwork（PDN 網路）

#### 作用
建立和管理整個 PDN 網路的結構。

#### 為什麼需要這個？
PDN 不是單個 pixel，而是由**很多 pixel 連接起來**的網路：
- 每個 pixel 是一個節點（Node）
- Pixel 之間的連接是分支（Branch）
- 我們需要計算整個網路的總電阻

#### 網路結構示意

```
假設 3×3 網格：

(0,0) ──Rx── (0,1) ──Rx── (0,2)
  │           │           │
  Ry         Ry          Ry
  │           │           │
(1,0) ──Rx── (1,1) ──Rx── (1,2)
  │           │           │
  Ry         Ry          Ry
  │           │           │
(2,0) ──Rx── (2,1) ──Rx── (2,2)
```

#### 代碼解析

**1. 資料結構**

```cpp
class PDNNetwork {
private:
    int _rows, _cols, _layers;  // 網格大小
    
    // 儲存所有 pixel model
    std::vector<std::vector<std::vector<PixelModel*>>> _pixelModels;
    // [layer][row][col] -> PixelModel*
    
    // 儲存網路節點
    std::vector<Node> _nodes;
    // Node 包含：id, voltage, currentLoad, pgNetName
    
    // 儲存網路分支
    std::vector<Branch> _branches;
    // Branch 包含：node1Id, node2Id, resistance, direction
};
```

**2. `buildNetwork()` 方法（關鍵！）**

```cpp
void PDNNetwork::buildNetwork() {
    // 1. 清空舊的網路
    _nodes.clear();
    _branches.clear();
    _nodeMap.clear();
    
    // 2. 為每個網格點創建節點
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                getNodeId(l, r, c, "VDD");  // 創建 VDD 節點
                getNodeId(l, r, c, "GND");  // 創建 GND 節點
            }
        }
    }
    
    // 3. 根據 pixel model 創建分支
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                PixelModel* model = _pixelModels[l][r][c];
                if (!model) continue;
                
                // 水平分支（X 方向）
                if (c < _cols - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r, c + 1, "VDD");
                    addBranch(node1, node2, model->getRx(), "x");
                }
                
                // 垂直分支（Y 方向）
                if (r < _rows - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r + 1, c, "VDD");
                    addBranch(node1, node2, model->getRy(), "y");
                }
                
                // 層間分支（Z 方向）
                if (l < _layers - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l + 1, r, c, "VDD");
                    addBranch(node1, node2, model->getRz(), "z");
                }
            }
        }
    }
}
```

**這段代碼做了什麼？**

1. **創建節點**：為每個網格點創建節點（代表電壓點）
2. **創建分支**：根據 pixel model 的 R 值創建連接
   - 如果 pixel (r, c) 的 Rx = 0.1Ω
   - 那麼節點 (r, c) 和 (r, c+1) 之間就有一個 0.1Ω 的分支
3. **建立網路**：把所有節點和分支組合起來，形成完整的 PDN 網路

**為什麼要這樣做？**
- 我們需要計算**整個網路**的總電阻
- 總電阻 = 所有分支電阻的總和（簡化模型）
- 當 R 值改變時，分支電阻改變，總電阻也改變

---

### 3. SensitivityAnalyzer（敏感度分析器）

#### 作用
執行敏感度分析的核心邏輯。

#### 研究方法的程式化

**方法 1：`sweepResistance()` - 掃描單一 R 值**

```cpp
std::vector<SensitivityResult> sweepResistance(
    const std::string& resistanceType,  // "Rx", "Ry", 或 "Rz"
    float baseRx, float baseRy, float baseRz,
    float startRatio, float endRatio, int steps)
```

**步驟詳解**：

```cpp
// 步驟 1: 設定掃描範圍
float stepSize = (endRatio - startRatio) / (steps - 1);
// 例如：從 0.5 到 2.0，20 步
// stepSize = (2.0 - 0.5) / 19 = 0.079

// 步驟 2: 對每個步驟
for (int i = 0; i < steps; i++) {
    // 2.1 計算新的 R 值
    float ratio = startRatio + i * stepSize;
    // i=0: ratio=0.5, i=1: ratio=0.579, ..., i=19: ratio=2.0
    
    float rx = baseRx;
    float ry = baseRy;
    float rz = baseRz;
    
    // 只改變目標 R 值
    if (resistanceType == "Rx") rx *= ratio;
    else if (resistanceType == "Ry") ry *= ratio;
    else if (resistanceType == "Rz") rz *= ratio;
    
    // 2.2 應用到所有 pixel model
    applyResistancesToNetwork(rx, ry, rz);
    // 這會更新網路中所有 pixel model 的 R 值
    
    // 2.3 重建網路
    _network->buildNetwork();
    // 因為 R 值改變了，分支電阻也改變了，需要重建
    
    // 2.4 計算指標
    SensitivityResult result;
    result.rx = rx;
    result.ry = ry;
    result.rz = rz;
    calculateMetrics(result);
    // 計算總電阻、平均分支電阻等
    
    // 2.5 記錄結果
    results.push_back(result);
}
```

**實際例子**：

假設我們要掃描 Rx，基礎值 = 0.1Ω：

| 步驟 | ratio | Rx (Ω) | 網路總電阻 (Ω) |
|------|-------|--------|----------------|
| 0    | 0.5   | 0.05   | 0.9            |
| 1    | 0.579 | 0.058  | 1.04           |
| ...  | ...   | ...    | ...            |
| 19   | 2.0   | 0.20   | 1.8            |

**方法 2：`calculateMetrics()` - 計算指標**

```cpp
void calculateMetrics(SensitivityResult& result) {
    auto branches = _network->getBranches();
    
    if (branches.empty()) return;
    
    // 計算總電阻
    double total = 0.0;
    double min = branches[0].resistance;
    double max = branches[0].resistance;
    
    for (const auto& branch : branches) {
        total += branch.resistance;  // 累加所有分支電阻
        if (branch.resistance < min) min = branch.resistance;
        if (branch.resistance > max) max = branch.resistance;
    }
    
    result.totalResistance = total;  // 總電阻
    result.avgBranchResistance = total / branches.size();  // 平均分支電阻
    result.maxBranchResistance = max;  // 最大分支電阻
    result.minBranchResistance = min;  // 最小分支電阻
}
```

**這段代碼做了什麼？**

1. **取得所有分支**：從網路中取得所有連接分支
2. **計算總電阻**：把所有分支電阻加起來
   - 這是簡化模型：實際 PDN 是並聯/串聯混合，但這裡簡化為總和
3. **計算統計值**：平均、最大、最小分支電阻

**為什麼要計算這些？**
- **總電阻**：主要指標，用來評估 R 值變化的影響
- **平均分支電阻**：了解整體電阻水平
- **最大/最小**：識別網路中的瓶頸

**方法 3：`printSensitivitySummary()` - 計算敏感度**

```cpp
void printSensitivitySummary(const std::vector<SensitivityResult>& results,
                            const std::string& resistanceType) {
    if (results.size() >= 2) {
        // 計算 R 值的變化
        double resistanceChange = results.back().rx - results.front().rx;
        if (resistanceType == "Ry") {
            resistanceChange = results.back().ry - results.front().ry;
        } else if (resistanceType == "Rz") {
            resistanceChange = results.back().rz - results.front().rz;
        }
        
        // 計算總電阻的變化
        double totalResChange = results.back().totalResistance - 
                               results.front().totalResistance;
        
        // 計算敏感度
        double sensitivity = (resistanceChange != 0) ? 
                            totalResChange / resistanceChange : 0.0;
        
        std::cout << "Sensitivity: " << sensitivity 
                  << " (ΔTotalR / Δ" << resistanceType << ")" << std::endl;
    }
}
```

**敏感度計算公式**：

```
敏感度 = (總電阻變化) / (R 值變化)
       = ΔTotalR / ΔR
```

**例子**：
- Rx 從 0.05Ω → 0.20Ω（變化 0.15Ω）
- 總電阻從 0.9Ω → 1.8Ω（變化 0.9Ω）
- 敏感度 = 0.9 / 0.15 = 6.0

**解讀**：
- 敏感度 = 6.0 表示：Rx 每增加 0.01Ω，總電阻增加 0.06Ω
- 敏感度越高，表示 R 值對結果影響越大

---

## 四、主程式流程（main 函數）

```cpp
int main() {
    // === 步驟 1: 建立測試網路 ===
    int rows = 10, cols = 10, layers = 1;
    PDNNetwork network(rows, cols, layers);
    
    // 創建 pixel models
    float baseRx = 0.1f, baseRy = 0.1f, baseRz = 0.5f;
    
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            // 為每個網格點創建一個 pixel model
            auto model = std::make_unique<PixelModel>(
                "model_" + std::to_string(r) + "_" + std::to_string(c),
                baseRx, baseRy, baseRz, 1.0f
            );
            network.setPixelModel(0, r, c, model.get());
        }
    }
    
    // 建立網路結構
    network.buildNetwork();
    
    // === 步驟 2: 執行敏感度分析 ===
    SensitivityAnalyzer analyzer(&network);
    
    // 掃描 Rx
    auto rxResults = analyzer.sweepResistance("Rx", baseRx, baseRy, baseRz, 
                                              0.5f, 2.0f, 20);
    // 從 0.5 倍到 2.0 倍，20 個步驟
    
    // 掃描 Ry
    auto ryResults = analyzer.sweepResistance("Ry", baseRx, baseRy, baseRz, 
                                              0.5f, 2.0f, 20);
    
    // 掃描 Rz
    auto rzResults = analyzer.sweepResistance("Rz", baseRx, baseRy, baseRz, 
                                              0.5f, 2.0f, 20);
    
    // === 步驟 3: 輸出結果 ===
    analyzer.printSensitivitySummary(rxResults, "Rx");
    analyzer.writeResultsToCSV(rxResults, "sensitivity_rx.csv");
    // ... 其他輸出
}
```

**完整流程圖**：

```
開始
  ↓
建立 10×10 網格
  ↓
為每個網格點創建 PixelModel（Rx=0.1, Ry=0.1, Rz=0.5）
  ↓
建立網路（創建節點和分支）
  ↓
┌─────────────────────────────────┐
│ 敏感度分析循環                    │
│                                 │
│ for each step (0.5x to 2.0x):  │
│   1. 計算新的 R 值               │
│   2. 更新所有 pixel model       │
│   3. 重建網路                    │
│   4. 計算總電阻                  │
│   5. 記錄結果                    │
└─────────────────────────────────┘
  ↓
計算敏感度（ΔTotalR / ΔR）
  ↓
輸出 CSV 檔案
  ↓
結束
```

---

## 五、研究方法的邏輯

### 為什麼這樣設計？

**1. 為什麼要掃描 R 值？**
- 我們不知道 R 值的準確值（可能有誤差）
- 我們想知道：如果 R 值不準確，結果會差多少？
- 通過系統性地改變 R 值，我們可以觀察結果的變化

**2. 為什麼計算總電阻？**
- 總電阻是 PDN 網路的重要指標
- 總電阻越大 → IR drop 越大（V = I × R）
- 如果 R 值不準確，總電阻也不準確，IR drop 分析就不準確

**3. 為什麼計算敏感度？**
- 敏感度量化了「R 值變化對結果的影響」
- 如果敏感度高 → R 值很重要，需要精確建模
- 如果敏感度低 → R 值影響較小，可以容忍誤差

### 研究假設

1. **簡化模型**：總電阻 = 所有分支電阻的總和
   - 實際 PDN 是並聯/串聯混合，但這裡簡化為總和
   - 這是為了快速分析，實際應用需要更複雜的模型

2. **均勻 R 值**：所有 pixel model 使用相同的 R 值
   - 實際 PDN 不同位置可能有不同的 R 值
   - 但為了分析「R 值的重要性」，這個假設是合理的

3. **線性關係**：假設 R 值變化與總電阻變化是線性的
   - 在簡化模型中，這個假設成立
   - 但在複雜網路中可能不成立

---

## 六、輸出結果的解讀

### CSV 檔案格式

```csv
Rx,Ry,Rz,TotalResistance,AvgBranchResistance,MaxBranchResistance,MinBranchResistance,PathResistance
0.05,0.10,0.50,0.9,0.075,0.05,0.10,1.35
0.10,0.10,0.50,1.2,0.100,0.10,0.10,1.80
0.15,0.10,0.50,1.5,0.125,0.15,0.10,2.25
...
```

### 如何解讀？

**1. 觀察趨勢**
- 當 Rx 增加時，TotalResistance 是否也增加？
- 如果是 → 正相關，R 值確實影響結果

**2. 計算敏感度**
```
敏感度 = (最後的 TotalR - 第一個 TotalR) / (最後的 R - 第一個 R)
       = (1.8 - 0.9) / (0.20 - 0.05)
       = 0.9 / 0.15
       = 6.0
```

**3. 評估重要性**
- 敏感度 = 6.0 → 中等影響
- 如果敏感度 = 100 → 高影響，R 值非常重要
- 如果敏感度 = 0.1 → 低影響，R 值影響很小

---

## 七、方法的優點和限制

### 優點

1. **簡單直接**：容易理解和實作
2. **快速分析**：不需要完整的 SPICE 模擬
3. **量化結果**：敏感度提供具體的數字
4. **系統性**：覆蓋整個參數空間

### 限制

1. **簡化模型**：總電阻 = 總和，不是實際的並聯/串聯
2. **沒有 IR drop**：只計算電阻，沒有實際的電壓降
3. **均勻假設**：所有位置使用相同的 R 值
4. **線性假設**：假設線性關係，可能不適用於複雜網路

### 改進方向

1. **整合 SPICE**：加入實際的 IR drop 計算
2. **非均勻 R 值**：支援不同位置使用不同的 R 值
3. **多變量分析**：同時改變多個參數
4. **統計分析**：加入不確定性分析

---

## 八、總結

### 這個項目做了什麼？

1. **建立 PDN 網路模型**：用 PixelModel 和 PDNNetwork 表示 PDN
2. **執行敏感度分析**：系統性地改變 R 值，觀察結果變化
3. **計算敏感度**：量化 R 值對結果的影響
4. **輸出結果**：CSV 檔案和摘要報告

### 研究目標達成？

**問題**：「R 值對結果的影響有多大？」

**答案**：通過敏感度分析，我們可以：
- 量化影響程度（敏感度值）
- 比較不同 R 值的重要性（Rx vs Ry vs Rz）
- 判斷是否需要精確的 R 值建模

**結論**：如果敏感度高 → 找到「好的 R 值」確實很重要！

