# PDN Modeling Extraction Summary

## 概述

這個專案從 NTU-3DIC 中提取了 PDN (Power Delivery Network) 建模的核心部分，專門用於分析 **pixel model 中 R 值對結果的影響**（敏感度分析）。

## 提取的組件

### 1. PixelModel (像素模型)
- **位置**: `include/PixelModel.h`, `src/PixelModel.cpp`
- **功能**: 表示 PDN 網格中每個像素的電阻模型
- **關鍵參數**:
  - `Rx`: 水平方向電阻
  - `Ry`: 垂直方向電阻  
  - `Rz`: 層間垂直電阻
  - `Voltage`: 電壓值

### 2. PDNNetwork (PDN 網路)
- **位置**: `include/PDNNetwork.h`, `src/PDNNetwork.cpp`
- **功能**: 簡化的 PDN 網路結構，支援網格建模
- **特性**:
  - 網格基礎的節點和分支結構
  - 支援多層網路
  - 可以計算路徑電阻

### 3. SensitivityAnalyzer (敏感度分析工具)
- **位置**: `src/sensitivity_analysis.cpp`
- **功能**: 掃描 R 值並分析其對網路特性的影響
- **分析類型**:
  - 單一 R 值掃描 (Rx, Ry, 或 Rz)
  - 所有 R 值同時掃描
  - 輸出 CSV 格式結果

## 與原始 NTU-3DIC 的差異

### 簡化的部分
1. **移除複雜依賴**: 不依賴完整的 3dblox 結構、TdbxMgr 等
2. **簡化網路模型**: 專注於電阻網路，不包含完整的 SPICE 模擬
3. **獨立執行**: 可以獨立編譯和執行，不需要完整的 NTU-3DIC 環境

### 保留的核心功能
1. **Pixel Model 概念**: 保留了 Rx, Ry, Rz 的核心概念
2. **網格建模**: 保留了網格基礎的 PDN 建模方法
3. **電阻計算**: 保留了電阻網路的計算邏輯

## 使用目的

這個工具的主要目的是：
- **測試 R 值敏感度**: 了解 Rx, Ry, Rz 對網路總電阻的影響程度
- **評估 R 值重要性**: 幫助判斷找到「好的 R 值」的重要性
- **快速分析**: 不需要完整的 SPICE 模擬，可以快速掃描參數空間

## 編譯和使用

### 編譯
```bash
cd /home/piohuang/PDN_Sensitivity_Analysis
mkdir build && cd build
cmake ..
make
```

或使用 Makefile:
```bash
make
```

### 執行
```bash
./pdn_sensitivity
```

### 輸出
- `sensitivity_rx.csv`: Rx 的敏感度分析結果
- `sensitivity_ry.csv`: Ry 的敏感度分析結果
- `sensitivity_rz.csv`: Rz 的敏感度分析結果
- `sensitivity_all.csv`: 所有 R 值同時變化的結果

## 擴展建議

如果需要更完整的分析，可以考慮：
1. **整合 SPICE 模擬**: 加入實際的 IR drop 計算
2. **多層網路**: 支援更複雜的多層結構
3. **非均勻 R 值**: 支援不同位置使用不同的 R 值
4. **視覺化**: 加入結果視覺化功能

## 檔案結構

```
PDN_Sensitivity_Analysis/
├── include/
│   ├── PixelModel.h          # Pixel model 定義
│   └── PDNNetwork.h          # PDN 網路定義
├── src/
│   ├── PixelModel.cpp        # Pixel model 實作
│   ├── PDNNetwork.cpp        # PDN 網路實作
│   └── sensitivity_analysis.cpp  # 敏感度分析主程式
├── CMakeLists.txt            # CMake 建置檔
├── Makefile                  # Makefile 建置檔
├── README.md                 # 使用說明
└── EXTRACTION_SUMMARY.md     # 本文件
```

## 注意事項

1. 這是簡化版本，專注於 R 值敏感度分析
2. 實際的 IR drop 計算需要完整的 SPICE 模擬（原始 NTU-3DIC 有）
3. 可以根據需要擴展功能
4. Pixel model 的記憶體管理在 main 函數中處理，實際使用時可能需要更完善的設計

