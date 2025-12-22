# make

cd build
cmake ..
make

# execute

./pdn_load_example \
--flp ../voltspot/example.flp \
--config ../voltspot/pdn.config \
--padloc ../voltspot/example.vgrid.padloc \
--ptrace ../voltspot/example.ptrace \
--mlcf ../voltspot/example.mlcf

# Step 1-2: 解析 FLP 和 config 文件

從 FLP 文件獲取晶片尺寸
從 config 文件獲取 PDN 參數

# Step 3: 建立 virtual grid（已整合）

計算 virtual grid 大小（rows, cols）
計算 grid spacing（dx, dy）
計算 edge resistances（rx, ry）

# Step 4: 建立 PDNNetwork 並設置 pixel models

使用 virtual grid 的大小建立 network
為每個 grid cell 創建 pixel model（使用計算出的 rx, ry）

# Step 5: 建立網絡結構

連接所有節點和分支

# Step 6: 載入 voltage sources（可選）

# Step 7: 載入 current loads（可選）

# Step 8: 導出結果

# Step 9: visualize

python3 tools/visualize_pdn_loads.py \
--branches out_voltspot/pdn_with_loads_branches.csv \
--nodes out_voltspot/pdn_with_loads_branches_nodes.csv \
--out out_voltspot/pdn_visualization.png
