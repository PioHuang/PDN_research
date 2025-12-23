### PDN Sensitivity Analysis (VoltSpot-compatible steady-state IR drop)

這個專案把 VoltSpot 的「steady-state PDN」抽成一個 C++ flow：把 PDN 建成電阻網路（電路），把 pads（供電接點）與 loads（負載電流）放上去，解線性方程得到每個 grid node 的 VDD/GND 電壓，最後用 VoltSpot 相同的定義輸出 `ir_drop.gridIR`。

---

### TL;DR

```bash
cd build
cmake ..
make

./pdn_load_example \
  --flp ../voltspot/example.flp \
  --config ../voltspot/pdn.config \
  --padloc ../voltspot/example.vgrid.padloc \
  --ptrace ../voltspot/example.ptrace \
  --mlcf ../voltspot/example.mlcf

python3 ../tools/visualize_pdn_loads.py \
  --branches ../out_voltspot/pdn_with_loads_branches.csv \
  --nodes ../out_voltspot/pdn_with_loads_branches_nodes.csv \
  --out ../out_voltspot/pdn_visualization.png

python3 ../tools/visualize_ir_drop.py \
  --gridir ../out_voltspot/ir_drop.gridIR \
  --out ../out_voltspot/ir_drop_visualization.png
```

---

### 名詞

- **grid / virtual grid**：把整個 chip 切成 `rows x cols` 的格點。每個格點就是一個「node」。
- **physical layer（網格層數）**：金屬層/PDN 層的意思。本專案 `PDNNetwork(rows, cols, layers)` 的 `layers` 是這個。
  - 範例目前是 **1 層**（Layer 0）。
- **rail / net（VDD / GND）**：每個 physical layer 上有兩張網：**VDD rail** + **GND rail**。
  - 所以未知數不是 `rows*cols`，而是 **`2 * layers * rows * cols`**。
  - 例子：`rows=73, cols=73, layers=1` → unknowns = `2*1*73*73 = 10658`（輸出也看得到 Nodes: 10658）。
- **branch**：兩個 node 之間的電阻連線（上下左右、或跨層）。
- **pad（供電接點）**：晶片連到外部供電/封裝的接點（C4 bump / power pad / ground pad）。
  - **VDD pad**：把附近電壓「拉向」VDD
  - **GND pad**：把附近電壓「拉向」GND
  - 重點：pad 不是理想電壓源，而是透過 **pad 電阻 `PDN_padR`** 連到供應電壓，電壓不會被硬 clamp。
- **load（負載）**：從 power trace 轉成的電流消耗（A）。VoltSpot 模型是「負載電流在 VDD 與 GND 之間流」：
  - VDD node：注入 `-Iload`
  - GND node：注入 `+Iload`

---

### flow 介紹

執行檔：`src/examples/pdn_load_example.cpp`

1. **讀 FLP**（`--flp`）
   - 取得 chip 寬高 + 每個 unit 的矩形位置（m）。
2. **讀 VoltSpot config**（`--config`）
   - 重要參數：`PDN_padpitch`, `PDN_grid_intv`, `PDN_padR`, `PDN_padconfig`, `vdd`, `gnd`。
3. **建 virtual grid**（`VoltSpotVirtualGrid`）
   - 算出 `rows/cols`, `dx/dy`, 等效 `rx/ry`。
4. **建 PDNNetwork**（電阻網路骨架）
   - 在每個 grid (r,c,l) 建 node，建立相鄰連線電阻（branch）。
5. **放 pads（供電點）**
   - VoltSpot 行為（你要對齊正解一定要一致）：
     - **`PDN_padconfig == 0`**：從 `--padloc` 讀 V/G pad 位置（例如 `voltspot/example.vgrid.padloc`）
     - **`PDN_padconfig != 0`**：忽略 `--padloc`，把 **所有 pad seats** 都放滿 P/G（本專案已實作）
   - 例子 `PDN_grid_intv=2`、grid 73x73 → pad grid 為 37x37 → pad seats 共 1369。
6. **放 loads（負載電流）**
   - 從 `--ptrace` 讀功耗 W，轉成 \(I=P/V\)。
   - 依 unit 與 grid cell 的重疊面積分配到格點，並同時寫到 VDD/GND rail（VDD:-I、GND:+I 的模型靠 solver RHS 來完成）。
7. **解 steady-state（Solver）**
   - 組 conductance matrix \(G\) 與 RHS \(b\)，解 \(Gx=b\)。
8. **輸出**
   - `out_voltspot/ir_drop.gridIR`：VoltSpot 格式（`#Layer:N` + `col row drop%`）
   - `out_voltspot/pdn_with_loads_branches*.csv`：給視覺化用
   - `tools/visualize_*`：把分佈畫成 png

---

### 座標系

VoltSpot 對外輸出（padloc/gridIR）用的座標習慣是：**(0,0) 在左下角**。  
但實作上常會用 `nr-i-1` 這種方式把 internal row index 翻轉後再輸出。

本專做兩件事跟 VoltSpot 一致：

- **IR drop 的 `gridIR` 輸出 row 方向一致**（VoltSpot 的 row 翻轉規則）
- **load mapping / PDN visualization 的 row 方向一致**（避免「內部 row」與「輸出 row」混用造成整張圖翻面）

---

### Padconfig

- **pad 分佈模式**：`PDN_padconfig=1`（補滿所有 pad seats）vs `PDN_padconfig=0`（用稀疏 padloc）
- **padR **：VoltSpot pdn.config 裡 `PDN_padR=10e-3`

---

### Solver 在做什麼？（大 KCL）

可以把 solver 想成：對每個 node 寫 KCL，最後變成一個超大的線性方程 \(Gx=b\)。

下面用一個最小可算的例子（**只看 VDD rail、兩個 node**）示範：

- node0 ↔ node1 有電阻 (R = 1, g = 1)
- node0 是 VDD pad，pad 電阻 (R_p = 0.5, g_p = 2)
- 供應電壓 (1.0 V)
- node1 有負載電流 (0.2 A)

未知數：

```text
x = [V0, V1]^T
```

KCL：

```text
node0: g(V0 - V1) + gp(V0 - Vdd) = 0
       (g+gp)V0 - gV1 = gp*Vdd
       (1+2)V0 - 1*V1 = 2*1   =>  3V0 - V1 = 2

node1: g(V1 - V0) = -Iload
       -gV0 + gV1 = -Iload
       -1*V0 + 1*V1 = -0.2    =>  -V0 + V1 = -0.2
```

組成矩陣：

```text
|  3  -1 | | V0 | = |  2   |
| -1   1 | | V1 |   | -0.2 |
```

解：

```text
由第二式：V1 = V0 - 0.2
代回第一式：3V0 - (V0 - 0.2) = 2
          2V0 + 0.2 = 2  => V0 = 0.9
所以 V1 = 0.7
```

- **對角線**：
  - node0: g+g_p = 1+2 = 3 (鄰居的電導總和 +（pad 的 \(1/R_p\))
  - node1: g = 1
- **非對角線**：跟鄰居連線的 \(-g\)
- **RHS**：
  - node0: 有 pad --> +Vdd/R_p = 1/0.5 = 2
  - node1: 有 load --> -Iload = -0.2

> 真實情況只是把「兩個 node」擴成「73x73 的所有 node」，再把 VDD/GND rail 都一起解（總 unknowns 變 10658），矩陣更大而已。

---

### 輸出檔案說明

- **`out_voltspot/ir_drop.gridIR`**：VoltSpot 格式 IR drop（每層一段 `#Layer:N`，每行 `col row drop%`）
- **`out_voltspot/ir_drop_visualization.png`**：IR drop heatmap
- **`out_voltspot/pdn_with_loads_branches.csv`**：branch
- **`out_voltspot/pdn_with_loads_branches_nodes.csv`**：node（座標、pgNet、voltage、currentLoad、isPad）
- **`out_voltspot/pdn_visualization.png`**：pads / loads 分佈圖

---

### FAQ

- **GND pad 不能用 `voltage==0` 去判斷是不是 pad**：所以本專案用 `isPad` 明確標記。
- **load 不能只打到 VDD**：VoltSpot steady-state 的負載是 VDD:-I、GND:+I（兩邊都要一致）。
- **padconfig 跑錯會讓 max 差很多**：`PDN_padconfig=1`（補滿）跟 `=0`（padloc 稀疏）是完全不同的供電能力。
