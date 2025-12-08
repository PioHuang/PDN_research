#!/bin/bash
# Demo script for PDN visualization

echo "=== PDN Visualization Demo ==="
echo ""

# Check if Python and matplotlib are available
if ! python3 -c "import matplotlib" 2>/dev/null; then
    echo "Error: matplotlib is not installed"
    echo "Install it with: pip install matplotlib numpy"
    exit 1
fi

echo "1. Generating grid visualization (5x5)..."
python3 visualize_pdn.py --rows 5 --cols 5 --output demo_grid.png --no-show

echo "2. Generating grid visualization (10x10)..."
python3 visualize_pdn.py --rows 10 --cols 10 --output demo_grid_10x10.png --no-show

echo "3. Generating Rx heatmap..."
python3 visualize_pdn.py --type heatmap --direction x --rows 10 --cols 10 \
    --output demo_rx_heatmap.png --no-show

echo "4. Generating Ry heatmap..."
python3 visualize_pdn.py --type heatmap --direction y --rows 10 --cols 10 \
    --output demo_ry_heatmap.png --no-show

echo "5. Generating graph visualization..."
python3 visualize_pdn.py --type graph --rows 5 --cols 5 \
    --output demo_graph.png --no-show

echo ""
echo "=== All visualizations generated! ==="
echo "Generated files:"
ls -lh demo_*.png 2>/dev/null

echo ""
echo "To view the visualizations, open the PNG files or run:"
echo "  python3 visualize_pdn.py --rows 5 --cols 5"

