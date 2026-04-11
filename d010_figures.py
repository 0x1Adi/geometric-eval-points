#!/usr/bin/env python3
"""
D010 Figure Generator — Publication-quality plots for paper
============================================================
Run AFTER: sage d010_comprehensive.sage
Usage:     python3 d010_figures.py
Output:    d010_results/figures/*.pdf + *.png

Author: Aditya | ORCID: 0009-0004-8454-9777
"""
import csv, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

DATA = "d010_results"
FIGS = f"{DATA}/figures"
os.makedirs(FIGS, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.labelsize': 10, 'axes.titlesize': 11,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})
COLORS = ['#0072B2', '#D55E00', '#009E73', '#E69F00', '#CC79A7', '#56B4E9', '#F0E442']

def load_csv(name):
    with open(f"{DATA}/{name}") as f:
        reader = csv.DictReader(f)
        return list(reader)

def savefig(fig, name):
    fig.savefig(f"{FIGS}/{name}.pdf")
    fig.savefig(f"{FIGS}/{name}.png", dpi=300)
    plt.close(fig)
    print(f"  Saved {name}.pdf + .png")

# ============================================================
# FIGURE 1: ||w||_1 vs t (log scale, 3 methods)
# ============================================================
print("Generating Figure 1: Norm vs threshold t")
data = load_csv("fig1_norm_vs_t.csv")
ts = [int(r['t']) for r in data]
consec = [float(r['consecutive']) for r in data]
geo = [float(r['geometric_b10']) for r in data]
anchor = [float(r['anchor_2_256']) for r in data]

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.semilogy(ts, consec, 'o-', color=COLORS[0], markersize=3, label='Consecutive {1,2,...,t}', linewidth=1.5)
ax.semilogy(ts, geo, 's-', color=COLORS[1], markersize=3, label='Geometric {10,10²,...,10ᵗ}', linewidth=1.5)
ax.semilogy(ts, anchor, '^-', color=COLORS[2], markersize=3, label='Anchor {1, M/t, 2M/t,...}', linewidth=1.5)
ax.axhline(y=1.25, color=COLORS[1], linestyle='--', alpha=0.5, linewidth=0.8)
ax.axhline(y=1.0, color=COLORS[2], linestyle='--', alpha=0.5, linewidth=0.8)
ax.text(32, 1.35, 'b/(b−2) = 1.25', fontsize=7, color=COLORS[1])
ax.text(32, 1.05, '1.0', fontsize=7, color=COLORS[2])
ax.set_xlabel('Threshold t (number of parties)')
ax.set_ylabel('||w||₁ (Lagrange weight L1-norm)')
ax.set_title('Noise amplification: Consecutive vs Geometric vs Anchor')
ax.legend(loc='upper left', frameon=False)
ax.set_xlim(2, 34)
ax.set_ylim(0.8, 1e11)
savefig(fig, "fig1_norm_vs_t")

# ============================================================
# FIGURE 2: Individual weight decay
# ============================================================
print("Generating Figure 2: Weight decay pattern")
data = load_csv("claim_2_5_weight_decay.csv")
indices = [int(r['index']) for r in data]
weights = [float(r['abs_weight']) for r in data]

fig, ax = plt.subplots(figsize=(4.5, 3))
ax.semilogy(indices, weights, 'o-', color=COLORS[0], markersize=5, linewidth=1.5)
ax.set_xlabel('Weight index i')
ax.set_ylabel('|wᵢ| (absolute weight)')
ax.set_title('Individual weight magnitudes (geometric b=10, t=15)')
ax.set_xticks(range(1, 16, 2))
# Annotate first two
ax.annotate(f'|w₁| = {weights[0]:.4f}', xy=(1, weights[0]),
            xytext=(3, weights[0]*2), fontsize=7,
            arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
ax.annotate(f'|w₂| ≈ 10⁻¹', xy=(2, weights[1]),
            xytext=(4, weights[1]*5), fontsize=7,
            arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))
savefig(fig, "fig2_weight_decay")

# ============================================================
# FIGURE 3: Noise comparison (grouped bar chart)
# ============================================================
print("Generating Figure 3: Noise amplification comparison")
data = load_csv("claim_3_1_noise_simulation.csv")
ts = [int(r['t']) for r in data]
c_noise = [float(r['consecutive_noise']) for r in data]
g_noise = [float(r['geometric_noise']) for r in data]
a_noise = [float(r['anchor_noise']) for r in data]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3), gridspec_kw={'width_ratios': [1, 1]})

# Left: log-scale bar chart
x = np.arange(len(ts))
w = 0.25
ax1.bar(x - w, c_noise, w, label='Consecutive', color=COLORS[0])
ax1.bar(x, g_noise, w, label='Geometric (b=10)', color=COLORS[1])
ax1.bar(x + w, a_noise, w, label='Anchor', color=COLORS[2])
ax1.set_yscale('log')
ax1.set_xlabel('Threshold t')
ax1.set_ylabel('Total noise after recombination')
ax1.set_xticks(x)
ax1.set_xticklabels([str(t) for t in ts])
ax1.legend(fontsize=6, frameon=False)
ax1.set_title('(a) Noise amplification (log scale)')

# Right: improvement factor
imp_g = [float(r['improvement_geometric']) for r in data]
imp_a = [float(r['improvement_anchor']) for r in data]
ax2.semilogy(ts, imp_g, 's-', color=COLORS[1], label='Geometric improvement', markersize=5)
ax2.semilogy(ts, imp_a, '^-', color=COLORS[2], label='Anchor improvement', markersize=5)
ax2.set_xlabel('Threshold t')
ax2.set_ylabel('Improvement factor over consecutive')
ax2.legend(fontsize=7, frameon=False)
ax2.set_title('(b) Improvement factor')
fig.tight_layout()
savefig(fig, "fig3_noise_comparison")

# ============================================================
# FIGURE 4: Convergence to b/(b-2)
# ============================================================
print("Generating Figure 4: Convergence to b/(b-2)")
data = load_csv("fig4_convergence.csv")
fig, ax = plt.subplots(figsize=(5, 3.5))

for idx, b_val in enumerate([3, 5, 10, 50, 100]):
    subset = [r for r in data if int(r['b']) == b_val]
    ts_b = [int(r['t']) for r in subset]
    norms_b = [float(r['norm']) for r in subset]
    pred = float(subset[0]['predicted'])
    ax.plot(ts_b, norms_b, 'o-', color=COLORS[idx], markersize=3,
            label=f'b={b_val} (→{pred:.3f})', linewidth=1.2)
    ax.axhline(y=pred, color=COLORS[idx], linestyle=':', alpha=0.3, linewidth=0.7)

ax.set_xlabel('Threshold t')
ax.set_ylabel('||w||₁')
ax.set_title('Convergence of geometric norm to b/(b−2)')
ax.legend(fontsize=7, frameon=False, ncol=2)
ax.set_xlim(3, 29)
ax.set_ylim(0.95, 3.5)
savefig(fig, "fig4_convergence")

# ============================================================
# FIGURE 5: Point family comparison (bar chart)
# ============================================================
print("Generating Figure 5: Point family comparison")
data = load_csv("claim_2_4_family_comparison.csv")
names = [r['family'] for r in data]
norms = [float(r['norm']) for r in data]

fig, ax = plt.subplots(figsize=(6, 3.5))
colors_bar = [COLORS[0]]*5 + [COLORS[1]]*2 + [COLORS[3]] + [COLORS[2]]
bars = ax.barh(range(len(names)), norms, color=colors_bar, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel('||w||₁ (L1-norm)')
ax.set_title('Lagrange weight norm by evaluation point family (t=10)')
ax.set_xscale('log')
# Annotate values
for i, (n, bar) in enumerate(zip(norms, bars)):
    ax.text(n * 1.3, i, f'{n:.2f}', va='center', fontsize=7)
ax.invert_yaxis()
savefig(fig, "fig5_family_comparison")

# ============================================================
# FIGURE 6: ML-DSA feasibility
# ============================================================
print("Generating Figure 6: ML-DSA feasibility")
data = load_csv("claim_3_2_mldsa_feasibility.csv")

fig, ax = plt.subplots(figsize=(5, 3))
for idx, variant in enumerate(["ML-DSA-44", "ML-DSA-65", "ML-DSA-87"]):
    subset = [r for r in data if r['variant'] == variant]
    ts_v = [int(r['t']) for r in subset]
    headroom = [float(r['headroom']) for r in subset]
    ax.plot(ts_v, headroom, 'o-', color=COLORS[idx], markersize=4,
            label=variant, linewidth=1.5)

ax.axhline(y=1, color='red', linestyle='--', alpha=0.5, linewidth=1, label='FIPS bound (headroom=1)')
ax.set_xlabel('Threshold t')
ax.set_ylabel('Headroom factor (bound / amplified noise)')
ax.set_title('ML-DSA FIPS 204 feasibility with geometric points (b=10)')
ax.set_yscale('log')
ax.legend(fontsize=7, frameon=False)
savefig(fig, "fig6_mldsa_feasibility")

# ============================================================
# FIGURE 7: Practical parameter space
# ============================================================
print("Generating Figure 7: Parameter space")
data = load_csv("claim_3_3_practical_params.csv")

fig, ax = plt.subplots(figsize=(5, 3.5))
for idx, lq in enumerate([200, 400, 800]):
    subset = [r for r in data if int(r['log2_q']) == lq]
    bs = [int(r['b']) for r in subset]
    max_ts = [int(r['max_threshold']) for r in subset]
    ax.plot(bs, max_ts, 'o-', color=COLORS[idx], markersize=5,
            label=f'q = 2^{lq}', linewidth=1.5)

ax.set_xlabel('Base b')
ax.set_ylabel('Maximum threshold t')
ax.set_title('Maximum parties vs base b for different FHE moduli')
ax.set_xscale('log')
ax.legend(fontsize=8, frameon=False)
ax.axhline(y=13, color='gray', linestyle=':', alpha=0.5)
ax.text(1200, 15, 'Zama TKMS limit (13)', fontsize=7, color='gray')
savefig(fig, "fig7_parameter_space")

# ============================================================
# FIGURE 8: Comparison table (visual)
# ============================================================
print("Generating Figure 8: Approach comparison table")
fig, ax = plt.subplots(figsize=(7, 3))
ax.axis('off')
table_data = [
    ['Standard Shamir', 'O(2ᵗ)', '1', 'O(t) bits', '~16'],
    ['{0,1}-LSSS', 'O(t)', 'O(N⁴·²)', 'O(log N)', 'Any'],
    ['Denom. clearing', 'O((N!)²)', '1', 'O(N log N)', 'Any'],
    ['Noise flooding', 'O(1)+2^λ', '1', 'O(λ)', 'Any'],
    ['Masking', 'O(1)', '1+PRF', 'None', 'Any'],
    ['Geometric (ours)', '1.25', '1', '0.32 bits', '120'],
    ['Anchor (ours)', '≈1.0', '1', '≈0 bits', '120+'],
]
col_labels = ['Method', 'Noise factor', 'Shares/party', 'Modulus overhead', 'Max t']
table = ax.table(cellText=table_data, colLabels=col_labels, loc='center',
                 cellLoc='center', colColours=['#e6e6e6']*5)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1.0, 1.4)
# Highlight our rows
for j in range(5):
    table[6, j].set_facecolor('#d4edda')
    table[7, j].set_facecolor('#d4edda')
ax.set_title('Comparison of threshold noise reduction approaches', fontsize=10, pad=20)
savefig(fig, "fig8_comparison_table")

print(f"\nAll {8} figures saved to {FIGS}/")
print("Formats: PDF (vector) + PNG (300 DPI)")
