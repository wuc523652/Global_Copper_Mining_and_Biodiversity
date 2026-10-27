"""Final uncapped screening-rule comparison, with all 891 affected species."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

MAIN = "Main_SreekarWatson2026"
ALT = "Traditional_RodriguesLewin"
WIDTH_MM = 183.0
SCHEMES = {
    MAIN: {"low": 200.0, "floor": 0.30, "color": "#D55E00", "label": "Main rule"},
    ALT: {"low": 1000.0, "floor": 0.10, "color": "#0072B2", "label": "Alternative rule"},
}

def screening_threshold(area, scheme, *, apply_cap):
    area = np.asarray(area, dtype=float)
    if not np.isfinite(area).all() or not (area > 0).all():
        raise ValueError("Range areas must be finite and strictly positive.")
    fraction = np.clip(
        (np.log10(area) - np.log10(scheme["low"]))
        / (np.log10(250000.0) - np.log10(scheme["low"])),
        0.0,
        1.0,
    )
    target = 1.0 - fraction * (1.0 - scheme["floor"])
    if apply_cap:
        target = np.minimum(target, 1000000.0 / area)
    return 100.0 * (1.0 - target)

def style():
    # Preserve existing typography rather than applying an unrelated redesign.
    mpl.rcParams.update({
        "font.family": ["Times New Roman", "Arial"],
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.2,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7.0,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "legend.fontsize": 5.8,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "legend.frameon": False,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })

def composition_panel(ax, data, field, order, colors, labels, rule_labels):
    bottom = np.zeros(2)
    for group in order:
        values = []
        for prefix in [MAIN, ALT]:
            chosen = data.loc[data[f"{prefix}_Priority"]]
            values.append(100.0 * chosen[field].eq(group).sum() / len(chosen))
        ax.bar(np.arange(2), values, bottom=bottom, width=0.62,
               color=colors[group], edgecolor="white", linewidth=0.4,
               label=labels.get(group, group))
        bottom += values
    assert np.allclose(bottom, 100.0)
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(rule_labels)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Composition of conservation-priority set (%)")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=5.2)

def build_plot(data, comparison):
    style()
    rule_labels = ["Main rule", "Alternative rule"]
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_MM / 25.4, 6.65))
    ax_a, ax_b, ax_c, ax_d = axes.flat
    ax_a.scatter(data["Range_Area_km2"], data["Impact_Percent"], s=7,
                 color="#777777", alpha=0.28, linewidth=0)
    curve_x = np.logspace(
        np.log10(max(1, data["Range_Area_km2"].min() * 0.75)),
        np.log10(max(15000000, data["Range_Area_km2"].max() * 1.2)), 500)
    for i, scheme in enumerate(SCHEMES.values()):
        curve_y = screening_threshold(curve_x, scheme, apply_cap=False)
        assert np.all(np.diff(curve_x) > 0)
        assert np.all(np.diff(curve_y) >= -1e-10)
        label = scheme["label"]
        ax_a.plot(curve_x, curve_y, color=scheme["color"], linewidth=1.5, label=label)
    ax_a.set_xscale("log")
    ax_a.tick_params(axis="x", labelsize=7.2)
    ax_a.set_ylim(-2, 102)
    ax_a.set_xlabel("Global habitat area (km²)")
    ax_a.set_ylabel("Proportion of habitat affected (%)")
    ax_a.legend(loc="upper right", fontsize=5.2)

    stable = int(comparison["Stable_Core_N"])
    ax_b.bar([0, 1], [stable, stable], color="#4d4d4d", width=0.62, label="Shared")
    ax_b.bar(0, int(comparison["Current_Only_N"]), bottom=stable,
             color="#D55E00", width=0.62,
             label="Main rule only")
    ax_b.bar(1, int(comparison["Traditional_Only_N"]), bottom=stable,
             color="#0072B2", width=0.62,
             label="Alternative rule only")
    ax_b.set_xticks([0, 1])
    ax_b.set_xticklabels(rule_labels)
    ax_b.set_ylabel("Conservation-priority species (n)")
    ax_b.text(0.97, 0.97, f"Jaccard = {float(comparison['Jaccard_Index']):.3f}",
              transform=ax_b.transAxes, ha="right", va="top")
    ax_b.legend(loc="lower left", bbox_to_anchor=(0.0, -0.24), ncol=1, fontsize=5.2)
    composition_panel(
        ax_c, data, "Taxon", ["Amphibians", "Birds", "Fish", "Mammals", "Reptiles"],
        {"Amphibians": "#D55E00", "Birds": "#E69F00", "Fish": "#56B4E9",
         "Mammals": "#009E73", "Reptiles": "#7A5195"},
        {"Fish": "Freshwater\nfishes"}, rule_labels)
    composition_panel(
        ax_d, data, "Current_Range_Class", ["Narrow", "Intermediate", "Widespread"],
        {"Narrow": "#e28740", "Intermediate": "#84cfbb", "Widespread": "#7fc7ef"},
        {"Narrow": "Narrow-habitat", "Intermediate": "Intermediate-\nhabitat"}, rule_labels)
    for label, ax in zip("abcd", axes.flat):
        ax.text(-0.16, 1.04, label, transform=ax.transAxes, ha="left", va="top",
                 fontsize=8.0, fontweight="bold", clip_on=False, zorder=100)
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.3)
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.10, right=0.90, top=0.96, bottom=0.09,
                        wspace=0.43, hspace=0.35)
    return fig
