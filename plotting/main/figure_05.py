"""Figure 5: country associations and mean mines per priority species.

All 30 countries retain their source counts. Bubble area is proportional to
mine-species associations; panel b reports associations / distinct species.
Source data are read only when the figure is explicitly built.
"""

import sys
sys.dont_write_bytecode = True

from common import configure_style, table
from common import FIGURES


COLORS = {"focus": "#C15D24", "other": "#7294A4", "quiet": "#DFE8EC", "ink": "#263842"}
FOCUS = {"ECU", "VNM", "PER", "CHL", "COD"}
AREA_SCALE = 3.5


def get_data():
    import numpy as np

    df = table(5, 'Country_Summary').rename(columns={
        'Mine host iso3': 'iso3', 'Mine host country': 'country',
        'Associated mine count': 'mines', 'Priority species count': 'species',
        'Mine species association count': 'associations', 'Associations per mine': 'associations_per_mine',
    })
    for column in ["mines", "species", "associations"]:
        values = df[column].to_numpy()
        if not np.isfinite(values).all() or not (values > 0).all():
            raise ValueError(f"Invalid positive counts: {column}")
        if not np.equal(values, np.floor(values)).all():
            raise ValueError(f"Non-integer count: {column}")
    if df["iso3"].isna().any() or not df["iso3"].is_unique:
        raise ValueError("Country identifiers must be complete and unique")
    if len(df) != 30 or df["mines"].sum() != 224 or df["associations"].sum() != 548:
        raise ValueError("Country totals differ from the verified source")
    if not np.allclose(df["associations_per_mine"], df["associations"] / df["mines"]):
        raise ValueError("Existing per-mine summary is inconsistent")
    df["mines_per_species"] = df["associations"] / df["species"]
    df["display_name"] = df["country"].replace({
        "Democratic Republic of the Congo": "DR Congo",
        "Dominican Republic": "Dominican Rep.",
    })
    df["color"] = np.where(df["iso3"].isin(FOCUS), COLORS["focus"], COLORS["other"])
    return df


def label_point(ax, row, xytext, ha="left"):
    ax.annotate(
        row["display_name"],
        xy=(row["mines"], row["species"]), xytext=xytext,
        textcoords="data", ha=ha, va="center", fontsize=6.2,
        color=COLORS["focus"] if row["iso3"] in FOCUS else COLORS["ink"],
        fontweight="bold" if row["iso3"] in FOCUS else "normal",
        arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "#8B9BA3", "shrinkA": 1, "shrinkB": 5},
        zorder=6,
    )


def draw_panels(df):
    from decimal import Decimal, ROUND_HALF_UP

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    # Equal 90-mm data-frame heights; the shared layout pass trims the canvas.
    spec = FIGURES[5]
    # Match the accepted layout renderer; exports set their own 300/600 dpi.
    fig = plt.figure(figsize=(spec.width_mm / 25.4, spec.height_mm / 25.4),
                     dpi=100, facecolor="white")
    ax = fig.add_axes([20 / 200, 15 / 110, 90 / 200, 90 / 110])
    rank = fig.add_axes([150 / 200, 15 / 110, 40 / 200, 90 / 110])

    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E7ECEF", linewidth=0.5)
    ax.set_xlim(0, 57)
    ax.set_ylim(0, 34)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30])
    ax.tick_params(length=2.5, pad=3)
    for is_focus in [False, True]:
        subset = df.loc[df["iso3"].isin(FOCUS).eq(is_focus)].sort_values("associations", ascending=False)
        ax.scatter(
            subset["mines"], subset["species"],
            s=AREA_SCALE * subset["associations"], c=subset["color"],
            alpha=0.88 if is_focus else 0.72,
            edgecolors="#84431F" if is_focus else "#557987",
            linewidths=0.5, zorder=4 if is_focus else 3,
        )
    label_positions = {
        "ECU": ((9, 30.8), "left"),
        "PER": ((41, 26.0), "center"),
        "CHL": ((51, 15.3), "center"),
        "CHN": ((25, 18.2), "left"),
        "MEX": ((13, 20.4), "left"),
        "HND": ((5.0, 16.8), "left"),
        "CUB": ((5.0, 13.4), "left"),
        "VNM": ((5.0, 10.8), "left"),
        "DOM": ((7.0, 5.7), "left"),
        "PHL": ((18, 8.8), "left"),
        "COD": ((20, 3.6), "left"),
    }
    for iso, (position, align) in label_positions.items():
        label_point(ax, df.loc[df["iso3"].eq(iso)].iloc[0], position, align)
    ax.set_xlabel("Mines associated with priority species (n)", labelpad=7)
    ax.set_ylabel("Associated conservation-priority species (n)", labelpad=7)

    # Legend markers use the same area transformation as all plotted observations.
    handles = [Line2D([], [], linestyle="", marker="o", markersize=np.sqrt(AREA_SCALE * n),
                      markerfacecolor=COLORS["other"], markeredgecolor="#557987",
                      markeredgewidth=0.5, alpha=0.72, label=str(n)) for n in [10, 50, 150]]
    ax.legend(handles=handles, title="Mine–species\nassociations (n)", loc="lower right",
              bbox_to_anchor=(0.975, 0.055), ncol=1, fontsize=5.8, title_fontsize=5.8,
              borderaxespad=0, borderpad=0, handlelength=3.6, handletextpad=1.0,
              labelspacing=2.4)

    ranked = df.sort_values(["mines_per_species", "country"], ascending=[False, True]).reset_index(drop=True)
    y = np.arange(len(ranked))
    for i, row in ranked.iterrows():
        if row["iso3"] in FOCUS:
            rank.axhspan(i - 0.42, i + 0.42, color="#FBF1E9", lw=0, zorder=0)
    rank.set_axisbelow(True)
    rank.grid(axis="x", color="#E7ECEF", linewidth=0.5)
    rank.hlines(y, 0, ranked["mines_per_species"], colors=COLORS["quiet"], linewidth=1.1, zorder=1)
    rank.scatter(ranked["mines_per_species"], y, s=13, c=ranked["color"], edgecolors="none", zorder=3)
    rank.set_yticks(y)
    rank.set_yticklabels(ranked["display_name"], fontsize=5.8)
    rank.tick_params(axis="y", length=0, pad=5)
    for text, iso in zip(rank.get_yticklabels(), ranked["iso3"]):
        if iso in FOCUS:
            text.set_color(COLORS["focus"])
            text.set_fontweight("bold")
    for i, row in ranked.iterrows():
        ratio_label = format(Decimal(str(row["mines_per_species"])).quantize(
            Decimal("0.01"), rounding=ROUND_HALF_UP), "f").rstrip("0").rstrip(".")
        rank.text(row["mines_per_species"] + 0.43, i,
                  ratio_label,
                  va="center", fontsize=5.8,
                  color=COLORS["focus"] if row["iso3"] in FOCUS else "#61737D")
    rank.set_ylim(len(ranked) - 0.35, -0.8)
    rank.set_xlim(0, 16.1)
    rank.set_xticks([0, 5, 10, 15])
    rank.tick_params(axis="x", length=2.5, pad=3)
    rank.spines["left"].set_visible(False)
    rank.set_xlabel("Mean mines per priority species", labelpad=7)

    return fig, ranked


def build():
    configure_style()
    frame = get_data()
    fig, ranked = draw_panels(frame)
    return fig, {'countries': len(frame), 'mines': int(frame.mines.sum()),
                 'associations': int(frame.associations.sum())}


if __name__ == '__main__':
    from common import run
    run(5, build)
