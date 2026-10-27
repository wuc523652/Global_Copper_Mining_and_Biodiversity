"""Four-column alluvial panels with linear area-proportional ribbon widths.
"""

from __future__ import annotations

from functools import lru_cache
from math import floor, log10

from collections import defaultdict
from dataclasses import dataclass
from math import fsum, isclose, isfinite


SELECTED_OWNERS = ("CAN", "CHN", "GBR", "USA", "AUS")
SELECTED_HOSTS = ("PER", "CHL", "PHL", "ECU", "MEX", "CHN", "COD", "BRA", "ZAF")
CONTINENTS = ("North America", "South America", "Europe", "Asia", "Oceania", "Africa")
VALUE_FIELD = "attributed_exposure_area_km2"
COUNTERPART_THRESHOLDS_KM2 = {"owner": 10000.0, "host": 1000.0}
OTHERS_CODE = "__OTHERS__"


@dataclass(frozen=True)
class Node:
    key: tuple[int, str]
    name: str
    continent: str
    value: float


@dataclass(frozen=True)
class Link:
    source: tuple[int, str]
    target: tuple[int, str]
    value: float


@dataclass(frozen=True)
class FlowModel:
    perspective: str
    selected: tuple[str, ...]
    columns: tuple[tuple[Node, ...], ...]
    links: tuple[Link, ...]
    total: float
    source_rows: int
    selected_rows: int
    counterpart_threshold: float
    original_counterparts: tuple[Node, ...]
    grouped_counterparts: tuple[Node, ...]

    def summary(self):
        return {
            "perspective": self.perspective,
            "selected_left_iso3": list(self.selected),
            "source_rows": self.source_rows,
            "selected_cross_border_pairs": self.selected_rows,
            "column_node_counts": [len(column) for column in self.columns],
            "column_area_totals_km2": [fsum(n.value for n in col) for col in self.columns],
            "displayed_area_km2": self.total,
            "country_totals_km2": {
                side: {n.key[1]: n.value for n in self.columns[index]}
                for side, index in [("left", 1), ("right", 2)]
            },
            "continent_totals_km2": {
                side: {n.name: n.value for n in self.columns[index]}
                for side, index in [("left", 0), ("right", 3)]
            },
            "counterpart_policy": "right-country panel totals strictly below the threshold grouped as Others; Others placed last",
            "counterpart_threshold_km2": self.counterpart_threshold,
            "original_counterpart_count": len(self.original_counterparts),
            "individually_displayed_counterpart_count": len(self.original_counterparts) - len(self.grouped_counterparts),
            "grouped_counterpart_count": len(self.grouped_counterparts),
            "others_total_km2": fsum(n.value for n in self.grouped_counterparts),
            "others_members": [{"iso3": n.key[1], "country": n.name, "continent": n.continent,
                                "area_km2": n.value} for n in self.grouped_counterparts],
            "others_continent_policy": "original country-continent contributions retained; Others is not assigned to a single continent",
            "area_scope": "selected panel only, not the global total",
        }


def _text(record, field):
    value = str(record.get(field, "")).strip()
    if not value or value.lower() in {"nan", "none", "null"}:
        raise ValueError(f"Missing country metadata: {field}")
    return value


def build_model(records, perspective, *, counterpart_threshold=None):
    """Aggregate existing bilateral results; do not recalculate attribution."""
    if perspective not in {"owner", "host"}:
        raise ValueError("Perspective must be owner or host")
    threshold = (COUNTERPART_THRESHOLDS_KM2[perspective] if counterpart_threshold is None
                 else float(counterpart_threshold))
    if not isfinite(threshold) or threshold < 0:
        raise ValueError("Counterpart threshold must be finite and nonnegative")
    left, right = ("investor", "resource") if perspective == "owner" else ("resource", "investor")
    selected = SELECTED_OWNERS if perspective == "owner" else SELECTED_HOSTS
    names, continents = {}, {}
    values = defaultdict(list)
    pairs, present = set(), set()
    source_rows = selected_rows = 0
    for record in records:
        source_rows += 1
        value = float(record[VALUE_FIELD])
        domestic = float(record["is_domestic"])
        if not isfinite(value) or value < 0 or domestic not in {0.0, 1.0}:
            raise ValueError(f"Invalid area or domestic flag at source row {source_rows}")
        left_code, right_code = (_text(record, f"{side}_iso3").upper() for side in (left, right))
        if bool(domestic) != (left_code == right_code):
            raise ValueError(f"Domestic flag conflicts with country codes at row {source_rows}")
        if domestic or value == 0 or left_code not in selected:
            continue
        pair = (left_code, right_code)
        if pair in pairs:
            raise ValueError(f"Duplicate bilateral country pair: {pair}")
        pairs.add(pair)
        present.add(left_code)
        selected_rows += 1
        for index, side, code in [(1, left, left_code), (2, right, right_code)]:
            key = (index, code)
            name = _text(record, f"{side}_country")
            continent = _text(record, f"{side}_continent")
            if continent not in CONTINENTS:
                raise ValueError(f"Unknown continent for {key}: {continent}")
            if key in names and (names[key], continents[key]) != (name, continent):
                raise ValueError(f"Inconsistent country metadata for {key}")
            names[key], continents[key] = name, continent
        source, target = (1, left_code), (2, right_code)
        l_cont, r_cont = (0, continents[source]), (3, continents[target])
        for key in (l_cont, r_cont):
            names[key] = continents[key] = key[1]
        for edge in [(l_cont, source), (source, target), (target, r_cont)]:
            values[edge].append(value)
    if set(selected) != present:
        raise ValueError(f"Selected countries without positive cross-border area: {sorted(set(selected) - present)}")
    # Apply the threshold to each right country's total across the selected left
    # countries, never to individual bilateral links or to rounded display text.
    right_totals = defaultdict(list)
    for (source, target), contributions in values.items():
        if source[0] == 1:
            right_totals[target].extend(contributions)
    original = tuple(sorted((Node(key, names[key], continents[key], fsum(v))
                             for key, v in right_totals.items()), key=lambda n: (-n.value, n.key[1])))
    grouped = tuple(n for n in original if n.value < threshold)
    grouped_keys = {n.key for n in grouped}
    if grouped:
        others = (2, OTHERS_CODE)
        collapsed = defaultdict(list)
        for (source, target), contributions in values.items():
            source = others if source in grouped_keys else source
            target = others if target in grouped_keys else target
            collapsed[(source, target)].extend(contributions)
        values = collapsed
        for key in grouped_keys:
            del names[key], continents[key]
        names[others], continents[others] = "Others", "Mixed"
    links = tuple(Link(s, t, fsum(v)) for (s, t), v in sorted(values.items()))
    incoming, outgoing = defaultdict(list), defaultdict(list)
    for link in links:
        outgoing[link.source].append(link.value)
        incoming[link.target].append(link.value)
    nodes = []
    for key in sorted(names):
        if key[0] in (1, 2) and not isclose(fsum(incoming[key]), fsum(outgoing[key]), rel_tol=1e-12):
            raise ValueError(f"Area conservation failed at {key}")
        total = fsum(outgoing[key] if key[0] < 3 else incoming[key])
        nodes.append(Node(key, names[key], continents[key], total))
    columns = tuple(tuple(sorted((n for n in nodes if n.key[0] == i),
                                 key=lambda n: (n.key[1] == OTHERS_CODE, -n.value, n.key[1]))) for i in range(4))
    total = fsum(node.value for node in columns[0])
    if not all(isclose(fsum(n.value for n in col), total, rel_tol=1e-12) for col in columns):
        raise ValueError("Four-column totals differ")
    return FlowModel(perspective, selected, columns, links, total, source_rows, selected_rows,
                     threshold, original, grouped)


def major_keys(model):
    """All selected left countries; five largest owners on the host view."""
    result = {node.key for node in model.columns[1]}
    if model.perspective == "host":
        individual = [node for node in model.columns[2] if node.key[1] != OTHERS_CODE]
        result.update(node.key for node in individual[:5])
    return result


def allocate_layout(model, height_points, y_min=0.08, y_max=0.96):
    """Fit fixed-size labels using blank slots and one shared linear area scale.

    Coordinates increase downwards. A node remains centered on its label slot;
    source and destination strips partition its actual area-proportional bar.
    """
    if not isfinite(height_points) or height_points <= 0 or not 0 <= y_min < y_max <= 1:
        raise ValueError("Invalid panel dimensions")
    major = major_keys(model)
    minimum = {n.key: (12.0 if n.key[0] == 0 else 20.0 if n.key[0] == 3
                      else 17.0 if n.key in major else 6.2) / height_points
               for col in model.columns for n in col}
    gap = 0.30 / height_points
    available = y_max - y_min

    def used(column, scale):
        return fsum(max(n.value * scale, minimum[n.key]) for n in column) + gap * (len(column) - 1)

    if any(used(col, 0) >= available for col in model.columns):
        raise ValueError("Panel is too short for 5.5 pt labels; increase its height, do not drop countries")
    low, high = 0.0, available / model.total
    for _ in range(80):
        middle = (low + high) / 2
        if all(used(col, middle) <= available for col in model.columns):
            low = middle
        else:
            high = middle
    scale = low
    bounds, slots, rank = {}, {}, {}
    for column in model.columns:
        cursor = y_min
        for index, node in enumerate(column):
            width = node.value * scale
            height = max(width, minimum[node.key])
            slots[node.key] = (cursor, cursor + height)
            top = cursor + (height - width) / 2
            bounds[node.key] = (top, top + width)
            rank[node.key] = index
            cursor += height + gap
    source_bands, target_bands = {}, {}
    for endpoint, opposite, result in [("source", "target", source_bands), ("target", "source", target_bands)]:
        groups = defaultdict(list)
        for link in model.links:
            groups[getattr(link, endpoint)].append(link)
        for key, links in groups.items():
            cursor = bounds[key][0]
            for link in sorted(links, key=lambda edge: rank[getattr(edge, opposite)]):
                bottom = cursor + link.value * scale
                result[(link.source, link.target)] = (cursor, bottom)
                cursor = bottom
            if not isclose(cursor, bounds[key][1], abs_tol=1e-12):
                raise ValueError(f"Ribbon endpoints fail to partition node {key}")
    return {"area_scale": scale, "node_bounds": bounds, "label_slots": slots,
            "source_bands": source_bands, "target_bands": target_bands}



from common import INPUTS
ASSET_DIR = INPUTS / 'icons' if INPUTS is not None else None
STYLE = {'continent_colors': {'North America': '#4076B6',
                      'South America': '#20683A',
                      'Europe': '#F5A238',
                      'Asia': '#A61D31',
                      'Oceania': '#1D7978',
                      'Africa': '#652F7A'},
 'continent_icons': {'North America': 'image2.png',
                     'South America': 'image3.png',
                     'Oceania': 'image4.png',
                     'Asia': 'image5.png',
                     'Europe': 'image6.png',
                     'Africa': 'image7.png'},
 'country_icons': {'AUS': [{'file': 'flag_AUS_original.png',
                            'crop_pixels': [141, 123, 3595, 3577],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'PER': [{'file': 'flag_PER_original.png',
                            'crop_pixels': [142, 123, 3595, 3576],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'CHL': [{'file': 'flag_CHL_original.png',
                            'crop_pixels': [142, 123, 3596, 3577],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'CAN': [{'file': 'flag_CAN_original.png',
                            'crop_pixels': [140, 123, 3595, 3578],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'USA': [{'file': 'flag_USA_original.png',
                            'crop_pixels': [140, 123, 3595, 3578],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'MEX': [{'file': 'flag_MEX_original.png',
                            'crop_pixels': [142, 123, 3596, 3577],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'CHN': [{'file': 'flag_CHN_original.png',
                            'crop_pixels': [140, 124, 3595, 3579],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'CHE': [{'file': 'flag_CHE_original.png',
                            'crop_pixels': [142, 124, 3595, 3577],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'PHL': [{'file': 'flag_PHL_original.png',
                            'crop_pixels': [140, 123, 3595, 3578],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'COD': [{'file': 'flag_COD_original.png',
                            'crop_pixels': [795, 92, 3380, 2677],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'GBR': [{'file': 'flag_GBR_original.png',
                            'crop_pixels': [386, 116, 3639, 3369],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'ECU': [{'file': 'flag_ECU_original.png',
                            'crop_pixels': [142, 123, 3595, 3576],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'BRA': [{'file': 'flag_BRA_original.png',
                            'crop_pixels': [140, 123, 3595, 3578],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}],
                   'ZAF': [{'file': 'flag_ZAF_original.png',
                            'crop_pixels': [142, 124, 3595, 3577],
                            'mask': 'circle',
                            'box': [0, 0, 1, 1]}]},
 'font_families': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
 'country_font_pt': 6.5,
 'value_font_pt': 6.2,
 'minor_font_pt': 5.5,
 'continent_font_pt': 5.5,
 'node_x': [0.06, 0.3, 0.66, 0.88],
 'node_width': 0.006,
 'outer_link_alpha': 0.1,
 'country_label_aliases': {'Democratic Republic of the Congo': 'DR Congo',
                           'Dominican Republic': 'Dominican Rep.',
                           'United Arab Emirates': 'UAE',
                           'Papua New Guinea': 'Papua New Guinea'},
 'scientific_label_decimals': 2,
 'left_country_styles': {'host': {'MEX': {'country_name': 'Mexico',
                                          'node_base_color': '#4076B6',
                                          'node_opacity': 1.0,
                                          'ribbon_color': '#4076B6',
                                          'ribbon_opacity': 0.7},
                                  'USA': {'country_name': 'United States',
                                          'node_base_color': '#4076B6',
                                          'node_opacity': 0.8,
                                          'ribbon_color': '#6691C5',
                                          'ribbon_opacity': 0.7},
                                  'CAN': {'country_name': 'Canada',
                                          'node_base_color': '#4076B6',
                                          'node_opacity': 0.6,
                                          'ribbon_color': '#8CADD3',
                                          'ribbon_opacity': 0.7},
                                  'CHL': {'country_name': 'Chile',
                                          'node_base_color': '#20683A',
                                          'node_opacity': 1.0,
                                          'ribbon_color': '#20683A',
                                          'ribbon_opacity': 0.7},
                                  'AUS': {'country_name': 'Australia',
                                          'node_base_color': '#1D7978',
                                          'node_opacity': 1.0,
                                          'ribbon_color': '#1D7978',
                                          'ribbon_opacity': 0.7},
                                  'PER': {'country_name': 'Peru',
                                          'node_base_color': '#20683A',
                                          'node_opacity': 0.8,
                                          'ribbon_color': '#4D8661',
                                          'ribbon_opacity': 0.7},
                                  'ECU': {'country_name': 'Ecuador',
                                          'node_base_color': '#20683A',
                                          'node_opacity': 0.6,
                                          'ribbon_color': '#79A489',
                                          'ribbon_opacity': 1.0}},
                         'owner': {'CAN': {'country_name': 'Canada',
                                           'node_base_color': '#4076B6',
                                           'node_opacity': 1.0,
                                           'ribbon_color': '#4076B6',
                                           'ribbon_opacity': 0.7},
                                   'CHN': {'country_name': 'China',
                                           'node_base_color': '#A61D31',
                                           'node_opacity': 1.0,
                                           'ribbon_color': '#A61D31',
                                           'ribbon_opacity': 0.7},
                                   'GBR': {'country_name': 'United Kingdom',
                                           'node_base_color': '#F5A238',
                                           'node_opacity': 1.0,
                                           'ribbon_color': '#F5A238',
                                           'ribbon_opacity': 0.7},
                                   'USA': {'country_name': 'United States',
                                           'node_base_color': '#4076B6',
                                           'node_opacity': 0.8,
                                           'ribbon_color': '#6691C5',
                                           'ribbon_opacity': 0.7},
                                   'AUS': {'country_name': 'Australia',
                                           'node_base_color': '#1D7978',
                                           'node_opacity': 1.0,
                                           'ribbon_color': '#1D7978',
                                           'ribbon_opacity': 0.7}}},
 'right_country_colors': {'host': '#4B4B4B', 'owner': '#000000'},
 'default_country_node_opacity': 0.6,
 'default_central_link_opacity': 0.7,
 'country_tone_rule': 'Reuse extracted country styles where available; newly displayed countries use the '
                      'same continent hue at 60% node opacity and a white-tinted ribbon at 70% opacity.'}


def country_appearance(style, perspective, iso3, continent):
    """Pure color lookup; retain reference tones without reusing old values."""
    base = style["continent_colors"][continent]
    reference = style["left_country_styles"][perspective].get(iso3)
    if reference is not None:
        if reference["node_base_color"].upper() != base.upper():
            raise ValueError(f"Reference color conflicts with current continent: {iso3}")
        return {key: reference[key] for key in
                ("node_base_color", "node_opacity", "ribbon_color", "ribbon_opacity")}
    opacity = style["default_country_node_opacity"]
    channels = [int(base[i:i + 2], 16) for i in (1, 3, 5)]
    tinted = "#" + "".join(f"{round(value * opacity + 255 * (1 - opacity)):02X}" for value in channels)
    return {"node_base_color": base, "node_opacity": opacity, "ribbon_color": tinted,
            "ribbon_opacity": style["default_central_link_opacity"]}


def format_area(value):
    # Retain a visible positive value even when it rounds below one km².
    return f"{value:,.0f}" if value >= 1 else f"{value:.2g}"


def scientific_area(value, decimals=2):
    exponent = floor(log10(value))
    mantissa = round(value / 10 ** exponent, decimals)
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    superscript = str(exponent).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    # Unicode preserves the configured font size instead of shrinking mathtext.
    return f"{mantissa:.{decimals}f}×10{superscript}"


@lru_cache(maxsize=24)
def _icon(code):
    """Apply recorded original-asset crops only when drawing; preserve circles."""
    import numpy as np
    from PIL import Image, ImageChops, ImageDraw

    style = STYLE
    if code in style["continent_icons"]:
        layers = [{"file": style["continent_icons"][code], "box": [0, 0, 1, 1], "crop": {}}]
    else:
        layers = style["country_icons"][code]
    canvas = Image.new("RGBA", (256, 256), (255, 255, 255, 0))
    for layer in layers:
        if ASSET_DIR is None:
            raise ValueError("Supply --data-root or set COPPER_DATA_ROOT before rendering icons")
        path = ASSET_DIR / layer["file"]
        if not path.is_file():
            raise FileNotFoundError(f"Required original icon is missing: {path}")
        with Image.open(path) as source:
            rgba = source.convert("RGBA")
        width, height = rgba.size
        if "crop_pixels" in layer:
            bounds = tuple(layer["crop_pixels"])
            if (len(bounds) != 4 or not all(isinstance(v, int) for v in bounds)
                    or not 0 <= bounds[0] < bounds[2] <= width
                    or not 0 <= bounds[1] < bounds[3] <= height):
                raise ValueError(f"Invalid pixel crop for {path.name}")
        else:
            crop = layer.get("crop", {})
            bounds = (round(crop.get("l", 0) * width), round(crop.get("t", 0) * height),
                      round((1 - crop.get("r", 0)) * width), round((1 - crop.get("b", 0)) * height))
        x, y, w, h = layer["box"]
        cropped = rgba.crop(bounds).resize((round(w * 256), round(h * 256)), Image.Resampling.LANCZOS)
        if layer.get("mask") == "circle":
            if bounds[2] - bounds[0] != bounds[3] - bounds[1] or cropped.width != cropped.height:
                raise ValueError(f"Circular icon must have a square crop and display box: {path.name}")
            # Supersample transparency, retaining original colors and internal bevel.
            diameter = cropped.width * 4
            mask = Image.new("L", (diameter, diameter), 0)
            ImageDraw.Draw(mask).ellipse((0, 0, diameter - 1, diameter - 1), fill=255)
            mask = mask.resize(cropped.size, Image.Resampling.LANCZOS)
            cropped.putalpha(ImageChops.multiply(cropped.getchannel("A"), mask))
        elif layer.get("mask") is not None:
            raise ValueError(f"Unknown icon mask for {path.name}")
        canvas.alpha_composite(cropped, (round(x * 256), round(y * 256)))
    if code in style["continent_icons"]:
        # Fill only enclosed transparent regions. Padding lets the flood fill
        # reach all exterior corners even when a circular rim touches the edge.
        alpha = canvas.getchannel("A")
        barrier = alpha.point(lambda value: 255 if value >= 128 else 0)
        padded = Image.new("L", (258, 258), 0)
        padded.paste(barrier, (1, 1))
        ImageDraw.floodfill(padded, (0, 0), 128)
        interior = padded.crop((1, 1, 257, 257)).point(lambda value: 255 if value == 0 else 0)
        white = Image.new("RGBA", canvas.size, "white")
        white.putalpha(interior)
        canvas = Image.alpha_composite(white, canvas)
    return np.asarray(canvas)


def _ribbon(ax, x0, x1, source, target, color, alpha):
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch

    a, b = source
    c, d = target
    bend = (x1 - x0) * 0.5
    vertices = [(x0, a), (x0 + bend, a), (x1 - bend, c), (x1, c),
                (x1, d), (x1 - bend, d), (x0 + bend, b), (x0, b), (x0, a)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, zorder=1, rasterized=False))


def draw_alluvial(ax, bilateral, perspective):
    from matplotlib import font_manager
    from matplotlib.patches import Rectangle

    style = STYLE
    model = build_model(bilateral.to_dict("records"), perspective)
    width_pt = ax.get_position().width * ax.figure.get_figwidth() * 72
    height_pt = ax.get_position().height * ax.figure.get_figheight() * 72
    layout = allocate_layout(model, height_pt, y_min=0.055)
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in style["font_families"] if name in available_fonts), None)
    if family is None:
        raise RuntimeError("None of the configured Sankey fonts is available")
    major = major_keys(model)
    country_styles = {node.key: country_appearance(style, perspective, node.key[1], node.continent)
                      for node in model.columns[1]}
    nodes = {n.key: n for column in model.columns for n in column}
    xs, nw = list(style["node_x"]), style["node_width"]
    # The left edge of each 11 pt continent icon meets the map-frame edge.
    xs[0] = 5.5 / width_pt
    ax.set(xlim=(0, 1), ylim=(1, 0), aspect="auto")
    ax.axis("off")

    def text(x, y, value, size=5.5, **kwargs):
        return ax.text(x, y, value, fontsize=size, fontfamily=family, va="center",
                       color="#111111", zorder=10, clip_on=False, **kwargs)

    def icon(key, x, y, diameter_pt):
        dx, dy = diameter_pt / width_pt / 2, diameter_pt / height_pt / 2
        ax.imshow(_icon(key), extent=(x - dx, x + dx, y + dy, y - dy), origin="upper",
                  aspect="auto", interpolation="antialiased", zorder=8)

    for link in model.links:
        start, end = link.source[0], link.target[0]
        color_node = nodes[link.source] if start < 2 else nodes[link.target]
        color = style["continent_colors"][color_node.continent]
        alpha = style["outer_link_alpha"]
        if start == 1:
            appearance = country_styles[link.source]
            color, alpha = appearance["ribbon_color"], appearance["ribbon_opacity"]
        key = (link.source, link.target)
        _ribbon(ax, xs[start] + nw / 2, xs[end] - nw / 2,
                layout["source_bands"][key], layout["target_bands"][key], color, alpha)

    missing_flags = []
    for column in model.columns:
        for node in column:
            index, code = node.key
            top, bottom = layout["node_bounds"][node.key]
            cy = (top + bottom) / 2
            color = style["right_country_colors"][perspective] if index == 2 else style["continent_colors"][node.continent]
            opacity = country_styles[node.key]["node_opacity"] if index == 1 else 1.0
            ax.add_patch(Rectangle((xs[index] - nw / 2, top), nw, bottom - top,
                                   facecolor=color, alpha=opacity, edgecolor="none", zorder=4))
            if index in (0, 3):
                # Place continent icons only on the left.
                if index == 0:
                    icon_y = min(cy, top + 5.5 / height_pt)
                    icon(node.name, xs[0], icon_y, 11.0)
                    continue
                text(0.945, cy - 4.0 / height_pt, node.name.replace(" ", "\n"),
                     style["continent_font_pt"], ha="center", linespacing=1.0)
                text(0.945, cy + 6.0 / height_pt, scientific_area(node.value, style["scientific_label_decimals"]),
                     style["continent_font_pt"], ha="center", fontweight="bold")
                continue
            label = style["country_label_aliases"].get(node.name, node.name)
            if node.key in major:
                label_x, align = (xs[1] - 0.014, "right") if index == 1 else (xs[2] + 0.014, "left")
                text(label_x, cy - 4.5 / height_pt, label, style["country_font_pt"], ha=align)
                value_x = label_x
                text(value_x, cy + 4.5 / height_pt, format_area(node.value),
                     style["value_font_pt"], ha=align, fontweight="bold")
                if index == 1:
                    if code in style["country_icons"]:
                        icon(code, xs[1] - 0.105, cy + 4.5 / height_pt, 9.0)
                    else:
                        missing_flags.append(code)
            else:
                text(xs[2] + 0.054, cy, format_area(node.value), style["minor_font_pt"],
                     ha="right", fontweight="bold")
                text(xs[2] + 0.064, cy, label, style["minor_font_pt"], ha="left")

    headings = (["", "Owner-registration country", "Mine-host country", ""]
                if perspective == "owner" else
                ["", "Mine-host country", "Owner-registration countries and regions", ""])
    for x, heading in zip([xs[0], xs[1], xs[2], 0.945], headings):
        if heading:
            text(x, 0.020, heading, 6.2, ha="center", fontweight="bold")
    text(0.5, 0.989, "Equity-attributed affected habitat area (km²)", 5.5, ha="center")
    ax.set(xlim=(0, 1), ylim=(1, 0), aspect="auto")
    return {**model.summary(), "missing_flag_iso3": sorted(set(missing_flags))}
