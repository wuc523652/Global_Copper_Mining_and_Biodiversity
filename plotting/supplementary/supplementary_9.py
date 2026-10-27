"""Analytical workflow drawn as editable vector boxes, arrows, and text."""
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch
from matplotlib.path import Path as MplPath
WIDTH_MM, HEIGHT_MM = 183, 163
MM_TO_PT = 72 / 25.4

COLORS = {
    'input': '#FFFFFF', 'neutral': '#F4F5F6',
    'regional': '#E5EFE7', 'regional_out': '#BBD8C1',
    'spatial': '#E9F0F7', 'species': '#BDD3E8',
    'ownership': '#EEEAF3', 'ownership_out': '#C9BDDF',
    'portfolio': '#CDD7DF', 'edge': '#303438', 'text': '#15191C',
}

NODES = {'en': {'country': (6, 143, 39, 15, 'Country boundaries', 'GADM', 'input', True),
        'habitat': (48,
                    143,
                    39,
                    15,
                    'Species’ global\nhabitats',
                    '5,890 vertebrates',
                    'input',
                    True),
        'mine': (90, 143, 39, 15, 'Copper-mine\nboundaries', '1,278 mines', 'input', True),
        'records': (132,
                    143,
                    39,
                    15,
                    'Mine-project ownership\ninformation',
                    'Equity ownership structure',
                    'input',
                    True),
        'integration': (6,
                        116,
                        111,
                        16,
                        'Spatial data integration',
                        'Mollweide equal-area projection',
                        'neutral',
                        False),
        'rwr': (6,
                92,
                51,
                16,
                'Rarity-weighted richness',
                '10 × 10 km equal-area grid',
                'regional',
                False),
        'buffer': (66, 92, 51, 16, 'Potential influence zones', '30-km buffers', 'spatial', False),
        'regions': (6,
                    68,
                    51,
                    16,
                    'Regions of high\nconservation value',
                    'Global top 10%',
                    'regional_out',
                    False),
        'overlap': (66,
                    68,
                    51,
                    16,
                    'Spatial overlap analysis',
                    'Affected habitat area\nand proportion',
                    'spatial',
                    False),
        'holder': (126,
                   68,
                   51,
                   16,
                   'Attribution to project\nequity holders',
                   'Reported project equity shares',
                   'ownership',
                   False),
        'randomization': (6,
                          44,
                          51,
                          16,
                          'Country-stratified\nrandomization',
                          'Observed vs. expected\nmine occurrence',
                          'regional',
                          False),
        'priority': (66,
                     44,
                     51,
                     16,
                     'Conservation-priority\nspecies',
                     'Habitat-area-dependent thresholds',
                     'species',
                     False),
        'owner_country': (126,
                          44,
                          51,
                          16,
                          'Aggregation by\nowner-registration country',
                          '',
                          'ownership',
                          False),
        'portfolio': (66,
                      19,
                      51,
                      16,
                      'Project equity portfolios\nby owner-registration country',
                      'Regional and species-level criteria',
                      'portfolio',
                      False),
        'linkages': (126,
                     19,
                     51,
                     16,
                     'Domestic and cross-border\nresponsibility linkages',
                     '',
                     'ownership_out',
                     False)}}

def build_plot():
    lang = "en"
    font = 'Arial' if lang == 'en' else 'Microsoft YaHei'
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, WIDTH_MM), ylim=(0, HEIGHT_MM))
    ax.set_axis_off()
    labels = []
    ax.add_patch(Rectangle((2, 16), 179, 145, facecolor='none',
                           edgecolor='#6B7175', linewidth=0.8,
                           linestyle=(0, (5, 3)), zorder=0))

    def label(text, x, y, size, weight='normal', ha='center'):
        item = ax.text(x, y, text, ha=ha, va='center', fontsize=size,
                       fontfamily=font, fontweight=weight, color=COLORS['text'], zorder=4)
        return item

    for key, data in NODES[lang].items():
        x, y, w, h, title, detail, family, parallelogram = data
        if parallelogram:
            skew = 3
            patch = Polygon([(x, y), (x + w - skew, y), (x + w, y + h),
                             (x + skew, y + h)], closed=True)
        else:
            patch = Rectangle((x, y), w, h)
        patch.set(facecolor=COLORS[family], edgecolor=COLORS['edge'], linewidth=0.85, zorder=3)
        ax.add_patch(patch)
        title_size = 8.5
        detail_size = 7.0
        lines = [(s, title_size, 'bold' if key in ('regions', 'priority', 'portfolio', 'linkages') else 'normal')
                 for s in title.split('\n')]
        if detail:
            lines += [(s, detail_size, 'normal') for s in detail.split('\n')]
        heights = [size / MM_TO_PT * 1.2 for _, size, _ in lines]
        total = sum(heights) + (0.7 if detail else 0)
        pos = y + (h + total) / 2
        for i, (s, size, weight) in enumerate(lines):
            if i == len(title.split('\n')) and detail:
                pos -= 0.7
            pos -= heights[i] / 2
            item = label(s, x + w / 2, pos, size, weight)
            labels.append((key, item, (x + (4 if parallelogram else 2), y + 1,
                                        x + w - (4 if parallelogram else 2), y + h - 1)))
            pos -= heights[i] / 2

    def arrow(points):
        path = MplPath(points, [MplPath.MOVETO] + [MplPath.LINETO] * (len(points) - 1))
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>', mutation_scale=7,
                                    linewidth=0.85, color=COLORS['edge'],
                                    capstyle='butt', joinstyle='miter', zorder=2))

    # Three separate inputs enter the wide integration box vertically.
    for x in (25.5, 67.5, 109.5):
        arrow([(x, 143), (x, 132)])
    arrow([(31.5, 116), (31.5, 108)])
    arrow([(91.5, 116), (91.5, 108)])
    arrow([(31.5, 92), (31.5, 84)])
    arrow([(91.5, 92), (91.5, 84)])
    arrow([(31.5, 68), (31.5, 60)])
    arrow([(91.5, 68), (91.5, 60)])
    arrow([(117, 76), (126, 76)])
    arrow([(151.5, 143), (151.5, 84)])
    arrow([(151.5, 68), (151.5, 60)])
    arrow([(151.5, 44), (151.5, 35)])
    arrow([(6, 76), (3.5, 76), (3.5, 27), (66, 27)])
    arrow([(91.5, 44), (91.5, 35)])
    # The portfolio-to-linkages connection summarizes interpretation; area
    # attribution still follows the separate ownership column, as in the caption.
    arrow([(117, 27), (126, 27)])
    label('Area' if lang == 'en' else '面积', 121.5, 79.5, 6.5)
    label('Proportion' if lang == 'en' else '占比', 94, 64, 6.5, ha='left')

    ax.add_patch(Rectangle((6, 3), 171, 10, facecolor='#F5F5F5',
                           edgecolor='#848C92', linewidth=0.8,
                           linestyle=(0, (5, 3)), zorder=0))
    label('Sensitivity analyses' if lang == 'en' else '敏感性分析', 91.5, 10, 8.5, 'bold')
    label('Buffer distances (10–50 km) • Alternative screening rule • Complete ownership records'
          if lang == 'en' else '缓冲距离（10–50 km） • 替代筛选规则 • 完整所有权记录',
          91.5, 5.8, 7.0)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    failures = []
    boxes = []
    for key, item, boundary in labels:
        bb = item.get_window_extent(renderer).transformed(ax.transData.inverted())
        x0, y0, x1, y1 = boundary
        if bb.x0 < x0 - 0.1 or bb.x1 > x1 + 0.1 or bb.y0 < y0 - 0.1 or bb.y1 > y1 + 0.1:
            failures.append({'node': key, 'text': item.get_text(), 'bbox': list(bb.extents), 'boundary': boundary})
        boxes.append({'node': key, 'text': item.get_text(), 'font_pt': item.get_fontsize(),
                      'bbox_mm': list(bb.extents)})
    if failures:
        raise ValueError(json.dumps(failures, ensure_ascii=False))
    return fig
