"""Align panel furniture using physical translations, never resize data frames.

The existing builders remain the reference for content size, axes limits and
map geometry. This pass only translates complete axes, positions panel letters
and external legends, and trims peripheral white space.
"""

from __future__ import annotations

import numpy as np
from matplotlib.legend import Legend
from matplotlib.text import Text

OUTER_MM = 2.5
GUTTER_MM = 5.0
LEGEND_GAP_MM = 2.5
LABEL_SLOT_MM = 2.0
LABEL_GAP_MM = 2.0
LABEL_PT = 7.0


def _mm(fig, box):
    return np.array(box.extents, dtype=float) / fig.dpi * 25.4


def _union(boxes):
    values = np.array(boxes, dtype=float)
    return np.array([values[:, 0].min(), values[:, 1].min(),
                     values[:, 2].max(), values[:, 3].max()])


def _frame(fig, ax):
    return _mm(fig, ax.get_window_extent())


def _content(fig, group):
    renderer = fig.canvas.get_renderer()
    boxes = [_frame(fig, ax) for ax in group]
    for ax in group:
        box = ax.get_tightbbox(renderer)
        if box is not None:
            boxes.append(_mm(fig, box))
    return _union(boxes)


def _groups(fig, number):
    axes = fig.axes
    indices = {
        1: [[0], [1], [2, 3]], 2: [[0], [1]],
        3: [[0], [1], [2], [3], [4], [5]],
        4: [[0], [1], [2], [3], [4]], 5: [[0], [1]],
        6: [[0], [1], [2]], 7: [[0], [1, 2]],
    }[number]
    expected_count = sum(map(len, indices))
    if len(axes) != expected_count:
        raise ValueError(f'Figure {number}: axes structure changed ({len(axes)})')
    return [[axes[i] for i in members] for members in indices]


def _put_frames(fig, frames, offsets):
    width, height = fig.get_size_inches() * 25.4
    for ax, original in frames.items():
        dx, dy = offsets[ax]
        x0, y0, x1, y1 = original
        ax.set_position([(x0 + dx) / width, (y0 + dy) / height,
                         (x1 - x0) / width, (y1 - y0) / height])


def _external_legends(fig, groups, number):
    """Use rendered content edges, including tick labels, for legend clearance."""
    shared = groups[1][0].get_legend() if number == 3 else None
    if shared is not None:
        shared.set_in_layout(False)
    movable = []
    for group in groups:
        for ax in group:
            frame = _frame(fig, ax)
            legends = [item for item in ax.get_children()
                       if isinstance(item, Legend) and item is not shared]
            external = []
            for legend in legends:
                box = _mm(fig, legend.get_window_extent(fig.canvas.get_renderer()))
                if box[1] >= frame[3] - 0.1:
                    external.append((legend, box, 'above'))
                elif box[3] <= frame[1] + 0.1:
                    external.append((legend, box, 'below'))
            for legend, _, _ in external:
                legend.set_in_layout(False)
            content = _content(fig, [ax])
            for legend, box, side in external:
                dy = (content[3] + LEGEND_GAP_MM - box[1] if side == 'above'
                      else content[1] - LEGEND_GAP_MM - box[3])
                # Keep the anchor in axes coordinates so it follows translation.
                anchor = legend.get_bbox_to_anchor().transformed(ax.transAxes.inverted())
                coords = list(anchor.bounds)
                coords[1] += dy / (frame[3] - frame[1])
                legend.set_bbox_to_anchor(coords, transform=ax.transAxes)
                legend.set_in_layout(True)
                movable.append({'side': side, 'gap_mm': LEGEND_GAP_MM})
    fig.canvas.draw()
    return shared, movable


def align_figure(fig, number):
    """Apply the agreed spacing convention and verify all content dimensions."""
    if getattr(fig, '_spacing_alignment_done', False):
        raise ValueError('The spacing pass must run once on a fresh figure')
    fig.canvas.draw()
    groups = _groups(fig, number)
    before_canvas = list(fig.get_size_inches() * 25.4)
    for item in list(fig.findobj(match=Text)):
        if item.get_gid() == 'panel-label':
            item.remove()
    authorized_resize = None
    if number == 1:
        # The author explicitly requested a narrower country chart, including
        # both y-axis titles, inside the map's horizontal boundaries.
        fig.canvas.draw()
        map_frame = _frame(fig, groups[0][0])
        old_frame = _frame(fig, groups[2][0])
        content = _content(fig, groups[2])
        left = map_frame[0] + old_frame[0] - content[0]
        right = map_frame[2] - (content[2] - old_frame[2])
        width, height = fig.get_size_inches() * 25.4
        for ax in groups[2]:
            ax.set_position([left / width, old_frame[1] / height,
                             (right - left) / width, (old_frame[3] - old_frame[1]) / height])
        fig.canvas.draw()
        fitted = _content(fig, groups[2])
        if fitted[0] < map_frame[0] - 0.1 or fitted[2] > map_frame[2] + 0.1:
            raise ValueError('Figure 1c: full chart does not fit within the map boundaries')
        authorized_resize = {'panel': 'c', 'before_frame_mm': old_frame.tolist(),
                             'after_frame_mm': _frame(fig, groups[2][0]).tolist(),
                             'reason': 'author requested horizontal contraction including both axis titles'}
    all_axes = [*fig.axes, *(child for ax in fig.axes for child in ax.child_axes)]
    before = {ax: {'frame': _frame(fig, ax), 'xlim': tuple(ax.get_xlim()),
                   'ylim': tuple(ax.get_ylim())} for ax in all_axes}
    frames = {ax: before[ax]['frame'].copy() for ax in fig.axes}
    # Colourbar locators can otherwise re-impose an old GridSpec position.
    # Insets keep their parent-relative locators and therefore move with parents.
    for ax in fig.axes:
        ax.set_axes_locator(None)
    for item in list(fig.findobj(match=Text)):
        if item.get_gid() == 'panel-label':
            item.remove()
    fig.canvas.draw()
    shared, legend_moves = _external_legends(fig, groups, number)
    boxes = [_content(fig, group) for group in groups]
    shift = np.zeros((len(groups), 2))
    reserve = LABEL_SLOT_MM + LABEL_GAP_MM
    gutter = 3.0 if number == 6 else GUTTER_MM

    if number in (1, 6):
        # Full-width vertical stacks: preserve common data-frame x coordinates.
        for i in range(1, len(groups)):
            shift[i, 1] = boxes[i - 1][1] + shift[i - 1, 1] - gutter - boxes[i][3]
    elif number in (2, 5, 7):
        shift[1, 0] = boxes[0][2] + GUTTER_MM + reserve - boxes[1][0]
    elif number == 3:
        # Shared column translations retain the six identical data-frame sizes.
        for col in (1, 2):
            previous = max(boxes[i][2] + shift[i, 0] for i in (col - 1, col + 2))
            left = min(boxes[i][0] for i in (col, col + 3))
            shift[col, 0] = shift[col + 3, 0] = previous + GUTTER_MM + reserve - left
        dy = min(box[1] for box in boxes[:3]) - GUTTER_MM - max(box[3] for box in boxes[3:])
        shift[3:, 1] = dy
    # Figure 4 keeps every relative map-frame position exactly unchanged.
    offsets = {ax: shift[i].copy() for i, group in enumerate(groups) for ax in group}
    _put_frames(fig, frames, offsets)
    fig.canvas.draw()
    boxes = [_content(fig, group) for group in groups]

    labels = []
    if groups:
        if number in (1, 6):
            xs = [min(box[0] for box in boxes) - reserve] * len(groups)
            ys = [box[3] for box in boxes]
        elif number in (2, 5, 7):
            xs = [box[0] - reserve for box in boxes]
            ys = [max(box[3] for box in boxes)] * len(groups)
        elif number == 3:
            col_x = [min(boxes[i][0], boxes[i + 3][0]) - reserve for i in range(3)]
            xs = col_x * 2
            ys = [max(box[3] for box in boxes[:3])] * 3 + [max(box[3] for box in boxes[3:])] * 3
        else:
            left = min(boxes[0][0], boxes[2][0]) - reserve
            # Retain the original compact map block. Place regional letters
            # above the latitude ticks, clear of the preceding map's longitude
            # labels, rather than reserving a new gutter between map frames.
            right = _frame(fig, groups[1][0])[0] - 6.0
            middle = _frame(fig, groups[3][0])[0] - 6.0
            xs = [left, right, left, middle, right]
            ys = [max(box[3] for box in boxes[:2]) + 1.0] * 2 + [max(box[3] for box in boxes[2:]) + 1.0] * 3
        width, height = fig.get_size_inches() * 25.4
        for i, (x, y) in enumerate(zip(xs, ys)):
            labels.append(fig.text(x / width, y / height, chr(97 + i),
                                   fontsize=LABEL_PT, fontweight='bold',
                                   ha='left', va='top', gid='panel-label'))

    shared_anchor = None
    if shared is not None:
        # A single legend sits below the complete six-panel composition.
        left = min(_frame(fig, group[0])[0] for group in groups[:3])
        right = max(_frame(fig, group[0])[2] for group in groups[:3])
        x = (left + right) / 2
        y = min(box[1] for box in boxes[3:]) - LEGEND_GAP_MM
        width, height = fig.get_size_inches() * 25.4
        shared.set_bbox_to_anchor((x / width, y / height), transform=fig.transFigure)
        fig.canvas.draw()
        actual = _mm(fig, shared.get_window_extent(fig.canvas.get_renderer()))
        y += min(box[1] for box in boxes[3:]) - LEGEND_GAP_MM - actual[3]
        shared_anchor = np.array([x, y])
        shared.set_bbox_to_anchor((x / width, y / height), transform=fig.transFigure)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panel_boxes = [_content(fig, group) for group in groups]
    for i, label in enumerate(labels):
        panel_boxes[i] = _union([panel_boxes[i], _mm(fig, label.get_window_extent(renderer))])
    complete = panel_boxes.copy()
    if shared is not None:
        complete.append(_mm(fig, shared.get_window_extent(renderer)))
    bounds = _union(complete)
    translation = np.array([OUTER_MM - bounds[0], OUTER_MM - bounds[1]])
    size = bounds[2:] - bounds[:2] + 2 * OUTER_MM
    label_positions = [np.array(item.get_position()) * fig.get_size_inches() * 25.4 for item in labels]
    fig.set_size_inches(*(size / 25.4), forward=False)
    for ax in offsets:
        offsets[ax] += translation
    _put_frames(fig, frames, offsets)
    for item, position in zip(labels, label_positions):
        item.set_position((position + translation) / size)
    if shared is not None:
        shared.set_bbox_to_anchor(tuple((shared_anchor + translation) / size), transform=fig.transFigure)
    fig.canvas.draw()

    checks = []
    for index, ax in enumerate(all_axes):
        old = before[ax]
        after = _frame(fig, ax)
        old_size = old['frame'][2:] - old['frame'][:2]
        new_size = after[2:] - after[:2]
        if not np.allclose(old_size, new_size, rtol=0, atol=1e-6):
            raise ValueError(f'Figure {number}, axes {index}: content size changed: {old_size} -> {new_size}')
        if not np.allclose(ax.get_xlim(), old['xlim']) or not np.allclose(ax.get_ylim(), old['ylim']):
            raise ValueError(f'Figure {number}, axes {index}: data limits changed')
        checks.append({'axes_index': index, 'before_frame_mm': old['frame'].tolist(),
                       'after_frame_mm': after.tolist(), 'width_height_preserved': True,
                       'data_limits_preserved': True})
    if number == 4:
        displacements = np.array([item['after_frame_mm'][:2] for item in checks]) - np.array([item['before_frame_mm'][:2] for item in checks])
        if not np.allclose(displacements, displacements[0], rtol=0, atol=1e-6):
            raise ValueError('Figure 4: relative map-frame arrangement changed')
    renderer = fig.canvas.get_renderer()
    final_boxes = [_content(fig, group) for group in groups]
    for i, label in enumerate(labels):
        final_boxes[i] = _union([final_boxes[i], _mm(fig, label.get_window_extent(renderer))])
    all_boxes = final_boxes.copy()
    if shared is not None:
        all_boxes.append(_mm(fig, shared.get_window_extent(renderer)))
    final_bounds = _union(all_boxes)
    margins = [final_bounds[0], final_bounds[1], size[0] - final_bounds[2], size[1] - final_bounds[3]]
    if not np.allclose(margins, OUTER_MM, rtol=0, atol=0.08):
        raise ValueError(f'Figure {number}: outer margins changed after canvas trim: {margins}')
    fig._spacing_alignment_done = True
    return {'policy': 'physical translation only; no axes resizing or data changes',
            'before_canvas_mm': before_canvas, 'after_canvas_mm': size.tolist(),
            'author_requested_chart_resize': authorized_resize,
            'outer_margins_left_bottom_right_top_mm': margins,
            'panel_label_font_pt': LABEL_PT,
            'target_label_content_gap_mm': LABEL_GAP_MM,
            'target_panel_gutter_mm': gutter, 'external_legend_adjustments': legend_moves,
            'panel_content_bounds_mm': [box.tolist() for box in final_boxes],
            'content_geometry_checks': checks,
            'figure_4_relative_map_positions_preserved': number == 4,
            'figure_4_gutter_exception': 'original map block retained exactly' if number == 4 else None}
