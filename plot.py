import os
import logging
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from typing import Optional
from statistics import median
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Folder to store images
TEMP_DIR = 'temp'

# Ensure the temp directory exists
if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)
    logger.info(f"Created temporary directory: {TEMP_DIR}")
else:
    logger.debug(f"Temporary directory already exists: {TEMP_DIR}")

class Plot:
    def __init__(self, title: str, x_axis: str, y_axis: str, lo: Optional[float] = None, hi: Optional[float] = None):
        self.title = title
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.lo = lo
        self.hi = hi

    def _add_break(self, y_vals: pd.Series, threshold_ratio: float = 5.0) -> bool:
        if y_vals.empty:
            logger.warning("Y values are empty. Cannot determine if break is needed.")
            return False
        non_zero_vals = y_vals[y_vals > 0]
        if len(non_zero_vals) < 2:
            logger.info("Not enough non-zero Y values to assess break.")
            return False
        return non_zero_vals.max() / non_zero_vals.min() > threshold_ratio

    def _plot_break(self, ax_top, ax_bottom):
        kwargs = dict(marker=[(-1, -1), (1, 1)], markersize=6,
                      linestyle='none', color='k', mec='k', mew=1, clip_on=False)
        ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
        ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)

    def _finalize_plot(self, fig, plot_filename):
        path = os.path.join(TEMP_DIR, plot_filename)
        plot_path = os.path.abspath(path)
        plt.savefig(plot_path, format='png')
        plt.close(fig)
        return plot_path

    def _set_plain_yaxis_format(self, ax):
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter)
        ax.ticklabel_format(style='plain', axis='y')

    def generate_plots(self, data_frame, plot_filename: str):
        data_frame = data_frame.sort_values(by=self.x_axis)
        y_vals = data_frame[self.y_axis].dropna()

        if y_vals.empty:
            logger.warning("No data available for Y-axis. Skipping plot generation.")
            return None

        use_break = self._add_break(y_vals)

        if use_break:
            fig, (ax_top, ax_bottom) = plt.subplots(
                2, 1, sharex=True, figsize=(12, 6), dpi=100,
                gridspec_kw={'height_ratios': [6, 3], 'hspace': 0.001},
                constrained_layout=True
            )
            axes = (ax_top, ax_bottom)

            y_min = y_vals[y_vals > 0].min()
            y_max = y_vals.max()
            ax_bottom.set_ylim(0, y_min * 1.1)
            ax_top.set_ylim(y_min * 1.5, y_max * 1.1)

            self._plot_break(ax_top, ax_bottom)
            ax_for_labels = ax_top
        else:
            fig, ax = plt.subplots(figsize=(12, 4), dpi=100, constrained_layout=True)
            axes = (ax,)
            ax_for_labels = ax
            y_max = y_vals.max()
            ax.set_ylim(0, self.hi if self.hi is not None else y_max * 1.1)

        for ax in axes:
            ax.bar(data_frame[self.x_axis], data_frame[self.y_axis], color="#001675", edgecolor='w', alpha=0.6)
            ax.set_xticks([])
            ax.set_xticklabels([])
            self._set_plain_yaxis_format(ax)

        y_med = median(y_vals)
        ax_for_labels.axhline(y=y_med, color='#FF5733', linestyle='--', linewidth=1, label=f"Median: {y_med:.2f}")
        ax_for_labels.legend(fontsize="small")
        ax_for_labels.set_ylabel(self.y_axis, fontsize=12)

        return self._finalize_plot(fig, plot_filename)

    def generate_Seqplots(self, data_frame, plot_filename: str):
        data_frame = data_frame.sort_values(by=self.x_axis)
        y_vals = data_frame[self.y_axis].dropna()

        if y_vals.empty:
            logger.warning("No data available for Y-axis. Skipping SeqPlot generation.")
            return None

        use_break = self._add_break(y_vals)

        sample_map = {
            'Matched Normal': '#162456',
            'Tumour': '#E7180B',
        }

        marker_map = {
            'Matched Normal': 'o',
            'Tumour': '^',
        }

        use_sample_alignment = 'SampleID' in data_frame.columns and 'Sample Type' in data_frame.columns

        if use_sample_alignment:
            sample_ids = sorted(data_frame['SampleID'].unique())
            x_positions = {}
            pos = 0
            for sample_id in sample_ids:
                x_positions[(sample_id, 'Matched Normal')] = pos
                x_positions[(sample_id, 'Tumour')] = pos + 0.4
                pos += 1.0
            grouped = data_frame.groupby(['SampleID', 'Sample Type'])

        if use_break:
            fig, (ax_top, ax_bottom) = plt.subplots(
                2, 1, sharex=True, figsize=(12, 6), dpi=100,
                gridspec_kw={'height_ratios': [6, 3], 'hspace': 0.001},
                constrained_layout=True
            )
            axes = (ax_top, ax_bottom)

            y_min = y_vals[y_vals > 0].min()
            y_max = y_vals.max()
            ax_bottom.set_ylim(0, y_min * 1.1)
            ax_top.set_ylim(y_min * 1.5, y_max * 1.1)

            self._plot_break(ax_top, ax_bottom)
            ax_for_labels = ax_top
        else:
            fig, ax = plt.subplots(figsize=(12, 4), dpi=100, constrained_layout=True)
            axes = (ax,)
            ax_for_labels = ax
            y_max = y_vals.max()
            ax.set_ylim(0, self.hi if self.hi is not None else y_max * 1.1)

        for ax in axes:
            if use_sample_alignment:
                for (sample_id, sample_type), group in grouped:
                    if (sample_id, sample_type) not in x_positions:
                        continue
                    x_val = x_positions[(sample_id, sample_type)]
                    ax.scatter(
                        [x_val] * len(group),
                        group[self.y_axis],
                        color=sample_map.get(sample_type, '#E7180B'),
                        marker=marker_map.get(sample_type, 'o'),
                        label=sample_type,
                        edgecolors='w', alpha=0.6, s=60
                    )
            elif 'Sample Type' in data_frame.columns:
                for sample_type, group in data_frame.groupby('Sample Type'):
                    ax.scatter(
                        group[self.x_axis], group[self.y_axis],
                        color=sample_map.get(sample_type, '#E7180B'),
                        marker=marker_map.get(sample_type, 'o'),
                        label=sample_type,
                        edgecolors='w', alpha=0.6, s=60
                    )
            else:
                ax.scatter(
                    data_frame[self.x_axis], data_frame[self.y_axis],
                    color='#0F0094', edgecolors='w', alpha=0.6, s=60
                )

            ax.set_xticks([])
            ax.set_xticklabels([])
            self._set_plain_yaxis_format(ax)

        y_med = median(y_vals)
        ax_for_labels.axhline(y=y_med, color='#FF5733', linestyle='--', linewidth=1, label=f"Median: {y_med:.2f}")
        ax_for_labels.set_ylabel(self.y_axis, fontsize=12)

        handles, labels = ax_for_labels.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax_for_labels.legend(unique.values(), unique.keys(), fontsize=10, loc='best')

        return self._finalize_plot(fig, plot_filename)
