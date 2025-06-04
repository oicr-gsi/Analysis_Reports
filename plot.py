import os
import matplotlib.pyplot as plt
from typing import Optional
from statistics import median

# Folder to store images
TEMP_DIR = 'temp'

# Ensure the temp directory exists
if not os.path.exists(TEMP_DIR):
    os.makedirs(TEMP_DIR)

class Plot:
    def __init__(self, title: str, x_axis: str, y_axis: str, lo: Optional[float] = None, hi: Optional[float] = None):
        self.title = title
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.lo = lo
        self.hi = hi

    def generate_plots(self, data_frame, plot_filename: str):
        '''
        Generates a scatter plot and saves it as an image file.

        Parameters:
        - data_frame: DataFrame with data for this plot
        - plot_filename: The name of the plot image file to be saved

        Returns:
        - The path to the saved plot image
        '''
        fig, ax = plt.subplots(figsize=(12, 4), dpi=100)
        ax.set_title(self.title, fontsize=14)
        ax.set_xlabel(self.x_axis, fontsize=12)
        ax.set_ylabel(self.y_axis, fontsize=12)

        # Sort by x_axis for consistency
        data_frame = data_frame.sort_values(by=self.x_axis)

        # Create the scatter plot
        ax.scatter(data_frame[self.x_axis], data_frame[self.y_axis], color="#6495ED", edgecolors='w', alpha=0.6)

        # Draw median line
        y_val = data_frame[self.y_axis].dropna().tolist()
        if y_val:
            y_med = median(y_val)
            ax.axhline(y=y_med, color='#FF5733', linestyle='--', linewidth=1, label=f"Median: {y_med:.2f}")
            ax.legend(fontsize="small")

        # Set y-axis range
        if self.hi is not None and self.lo is not None:
            ax.set_ylim(self.lo, self.hi)
        elif self.lo is not None:
            ax.set_ylim(bottom=self.lo)
        elif self.hi is not None:
            ax.set_ylim(top=self.hi)

        # Remove a_axis ticks
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

        # Save the plot as a PNG file in the temp directory
        path = os.path.join(TEMP_DIR, plot_filename)
        plot_path = os.path.abspath(path)
        plt.tight_layout()
        plt.savefig(plot_path, format='png')
        plt.close(fig)

        return plot_path
    
    def generate_Seqplots(self, data_frame, plot_filename: str):
        '''
        Generates a scatter plot and saves it as an image file.

        Parameters:
        - data_frame: DataFrame with data for this plot
        - plot_filename: The name of the plot image file to be saved

        Returns:
        - The path to the saved plot image
        '''
        fig, ax = plt.subplots(figsize=(12, 4), dpi=100)
        ax.set_title(self.title, fontsize=14)
        ax.set_xlabel(self.x_axis, fontsize=12)
        ax.set_ylabel(self.y_axis, fontsize=12)

        # Sort by x_axis for consistency
        data_frame = data_frame.sort_values(by=self.x_axis)

        sample_map = {
            'Matched Normal': '#6495ED',  
            'Tumor': '#C70039',            
        }

        # Group and plot by Sample Type
        for sample_type, group in data_frame.groupby('Sample Type'):
            x_vals = group[self.x_axis]
            y_vals = group[self.y_axis]

            # Create the scatter plot
            ax.scatter(
                x_vals, y_vals, 
                color=sample_map[sample_type], 
                label=sample_type,
                edgecolors='w', 
                alpha=0.6
            )

            # Draw median line
            y_vals = y_vals.dropna().tolist()
            if y_vals:
                y_med = median(y_vals)
                ax.axhline(
                    y=y_med, 
                    color=sample_map[sample_type],
                    linestyle='--', 
                    linewidth=1, 
                    label=sample_type + f" Median: {y_med:.2f}"
                )

        # Set y-axis range
        if self.hi is not None and self.lo is not None:
            ax.set_ylim(self.lo, self.hi)
        elif self.lo is not None:
            ax.set_ylim(bottom=self.lo)
        elif self.hi is not None:
            ax.set_ylim(top=self.hi)

        # Remove a_axis ticks
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False, top=False)

        # Add legend
        ax.legend(fontsize=10, loc='best')

        # Save the plot as a PNG file in the temp directory
        path = os.path.join(TEMP_DIR, plot_filename)
        plot_path = os.path.abspath(path)
        plt.tight_layout()
        plt.savefig(plot_path, format='png')
        plt.close(fig)

        return plot_path
