from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any
import pandas as pd
import matplotlib.pyplot as plt
from numpy import median
import time
import os

class Plot:
    def __init__(self, title: str, x_axis: str, y_axis: str, hi: int=-1, lo: int=-1) -> None:
        self.data = {       # data to be plotted
            "x": [],        # x-axis data (usually Sample IDs)
            "y": [],        # y-axis data
        }
        self.title = title  # title of plot
        self.axis = {       # names of axis
            "x": x_axis,
            "y": y_axis,
        }
        self.hi= hi         # upper bound of y-axis
        self.lo= lo         # lower bound of y-axis

    def generate_plot(self, name: str) -> str:
        """
        (str) -> str

        Generates a plot for the name column.
        Saves the plot to "Analysis_Reports/plots" and returns the path to the plot.

        Parameters
        -----------
        - name (str): y-axis values to be plotted

        """

        plt.figure(figsize=(9,2.5), dpi=80)     # fits 4 plots per page
        plt.scatter(self.data["x"], self.data["y"])
        plt.ylabel(self.axis["y"])              

        #plot median line
        plt.axhline(y=median(self.data["y"]), color='r')

        # set y-axis range
        if self.hi >= 0 and self.lo >= 0:
            plt.ylim(self.lo, self.hi)
        elif self.lo >= 0:
            plt.ylim(ymin=self.lo)
        elif self.hi >= 0:
            plt.ylim(ymax=self.hi)

        ax = plt.gca()
        ax.get_xaxis().set_visible(False)  # don't show x-axis

        working_dir = os.path.join(os.getcwd(), "temp")
        os.makedirs(working_dir, exist_ok=True)

        current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
        outputfile = os.path.join(working_dir, '{0}.{1}._plot.png'.format(name, current_time))

        plt.savefig(
            outputfile,
            bbox_inches="tight"
        )
        
        #close plot to save memory since we don't need it anymore
        plt.close()
        return outputfile

    def load_context(self, process_col) -> Dict[str, str]:
        """
        (str) -> dict[str, str]

        Loads and returns the context for the plot process_col, used in jinja2 templating

        Parameters
        -----------
        - process_col (str): name of the process and column to be plotted

        """
        context = {
            "title": self.title,        # title of plot
            "fig_path": self.generate_plot(process_col),    # path to plot generated
        }
        return context
    
    def add_data(self, val:str, id:str) -> None:
        """
        (str, str) -> None

        Add val and id to the data of the plot.

        Parameters
        -----------
        - val (str): val to be added (y-value)
        - id (str): id that corresponds to val (x-value)

        """
        self.data["x"].append(id)
        self.data["y"].append(val)

#SeqPlot class defines the plot for Raw Sequence Data and Call Ready Alignments Data
class SeqPlot(Plot):
    def __init__(self, title: str, x_axis: str, y_axis: str, hi: int=-1, lo: int=-1) -> None:
        self.data = {       # data to be plotted
        }
        self.title = title  # title of plot
        self.axis = {       # names of axis
            "x": x_axis,
            "y": y_axis,
        }
        self.hi= hi         # upper bound of y-axis
        self.lo= lo         # lower bound of y-axis
    
    def add_data(self, sample_type: str,val:str, id:str) -> None:
        """
        (str, str) -> None

        Add val and id to the data of the plot.

        Parameters
        -----------
        - val (str): val to be added (y-value)
        - id (str): id that corresponds to val (x-value)

        """
        if sample_type not in self.data.keys():
            self.data[sample_type] = {"x": [], "y": []}
        self.data[sample_type]["x"].append(id)
        self.data[sample_type]["y"].append(val)
    
    def generate_plot(self, name: str) -> str:
        """
        (str) -> str

        Generates a plot for the name column.
        Saves the plot to "Analysis_Reports/plots" and returns the path to the plot.

        Parameters
        -----------
        - name (str): y-axis values to be plotted

        """
        plt.figure(figsize=(14,5), dpi=80)     # fits 4 plots per page
        #draw plots for different sample types
        for stype in self.data.keys():
            sc = plt.scatter(self.data[stype]["x"], self.data[stype]["y"], label=stype)
            plt.axhline(
                y=median(self.data[stype]["y"]),
                c=sc.get_facecolors()[0].tolist(), #get colour of scatter plot and set median line to be same colour
                label=f"{stype} median"
            )
            
        plt.ylabel(self.axis["y"])              
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #set y-axis range
        if self.hi >= 0 and self.lo >= 0:
            plt.ylim(self.lo, self.hi)
        elif self.lo >= 0:
            plt.ylim(ymin=self.lo)
        elif self.hi >= 0:
            plt.ylim(ymax=self.hi)

        ax = plt.gca()
        ax.get_xaxis().set_visible(False)  # don't show x-axis

        working_dir = os.path.join(os.getcwd(), "temp")
        os.makedirs(working_dir, exist_ok=True)
        current_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
        outputfile = os.path.join(working_dir, '{0}.{1}._plot.png'.format(name, current_time))

        plt.savefig(
            outputfile,
            bbox_inches="tight"
        )
        
        #close plot to save memory since we don't need it anymore
        plt.close()
        return outputfile
