###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import numpy as np

import matplotlib
import mpld3

from biolib.common import find_nearest

from refinem.plots.base_plot import BasePlot
from refinem.plots.mpld3_plugins import Tooltip


class GcPlots(BasePlot):
    """Create histogram and scatterplot showing GC distribution of scaffolds."""

    def __init__(self, options):
        """Initialize."""
        BasePlot.__init__(self, options)

    def plot(self, genome_scaffold_stats, highlight_scaffold_ids, link_scaffold_ids, mean_gc, gc_dist, percentiles_to_plot):
        """Setup figure for plots.

        Parameters
        ----------
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        mean_gc : float
          Mean GC of genome.
        gc_dist : d[gc][length][percentile] -> critical value
          GC distribution.
        percentiles_to_plot : iterable
          Percentile values to mark on plot.
        """

        # Set size of figure
        self.fig.clear()

        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.Reset(), mpld3.plugins.BoxZoom(), mpld3.plugins.Zoom())
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12, fmt='.1f'))

        self.fig.set_size_inches(self.options.width, self.options.height)

        axes_hist = self.fig.add_subplot(121)
        axes_scatter = self.fig.add_subplot(122)

        self.plot_on_axes(self.fig, genome_scaffold_stats,
                          highlight_scaffold_ids,
                          link_scaffold_ids,
                          mean_gc, gc_dist, percentiles_to_plot,
                          axes_hist, axes_scatter, True)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     genome_scaffold_stats,
                     highlight_scaffold_ids, link_scaffold_ids,
                     mean_gc, gc_dist, percentiles_to_plot,
                     axes_hist, axes_scatter, tooltip_plugin):
        """Create histogram and scatterplot.

        Parameters
        ----------
        figure : matplotlib.figure
          Figure on which to render axes.
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        mean_gc : float
          Mean GC of genome.
        gc_dist : d[gc][length][percentile] -> critical value
          GC distribution.
        percentiles_to_plot : iterable
          Percentile values to mark on plot.
        """

        # histogram plot
        if axes_hist:
            scaffold_gc = [stats.gc for stats in genome_scaffold_stats.values()]
            ylabel = '# scaffolds (out of %d)' % len(scaffold_gc)
            self.histogram(axes_hist, scaffold_gc, 20, 80, 2, '% GC', ylabel)

        # scatterplot
        xlabel = 'deviation in GC (mean GC = %.1f%%)' % mean_gc
        ylabel = 'Scaffold length (kbp)'

        scaffold_stats = {}
        for scaffold_id, stats in genome_scaffold_stats.iteritems():
            scaffold_stats[scaffold_id] = (stats.gc - mean_gc, stats.length / 1000.0)

        scatter, labels = self.scatter(axes_scatter,
                                         scaffold_stats,
                                         highlight_scaffold_ids,
                                         link_scaffold_ids,
                                         xlabel, ylabel)

        _, ymax = axes_scatter.get_ylim()
        xmin, xmax = axes_scatter.get_xlim()

        # draw vertical line at x=0
        axes_scatter.plot([0, 0], [0, ymax], linestyle='dashed', color=self.axes_colour, lw=1.0, zorder=0)

        # plot reference distributions
        closest_gc = find_nearest(np.array(gc_dist.keys()), mean_gc / 100)
        for percentile in percentiles_to_plot:
            # find closest distribution values
            temp_scaffold_len = gc_dist[closest_gc].keys()[0]
            d = gc_dist[closest_gc][temp_scaffold_len]
            gc_lower_bound_key = find_nearest(d.keys(), (100 - percentile) / 2.0)
            gc_upper_bound_key = find_nearest(d.keys(), (100 + percentile) / 2.0)

            xL = []
            xU = []
            y = []
            for window_size in gc_dist[closest_gc]:
                xL.append(gc_dist[closest_gc][window_size][gc_lower_bound_key] * 100)
                xU.append(gc_dist[closest_gc][window_size][gc_upper_bound_key] * 100)
                y.append(window_size / 1000.0)

            # sort by y-values
            sort_indexY = np.argsort(y)
            xL = np.array(xL)[sort_indexY]
            xU = np.array(xU)[sort_indexY]
            y = np.array(y)[sort_indexY]
            axes_scatter.plot(xL, y, 'r--', lw=1.0, zorder=0)
            axes_scatter.plot(xU, y, 'r--', lw=1.0, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axes_scatter.set_ylim([0, ymax])

        # ensure x-axis is set appropriately for sequences
        axes_scatter.set_xlim([xmin, xmax])

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter
