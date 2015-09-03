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

import matplotlib
import mpld3

from numpy import mean

from refinem.plots.base_plot import BasePlot
from refinem.plots.mpld3_plugins import Tooltip


class GcCovPlot(BasePlot):
    """Create a GC vs total coverage scatterplot."""

    def __init__(self, options):
        """Initialize."""
        BasePlot.__init__(self, options)

    def plot(self, genome_scaffold_stats, highlight_scaffold_ids, link_scaffold_ids, mean_gc, mean_coverage):
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
        mean_coverage : list of float
          Mean coverage profile of genome.
        """

        # Set size of figure
        self.fig.clear()

        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.Reset(), mpld3.plugins.BoxZoom(), mpld3.plugins.Zoom())
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12, fmt='.1f'))

        self.fig.set_size_inches(self.options.width, self.options.height)

        axis = self.fig.add_subplot(111)

        self.plot_on_axes(self.fig, genome_scaffold_stats,
                          highlight_scaffold_ids, link_scaffold_ids,
                          mean_gc, mean_coverage,
                          axis, True)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     genome_scaffold_stats,
                     highlight_scaffold_ids, link_scaffold_ids,
                     mean_gc, mean_coverage,
                     axis, tooltip_plugin):
        """Create GC vs. coverage scatterplot.

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
        mean_coverage : list of float
          Mean coverage profile of genome.
        axis : matplotlib.axis
          Axis on which to render scatterplot.
        """

        scaffold_stats = {}
        for scaffold_id, stats in genome_scaffold_stats.iteritems():
            scaffold_stats[scaffold_id] = (stats.gc, mean(stats.coverage))

        # scatterplot
        xlabel = 'GC (mean = %.1f%%)' % mean_gc
        ylabel = 'Coverage (mean = %.1f)' % mean(mean_coverage)

        scatter, labels = self.scatter(axis, scaffold_stats,
                                       highlight_scaffold_ids, link_scaffold_ids,
                                       xlabel, ylabel)

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter
