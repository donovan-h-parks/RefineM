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

import mpld3

from biolib.common import find_nearest
from biolib.genomic_signature import GenomicSignature

from refinem.plots.base_plot import BasePlot
from refinem.plots.mpld3_plugins import Tooltip


class TdPlots(BasePlot):
    """Create histogram and scatterplot showing tetranucleotide distribution (TD) of scaffolds."""

    def __init__(self, options):
        """Initialize."""
        BasePlot.__init__(self, options)

    def plot(self, genome_scaffold_stats,
             highlight_scaffold_ids, link_scaffold_ids,
             mean_signature, td_dist, percentiles_to_plot):
        """Setup figure for plots.

        Parameters
        ----------
        genome_scaffold_stats: d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        mean_signature : float
          Mean tetranucleotide signature of genome.
        td_dist : d[length][percentile] -> critical value
          TD distribution.
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

        self.plot_on_axes(self.fig,
                          genome_scaffold_stats,
                          highlight_scaffold_ids,
                          link_scaffold_ids,
                          mean_signature, td_dist, percentiles_to_plot,
                          axes_hist, axes_scatter, True)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     genome_scaffold_stats,
                     highlight_scaffold_ids, link_scaffold_ids,
                     mean_signature, td_dist, percentiles_to_plot,
                     axes_hist, axes_scatter, tooltip_plugin):
        """Create histogram and scatterplot.

        Parameters
        ----------
        figure : matplotlib.figure
          Figure on which to render axes.
        genome_scaffold_stats: d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        mean_signature : float
          Mean tetranucleotide signature of genome.
        td_dist : d[length][percentile] -> critical value
          TD distribution.
        percentiles_to_plot : iterable
          Percentile values to mark on plot.
        """

        # histogram plot
        genomic_signature = GenomicSignature(0)

        delta_tds = []
        for stats in genome_scaffold_stats.values():
            delta_tds.append(genomic_signature.manhattan(stats.signature, mean_signature))

        if axes_hist:
            axes_hist.hist(delta_tds, bins=20, color=(0.5, 0.5, 0.5))
            axes_hist.set_xlabel('tetranucleotide distance')
            axes_hist.set_ylabel('# scaffolds (out of %d)' % len(delta_tds))
            self.prettify(axes_hist)

        # scatterplot
        xlabel = 'tetranucleotide distance'
        ylabel = 'Scaffold length (kbp)'

        scaffold_stats = {}
        for i, (scaffold_id, stats) in enumerate(genome_scaffold_stats.iteritems()):
            scaffold_stats[scaffold_id] = (delta_tds[i], stats.length / 1000.0)

        scatter, labels = self.scatter(axes_scatter,
                                         scaffold_stats,
                                         highlight_scaffold_ids,
                                         link_scaffold_ids,
                                         xlabel, ylabel)

        _, ymax = axes_scatter.get_ylim()
        xmin, xmax = axes_scatter.get_xlim()

        # plot reference distributions
        for percentile in percentiles_to_plot:
            # find closest distribution values
            td_bound_key = find_nearest(td_dist[td_dist.keys()[0]].keys(), percentile)

            x = []
            y = []
            for window_size in td_dist:
                x.append(td_dist[window_size][td_bound_key])
                y.append(window_size / 1000.0)

            # sort by y-values
            sort_indexY = np.argsort(y)
            x = np.array(x)[sort_indexY]
            y = np.array(y)[sort_indexY]

            # make sure x-values are strictly decreasing as y increases
            # as this is conservative and visually satisfying
            for i in xrange(0, len(x) - 1):
                for j in xrange(i + 1, len(x)):
                    if x[j] > x[i]:
                        if j == len(x) - 1:
                            x[j] = x[i]
                        else:
                            x[j] = (x[j - 1] + x[j + 1]) / 2  # interpolate values from neighbours

                        if x[j] > x[i]:
                            x[j] = x[i]

            axes_scatter.plot(x, y, 'r--', lw=1.0, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axes_scatter.set_ylim([0, ymax])

        # ensure x-axis is set appropriately for sequences
        axes_scatter.set_xlim([xmin, xmax])

        # prettify scatterplot
        self.prettify(axes_scatter)

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter
