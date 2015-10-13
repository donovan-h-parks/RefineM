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

import itertools

import mpld3

import numpy as np

from refinem.plots.base_plot import BasePlot
from refinem.plots.mpld3_plugins import Tooltip


class CovPercPlots(BasePlot):
    """Histogram and scatterplot showing mean percent difference of coverage profiles of scaffolds."""

    def __init__(self, options):
        """Initialize."""
        BasePlot.__init__(self, options)

    def plot(self, genome_scaffold_stats,
             highlight_scaffold_ids, link_scaffold_ids,
             mean_coverage, cov_percs):
        """Setup figure for plots.

        Parameters
        ----------
        genome_scaffold_stats: d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        mean_coverage : list of float
          Mean coverage profile of genome.
        cov_percs : iterable
          Coverage percentile values to mark on plot.
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
                          mean_coverage, cov_percs,
                          axes_hist, axes_scatter, True)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     genome_scaffold_stats,
                     highlight_scaffold_ids,
                     link_scaffold_ids,
                     mean_coverage, cov_percs,
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
        mean_coverage : list of float
          Mean coverage profile of genome.
        cov_percs : iterable
          Coverage percentile values to mark on plot.
        """

        # calculate percent difference of coverage profiles for each scaffold
        mean_perc_diffs = []
        for stats in genome_scaffold_stats.values():

            mean_perc_diff = []
            for cov_genome, cov_scaffold in itertools.izip(mean_coverage, stats.coverage):
                if cov_genome == 0:
                    mean_perc_diff.append(0)
                elif len(mean_coverage) >= 2:
                    mean_perc_diff.append(abs(cov_scaffold - cov_genome) * 100 / cov_genome)
                else:
                    mean_perc_diff.append((cov_scaffold - cov_genome) * 100 / cov_genome)

            mean_perc_diffs.append(np.mean(mean_perc_diff))

        # histogram plot
        if axes_hist:
            axes_hist.hist(mean_perc_diffs, bins=20, color=(0.5, 0.5, 0.5))
            if len(mean_coverage) >= 2:
                axes_hist.set_xlabel('absolute mean percentage error of coverage')
            else:
                axes_hist.set_xlabel('mean percentage error of coverage')

            axes_hist.set_ylabel('# scaffolds (out of %d)' % len(mean_perc_diffs))
            self.prettify(axes_hist)

        # scatterplot
        if len(mean_coverage) >= 2:
            xlabel = 'absolute mean percentage error of coverage (mean coverage = %.1f)' % np.mean(mean_coverage)
        else:
            xlabel = 'mean percentage error of coverage (mean coverage = %.1f)' % np.mean(mean_coverage)
        ylabel = 'Scaffold length (kbp)'

        scaffold_stats = {}
        for i, (scaffold_id, stats) in enumerate(genome_scaffold_stats.iteritems()):
            scaffold_stats[scaffold_id] = (mean_perc_diffs[i], stats.length / 1000.0)

        scatter, labels = self.scatter(axes_scatter,
                                         scaffold_stats,
                                         highlight_scaffold_ids,
                                         link_scaffold_ids,
                                         xlabel, ylabel)

        _, ymax = axes_scatter.get_ylim()
        xmin, xmax = axes_scatter.get_xlim()

        # draw vertical line at x=0
        axes_scatter.plot([0, 0], [0, ymax], linestyle='dashed', color=self.axes_colour, lw=1.0, zorder=0)

        # draw vertical line for identifying outliers
        for cov_perc in cov_percs:
            axes_scatter.plot([cov_perc, cov_perc], [0, ymax], 'r--', lw=1.0, zorder=0)
            if len(mean_coverage) == 1:
                axes_scatter.plot([-cov_perc, -cov_perc], [0, ymax], 'r--', lw=1.0, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axes_scatter.set_ylim([0, ymax])

        # ensure x-axis is set appropriately for sequences
        axes_scatter.set_xlim([xmin, max(xmax, max(cov_percs))])
        if len(mean_coverage) == 1:
            axes_scatter.set_xlim([min(xmin, -max(cov_percs)), max(xmax, max(cov_percs))])

        # prettify scatterplot
        self.prettify(axes_scatter)

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter
