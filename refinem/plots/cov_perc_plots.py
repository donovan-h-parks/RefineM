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

from biolib.plots.abstract_plot import AbstractPlot

from refinem.plots.mpld3_plugins import Tooltip


class CovPercPlots(AbstractPlot):
    """Histogram and scatterplot showing mean percent deviant of coverage profiles of scaffolds."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, genome_scaffold_stats, highlight_scaffold_ids, mean_coverage, cov_percs):
        """Setup figure for plots.

        Parameters
        ----------
        genome_scaffold_stats: d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : iterable
            Scaffolds in genome to highlight.
        mean_coverage : float
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

        self.plot_on_axes(self.fig, genome_scaffold_stats, highlight_scaffold_ids,
                          mean_coverage, cov_percs,
                          axes_hist, axes_scatter, True)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     genome_scaffold_stats, highlight_scaffold_ids,
                     mean_coverage, cov_percs,
                     axes_hist, axes_scatter, tooltip_plugin):
        """Create histogram and scatterplot.

        Parameters
        ----------
        figure : matplotlib.figure
          Figure on which to render axes.
        genome_scaffold_stats: d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : iterable
            Scaffolds in genome to highlight.
        mean_coverage : float
          Mean coverage profile of genome.
        cov_percs : iterable
          Coverage percentile values to mark on plot.
        """

        # calculate percent deviant of coverage profiles for each scaffold
        mean_deviations = []
        for stats in genome_scaffold_stats.values():

            mean_deviation = []
            for cov_genome, cov_scaffold in itertools.izip(mean_coverage, stats.coverage):
                mean_deviation.append(abs(cov_genome - cov_scaffold) * 100 / cov_genome)

            mean_deviations.append(np.mean(mean_deviation))

        # histogram plot
        axes_hist.hist(mean_deviations, bins=20, color=(0.5, 0.5, 0.5))
        axes_hist.set_xlabel('mean % deviation of coverage')
        axes_hist.set_ylabel('# scaffolds of %d' % len(mean_deviations))

        # prettify histogram plot
        self.prettify(axes_hist)

        # delta-GC vs sequence length (in kbp) scatterplot
        # highlighted points are put last in the list so they are plotted on top
        sorted_mean_deviations = []
        scaffold_lens = []
        colours = []
        labels = []
        for i, (scaffold_id, stats) in enumerate(genome_scaffold_stats.iteritems()):
            if scaffold_id not in highlight_scaffold_ids:
                sorted_mean_deviations.append(mean_deviations[i])
                scaffold_lens.append(stats.length / 1000.0)
                colours.append((0.7, 0.7, 0.7))
                labels.append('<small>{title}</small>'.format(title=scaffold_id))

        for i, (scaffold_id, stats) in enumerate(genome_scaffold_stats.iteritems()):
            if scaffold_id in highlight_scaffold_ids:
                sorted_mean_deviations.append(mean_deviations[i])
                scaffold_lens.append(stats.length / 1000.0)
                colours.append((1.0, 0, 0))
                labels.append('<small>{title}</small>'.format(title=scaffold_id))

        scatter = axes_scatter.scatter(sorted_mean_deviations, scaffold_lens, c=colours, s=self.options.point_size, lw=0.5)
        axes_scatter.set_xlabel('mean %% deviation of coverage (mean coverage = %.1f)' % np.mean(mean_coverage))
        axes_scatter.set_ylabel('Scaffold length (kbp)')

        _, ymax = axes_scatter.get_ylim()
        xmin, xmax = axes_scatter.get_xlim()

        # ensure y-axis include zero and covers all sequences
        axes_scatter.set_ylim([0, ymax])

        # ensure x-axis is set appropriately for sequences
        axes_scatter.set_xlim([xmin, max(xmax, max(cov_percs))])

        # draw vertical line for identifying outliers
        for cov_perc in cov_percs:
            axes_scatter.plot([cov_perc, cov_perc], [0, ymax], 'r--', lw=1.0, zorder=0)

        # prettify scatterplot
        self.prettify(axes_scatter)

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter
