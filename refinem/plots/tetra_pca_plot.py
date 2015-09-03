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

import numpy as np

from refinem.plots.base_plot import BasePlot
from refinem.plots.mpld3_plugins import LinkedBrush, Tooltip

from biolib.pca import PCA


class TetraPcaPlot(BasePlot):
    """Create a scatterplot of the first 2 tetranucleotide principal components."""

    def __init__(self, options):
        """Initialize."""
        BasePlot.__init__(self, options)

        self.pca_computed = False
        self.pc = None
        self.variance = None

    def pca(self, genome_scaffold_stats):
        """Perform PCA.

        Principal components are given in self.pca,
        and the variance in self.variance.

        Parameters
        ----------
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        """

        data_matrix = []
        for stats in genome_scaffold_stats.values():
            cols = len(stats.signature)
            data_matrix.append(stats.signature)

        data_matrix = np.reshape(np.array(data_matrix), (len(data_matrix), cols))

        pca = PCA()
        self.pc, self.variance = pca.pca_matrix(data_matrix, 3, bCenter=True, bScale=False)

        # ensure pc matrix has at least 2 dimensions
        if self.pc.shape[1] == 1:
            self.pc = np.append(self.pc, np.zeros((self.pc.shape[0], 2)), 1)
            self.variance = np.append(self.variance[0], np.ones(2))
        elif self.pc.shape[1] == 2:
            self.pc = np.append(self.pc, np.zeros((self.pc.shape[0], 1)), 1)
            self.variance = np.append(self.variance[0:2], np.ones(1))

        self.pca_computed = True

    def plot(self, genome_scaffold_stats, highlight_scaffold_ids, link_scaffold_ids):
        """Setup figure for tetranucleotide PCA plots.

        Parameters
        ----------
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        """

        # Set size of figure
        self.fig.clear()

        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.Reset(), mpld3.plugins.BoxZoom(), mpld3.plugins.Zoom())
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12, fmt='.1f'))

        self.fig.set_size_inches(self.options.width, self.options.height)

        axis_pc1_pc2 = self.fig.add_subplot(221)
        axis_pc3_pc2 = self.fig.add_subplot(222)
        axis_pc1_pc3 = self.fig.add_subplot(223)
        axis_variance = self.fig.add_subplot(224)

        scatter = self.plot_on_axes(self.fig, 0, 1,
                          genome_scaffold_stats,
                          highlight_scaffold_ids,
                          link_scaffold_ids,
                          axis_pc1_pc2, True)

        self.plot_on_axes(self.fig, 2, 1,
                          genome_scaffold_stats,
                          highlight_scaffold_ids,
                          link_scaffold_ids,
                          axis_pc3_pc2, True)

        self.plot_on_axes(self.fig, 0, 2,
                          genome_scaffold_stats,
                          highlight_scaffold_ids,
                          link_scaffold_ids,
                          axis_pc1_pc3, True)

        self.plot_variance(axis_variance)

        mpld3.plugins.connect(self.fig, LinkedBrush(scatter))

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plot_on_axes(self, figure,
                     pc_xaxis, pc_yaxis,
                     genome_scaffold_stats,
                     highlight_scaffold_ids,
                     link_scaffold_ids,
                     axis, tooltip_plugin):
        """Create tetranucleotide PCA scatterplot.

        Parameters
        ----------
        figure : matplotlib.figure
          Figure on which to render axes.
        pc_xaxis : int
          Principal component to plot on x-axis (zero indexed).
        pc_yaxis : int
          Principal component to plot on y-axis (zero indexed).
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        axis : matplotlib.axis
          Axis on which to render scatterplot.
        """

        if not self.pca_computed:
            self.pca(genome_scaffold_stats)

        scaffold_stats = {}
        for i, scaffold_id in enumerate(genome_scaffold_stats):
            scaffold_stats[scaffold_id] = (self.pc[i][pc_xaxis], self.pc[i][pc_yaxis])

        # scatterplot
        xlabel = 'PC %d (%.1f%%)' % (pc_xaxis + 1, self.variance[pc_xaxis] * 100)
        ylabel = 'PC %d (%.1f%%)' % (pc_yaxis + 1, self.variance[pc_yaxis] * 100)

        scatter, labels = self.scatter(axis, scaffold_stats,
                                       highlight_scaffold_ids, link_scaffold_ids,
                                       xlabel, ylabel)

        # tooltips plugin
        if tooltip_plugin:
            tooltip = Tooltip(scatter, labels=labels, hoffset=5, voffset=-15)
            mpld3.plugins.connect(figure, tooltip)

        return scatter

    def plot_variance(self, axis):
        """Create plot of variance captured by each principal component.

        Parameters
        ----------
        axis : matplotlib.axis
          Axis on which to render scatterplot.
        """

        axis.plot(np.arange(len(self.variance), dtype=int) + 1, np.cumsum(self.variance))
        axis.set_xlabel('Principal component')
        axis.set_ylabel('Percentage of cumulative variance')
        # axis.vlines(3, 0, 1.0, linestyle='dashed', color=self.axesColour, zorder=0, lw=0.5)
        axis.set_ylim([0, 1.02])
        axis.set_xlim([0, len(self.variance)])

        axis.get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        xticks = axis.get_xticks()
        if 0 in xticks and 1 not in xticks:
            xticks = np.append(np.array([1]), xticks[1:])
        axis.set_xticks(xticks)
