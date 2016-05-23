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

from refinem.plots.gc_plots import GcPlots
from refinem.plots.td_plots import TdPlots
from refinem.plots.cov_perc_plots import CovPercPlots
from refinem.plots.mpld3_plugins import LinkedBrush, Tooltip

import matplotlib
import mpld3

from biolib.plots.abstract_plot import AbstractPlot


class DistributionPlots(AbstractPlot):
    def __init__(self, options):
        """Initialize plot."""
        AbstractPlot.__init__(self, options)

    def plot(self, genome_scaffold_stats,
             highlight_scaffold_ids, link_scaffold_ids,
             genome_stats,
             gc_dist, td_dist,
             gc_perc, td_perc, cov_perc):
        """Create figure with the GC, tetranucleotide signature, and coverage distributions.

        Parameters
        ----------
        genome_scaffold_stats : d[scaffold_id] -> namedtuple of scaffold stats
          Statistics for scaffolds in genome.
        highlight_scaffold_ids : d[scaffold_id] -> color
            Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
            Pairs of scaffolds to link together.
        genome_stats : float
          Mean statistics for genome.
        gc_dist : d[gc][length][percentile] -> critical value
          GC distribution.
        td_dist : d[length][percentile] -> critical value
          TD distribution.
        gc_perc : float
          GC percentile value to mark on plot.
        td_perc : float
          TD percentile value to mark on plot.
        cov_perc : float
            Mean percent deviation to mark on plot.
        """

        # Set size of figure
        self.fig.clear()

        mpld3.plugins.clear(self.fig)
        mpld3.plugins.connect(self.fig, mpld3.plugins.Reset(), mpld3.plugins.BoxZoom(), mpld3.plugins.Zoom())
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12, fmt='.1f'))

        self.fig.set_size_inches(self.options.width*(2.0/3), self.options.height)

        # create subplots depending on availability of coverage information
        if len(genome_stats.mean_coverage) >= 1:
            axes_hist_GC = self.fig.add_subplot(321)
            axes_scatter_GC = self.fig.add_subplot(322)
            axes_hist_TD = self.fig.add_subplot(323)
            axes_scatter_TD = self.fig.add_subplot(324)
            axes_hist_cov_perc = self.fig.add_subplot(325)
            axes_scatter_cov_perc = self.fig.add_subplot(326)
        else:
            axes_hist_GC = self.fig.add_subplot(221)
            axes_scatter_GC = self.fig.add_subplot(222)
            axes_hist_TD = self.fig.add_subplot(223)
            axes_scatter_TD = self.fig.add_subplot(224)

        gc_plots = GcPlots(self.options)
        scatter, _, _, _ = gc_plots.plot_on_axes(self.fig,
                                                    genome_scaffold_stats,
                                                    highlight_scaffold_ids,
                                                    link_scaffold_ids,
                                                    genome_stats.mean_gc,
                                                    gc_dist,
                                                    [gc_perc],
                                                    axes_hist_GC,
                                                    axes_scatter_GC,
                                                    True)

        td_plots = TdPlots(self.options)
        td_plots.plot_on_axes(self.fig,
                                genome_scaffold_stats,
                                highlight_scaffold_ids,
                                link_scaffold_ids,
                                genome_stats.mean_signature,
                                td_dist,
                                [td_perc],
                                axes_hist_TD,
                                axes_scatter_TD,
                                True)

        if len(genome_stats.mean_coverage) >= 1:
            cov_per_plots = CovPercPlots(self.options)
            cov_per_plots.plot_on_axes(self.fig,
                                    genome_scaffold_stats,
                                    highlight_scaffold_ids,
                                    link_scaffold_ids,
                                    genome_stats.mean_coverage,
                                    [cov_perc],
                                    axes_hist_cov_perc,
                                    axes_scatter_cov_perc,
                                    True)

        mpld3.plugins.connect(self.fig, LinkedBrush(scatter))

        self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()

    def save_html(self, output_html):
        """Save figure as HTML.

        Parameters
        ----------
        output_html : str
            Name of output file.
        """

        html_script = Tooltip.script_global
        html_body = Tooltip.html_body

        AbstractPlot.save_html(self, output_html, html_script, html_body)
