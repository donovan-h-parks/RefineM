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

from biolib.plots.abstract_plot import AbstractPlot

from refinem.plots.scatter import Scatter
from refinem.plots.gc_plots import GcPlots
from refinem.plots.td_plots import TdPlots
from refinem.plots.cov_perc_plots import CovPercPlots
from refinem.plots.gc_cov_plot import GcCovPlot
from refinem.plots.tetra_pca_plot import TetraPcaPlot
from refinem.plots.mpld3_plugins import LinkedBrush, Tooltip

from numpy import (mean as np_mean, abs as np_abs)


class CombinedPlots(AbstractPlot):
    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, genome_scaffold_stats,
             highlight_scaffold_ids, link_scaffold_ids,
             genome_stats,
             gc_dist, td_dist,
             gc_perc, td_perc, cov_perc):
        """Create figure containing distribution plots, tetranucleotide PCA plots, and GC vs. coverage plot.

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

        self.fig.set_size_inches(self.options.width, self.options.height)
        
        # create subplots depending on availability of coverage information
        if len(genome_stats.mean_coverage) >= 1:
            # note: the ordering here is specific and ensures
            # proper linked brushing
            axes_gc_dist = self.fig.add_subplot(241)
            axes_tetra_dist = self.fig.add_subplot(242)
            axes_coverage_dist = self.fig.add_subplot(243)
            axes_td_cov = self.fig.add_subplot(246)
            axes_pc1_cov = self.fig.add_subplot(247)
            axes_tetra_pc1_pc3 = self.fig.add_subplot(248)  
            axes_gc_coverage = self.fig.add_subplot(245)
            axes_tetra_pc1_pc2 = self.fig.add_subplot(244)            
        else:
            # note: the ordering here is specific and ensures
            # proper linked brushing
            axes_gc_dist = self.fig.add_subplot(231)
            axes_tetra_dist = self.fig.add_subplot(234)
            axes_deltaGC_td = self.fig.add_subplot(232)
            axes_pc1_td = self.fig.add_subplot(235)
            axes_tetra_pc1_pc2 = self.fig.add_subplot(233)
            axes_tetra_pc1_pc3 = self.fig.add_subplot(236)

        # create plots
        gc_plots = GcPlots(self.options)
        scatter, delta_gc, seq_len, label_plot_order = gc_plots.plot_on_axes(self.fig,
                                                                    genome_scaffold_stats,
                                                                    highlight_scaffold_ids,
                                                                    link_scaffold_ids,
                                                                    genome_stats.mean_gc,
                                                                    gc_dist,
                                                                    [gc_perc],
                                                                    None,
                                                                    axes_gc_dist,
                                                                    True)
                                                                                                     
        td_plots = TdPlots(self.options)
        _scatter, td, _, _ = td_plots.plot_on_axes(self.fig,
                                                    genome_scaffold_stats,
                                                    highlight_scaffold_ids,
                                                    link_scaffold_ids,
                                                    genome_stats.mean_signature,
                                                    td_dist,
                                                    [td_perc],
                                                    None,
                                                    axes_tetra_dist,
                                                    True)
                                
        if len(genome_stats.mean_coverage) >= 1:
            cov_per_plots = CovPercPlots(self.options)
            cov_per_plots.plot_on_axes(self.fig,
                                        genome_scaffold_stats,
                                        highlight_scaffold_ids,
                                        link_scaffold_ids,
                                        genome_stats.mean_coverage,
                                        [cov_perc],
                                        None,
                                        axes_coverage_dist,
                                        True)
                             
        if len(genome_stats.mean_coverage) >= 1:
            gc_cov_plot = GcCovPlot(self.options)
            _, gc, cov, _ = gc_cov_plot.plot_on_axes(self.fig,
                                                                     genome_scaffold_stats,
                                                                     highlight_scaffold_ids,
                                                                     link_scaffold_ids,
                                                                     genome_stats.mean_gc,
                                                                     genome_stats.mean_coverage,
                                                                     axes_gc_coverage,
                                                                     True) 

        tetra = TetraPcaPlot(self.options)
        _, pc1, _, _ = tetra.plot_on_axes(self.fig, 0, 1,
                                                          genome_scaffold_stats,
                                                          highlight_scaffold_ids,
                                                          link_scaffold_ids,
                                                          axes_tetra_pc1_pc2, True)

        tetra.plot_on_axes(self.fig, 0, 2,
                              genome_scaffold_stats,
                              highlight_scaffold_ids,
                              link_scaffold_ids,
                              axes_tetra_pc1_pc3, True)
             
        if len(genome_stats.mean_coverage) >= 1:
            td_cov_plot = Scatter(self.options)  
            td_cov_plot.plot_on_axes(self.fig,
                                        td, cov, label_plot_order,
                                        'tetranucleotide distance',
                                        'Coverage (mean = %.1f)' % np_mean(genome_stats.mean_coverage),
                                        highlight_scaffold_ids,
                                        link_scaffold_ids,
                                        axes_td_cov,
                                        True)
                                             
            pc1_cov_plot = Scatter(self.options)                                     
            pc1_cov_plot.plot_on_axes(self.fig,
                                        pc1, cov, label_plot_order,
                                        'PC 1 (%.1f%%)' % (tetra.variance[0] * 100),
                                        'Coverage (mean = %.1f)' % np_mean(genome_stats.mean_coverage),
                                        highlight_scaffold_ids,
                                        link_scaffold_ids,
                                        axes_pc1_cov,
                                        True)
        else:
            deltaGC_td_plot = Scatter(self.options)  
            deltaGC_td_plot.plot_on_axes(self.fig,
                                        np_abs(delta_gc), td, label_plot_order,
                                        '|delta GC| (mean = %.1f)' % genome_stats.mean_gc,
                                        'tetranucleotide distance',
                                        highlight_scaffold_ids,
                                        link_scaffold_ids,
                                        axes_deltaGC_td,
                                        True)
                                        
            pc1_td_plot = Scatter(self.options)                                     
            pc1_td_plot.plot_on_axes(self.fig,
                                        pc1, 
                                        td, 
                                        label_plot_order,
                                        'PC 1 (%.1f%%)' % (tetra.variance[0] * 100),
                                        'tetranucleotide distance',
                                        highlight_scaffold_ids,
                                        link_scaffold_ids,
                                        axes_pc1_td,
                                        True)

        mpld3.plugins.connect(self.fig, mpld3.plugins.LinkedBrush(scatter))

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
