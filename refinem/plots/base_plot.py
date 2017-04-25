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

from matplotlib.collections import LineCollection

from biolib.plots.abstract_plot import AbstractPlot

from refinem.plots.mpld3_plugins import Tooltip

from numpy import (mean as np_mean)


class BasePlot(AbstractPlot):
    """Base plotting function used by many plots."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def point_properties(self, scaffold_stats, highlight_scaffold_ids, link_scaffold_ids):
        """Get visual properties for each point to be plotted.

        This includes organizing points such that those to be
        highlighted or linked are rendered last. This helps ensure
        this points will be visible.

        Parameters
        ----------
        scaffold_stats : d[scaffold_id] = [x, y]
          Gives the x and y coordinate of each scaffold.
        highlight_scaffold_ids : d[scaffold_id] -> color
          Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
          Pairs of scaffolds to link together.
        """

        # create list of scaffolds to highlight or link
        scaffolds_special_case = set()
        for id1, _c1, id2, _c2 in link_scaffold_ids:
            scaffolds_special_case.add(id1)
            scaffolds_special_case.add(id2)

        if highlight_scaffold_ids:
            scaffolds_special_case.update(highlight_scaffold_ids.keys())

        # highlighted and linked points are put last in the list so they are plotted on top
        x = []
        y = []
        colours = []
        labels = []
        for scaffold_id, stats in scaffold_stats.iteritems():
            if scaffold_id not in scaffolds_special_case:
                x.append(stats[0])
                y.append(stats[1])
                colours.append((0.7, 0.7, 0.7))
                labels.append('<small>{title}</small>'.format(title=scaffold_id))

        for scaffold_id, stats in scaffold_stats.iteritems():
            if scaffold_id in highlight_scaffold_ids:
                x.append(stats[0])
                y.append(stats[1])
                colours.append(highlight_scaffold_ids[scaffold_id])
                labels.append('<small>{title}</small>'.format(title=scaffold_id))

        links = []
        link_colors = []
        for id1, c1, id2, c2 in link_scaffold_ids:
            stats1 = scaffold_stats.get(id1, None)
            stats2 = scaffold_stats.get(id2, None)

            if stats1 == None or stats2 == None:
                continue

            x.append(stats1[0])
            x.append(stats2[0])
            y.append(stats1[1])
            y.append(stats2[1])

            colours.append(c1)
            colours.append(c2)

            labels.append('<small>{title}</small>'.format(title=id1))
            labels.append('<small>{title}</small>'.format(title=id2))

            links.append((stats1, stats2))
            
            # set to average color of end points and add alpha channel
            c = np_mean([c1 + [0.5], c2 + [0.5]], axis=0) 
            link_colors.append(c)

        return x, y, colours, labels, links, link_colors

    def histogram(self, axis, values, xmin, xmax, step, xlabel, ylabel):
        """Create histogram.

        Parameters
        ----------
        axis : matplotlib.axis
          Axis on which to render histogram.
        values : iterable
          Values from which to create histogram.
        xmin : float
            Minimum bin value.
        xmax : float
            Maximum bin value.
        step : float
          Size of bin.
        xlabel : str
          Label for x-axis.
        ylabel : str
          Label for y-axis.
        """
        # histogram plot
        gc_min, gc_max = 20, 80
        gc_step = 2

        axis.hist(values, bins=xrange(gc_min, gc_max, gc_step), color=(0.5, 0.5, 0.5))
        axis.set_xlabel('% GC')
        axis.set_ylabel('# scaffolds (out of %d)' % len(values))
        axis.set_xlim([gc_min, gc_max])

        # prettify histogram plot
        self.prettify(axis)

    def scatter(self, 
                    axis,
                    scaffold_stats,
                    highlight_scaffold_ids,
                    link_scaffold_ids,
                    xlabel, ylabel):
        """Create scatterplot with points rearranged.
        
        Points are rearranged to ensure visibility of highlighted
        points and links.

        Parameters
        ----------
        axis : matplotlib.axis
          Axis on which to render histogram.
        scaffold_stats : d[scaffold_id] = [x, y]
          Gives the x and y coordinate of each scaffold.
        highlight_scaffold_ids : d[scaffold_id] -> color
          Scaffolds in genome to highlight.
        link_scaffold_ids : list of scaffold pairs
          Pairs of scaffolds to link together.
        xlabel : str
          Label for x-axis.
        ylabel : str
          Label for y-axis.
        """

        x, y, colours, labels, links, link_colors = self.point_properties(scaffold_stats,
                                                                               highlight_scaffold_ids,
                                                                               link_scaffold_ids)

        scatter = axis.scatter(x, y, c=colours, s=self.options.point_size, lw=0.5)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)

        if links:
            line_segments = LineCollection(links, colors=link_colors, alpha=0.5)
            axis.add_collection(line_segments)

        # prettify scatterplot
        self.prettify(axis)

        return scatter, x, y, labels
        
    def scatter_fixed_order(self, axis,
                                x, y, labels,
                                highlight_scaffold_ids,
                                link_scaffold_ids,
                                xlabel, ylabel):
        """Create scatterplot with points in a fixed order.
        
        The current method for doing plotting is a little hacky
        in order to ensure correct linked brushing. Linked brushing
        requires points along an axis to be in identical ordering. 
        The function scatter() is problematic for this since it 
        reorders the points so those being highlighted or
        linked are plotted last (i.e., on top).
        
        This function does no such reordering and essentially
        expects any required reordering to already have been done.
        """
        
        colours = []
        plot_labels = []
        link_start = {d[0]:d[1] for d in link_scaffold_ids}
        link_end = {d[2]:d[3] for d in link_scaffold_ids}
        for label in labels:
            plot_labels.append('<small>{title}</small>'.format(title=label))
            
            if label in highlight_scaffold_ids:
                colours.append(highlight_scaffold_ids[label])
            elif label in link_start:
                colours.append(link_start[label])
            elif label in link_end:
                colours.append(link_end[label])
            else:
                colours.append((0.7, 0.7, 0.7))

        scatter = axis.scatter(x, y, c=colours, s=self.options.point_size, lw=0.5)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)

        if link_scaffold_ids:
            pts = {label:(x[labels.index(label)],y[labels.index(label)]) for label in labels}
            links = []
            link_colors = []
            for id1, c1, id2, c2 in link_scaffold_ids:
                pts1 = pts.get(id1, None)
                pts2 = pts.get(id2, None)

                if pts1 == None or pts2 == None:
                    continue
                    
                # set to average color of end points and add alpha channel
                c = np_mean([c1 + [0.5], c2 + [0.5]], axis=0)
                links.append((pts1, pts2))
                link_colors.append(c)
            
            if links:
                line_segments = LineCollection(links, colors=link_colors, alpha=0.5)
                axis.add_collection(line_segments)

        # prettify scatterplot
        self.prettify(axis)

        return scatter, plot_labels
        
    def plot_order(self, labels):
        """Strip added label information and return raw labels in plot order.
        
        Parameters
        ----------
        labels : list of str
            Labels in the order they were plotted, likely decorated with HTML code.
            
        Returns
        -------
        list of str
            Raw labels in order they were plotted.
        """
        
        raw_labels = []
        for label in labels:
            label = label[label.find('>')+1:label.rfind('<')]
            raw_labels.append(label)
            
        return raw_labels

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
