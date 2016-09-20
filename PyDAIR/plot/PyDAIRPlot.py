from __future__ import division
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.cm as cmx
from PyDAIR.stats.PyDAIRStats import *
#from ggplot import *



class PyDAIRPlot:
    def __init__(self, stats, figure_style = 'fivethirtyeight'):
        self.__stats = stats
        self.__style = figure_style
        plt.style.use(figure_style)
    
       
    def __set_fig_labels(self, ax, main = '', xlab = '', ylab = ''):
        ax.set_title(main)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        return ax
    
    def __get_dev_params(self, fig_name = None, fig_format = None, fig_width = None, fig_height = None, fig_dpi = None,
                         def_fig_width = None, def_fig_height = None):
        if fig_name is None:
            fig_format = None
            fig_dpi = 96
        else:
            if fig_format is None:
                fig_format = os.path.splitext(fig_name)[1].replace('.', '')
            if fig_format not in ['pdf', 'eps', 'png', 'ps', 'svg', 'tiff', 'emf']:
                fig_format = 'pdf'
            fig_dpi = 300
        if os.path.splitext(fig_name)[1].replace('.', '') != fig_format:
            fig_name = fig_name + '.' + fig_format
           
        if fig_width is None:
            if def_fig_width is None:
                fig_width = 10
            else:
                fig_width = def_fig_width
        if fig_height is None:
            if def_fig_height is None:
                fig_height = 8
            else:
                fig_height = def_fig_height
        return [fig_name, fig_format, fig_width, fig_height, fig_dpi]
    
    
    
    
    
    def barplot_freq(self, gene, gene_names = None, prob = False, sort = True,
                     main = '', xlab = None, ylab = None,
                     fig_name = None, fig_format = None,
                     fig_width = None, fig_height = None, fig_dpi = None):
        if gene.lower() not in ['v', 'd', 'j']:
            raise ValueError('The \'gene\' should be one of \'v\', \'d\', and \'j\'.')
        
        
        xlab = 'Gene' if xlab is None else xlab
        if ylab is None:
            ylab = 'Frequency' if prob is False else ylab
            ylab = 'Probability' if prob is True else ylab
        
        sample_freqs = []
        sample_names = []
        for bsample in self.__stats.samples:
            sample_freqs.append(bsample.get_freq(gene))
            sample_names.append(bsample.name)
        
        freq_dataframe = pd.concat(sample_freqs, axis = 1).fillna(0)
        freq_dataframe.columns = sample_names
        if prob:
            freq_dataframe = freq_dataframe / freq_dataframe.sum(axis = 0)
        
        freq_dataframe.columns = sample_names
        if sort:
            freq_dataframe = freq_dataframe.ix[freq_dataframe.mean(axis = 1).sort_values(ascending = False).index]
        if gene_names is not None:
            freq_dataframe = freq_dataframe.reindex(gene_names).fillna(0)
        
        if self.__style == 'ggplot_dev':
            fig_name, fig_format, fig_width, fig_height, fig_dpi = self.__get_dev_params(
                    fig_name, fig_format, fig_width, fig_height, fig_dpi, int(0.8 * freq_dataframe.shape[0]), 6)
            freq_dataframe = pd.concat([pd.Series(freq_dataframe.index, name = 'Gene',
                                                  index = freq_dataframe.index), freq_dataframe], axis = 1)
            datadf = pd.melt(freq_dataframe, id_vars = ['Gene'], var_name  = 'Sample')
            p = ggplot(aes(x = 'Gene', y = 'value', fill = 'Sample'), data = datadf)
            p = p + geom_bar(position = 'dodge', stat = 'identity')
            #p = p + theme_bw()
            #p = p + scale_fill_brewer(type = 'qual', palette = 'Set1')
            #p = p + labs(title = main, x = xlab, y = ylab)
            p.save(fig_name)
            #ggsave(plot = p, filename = fig_name)
        else:
            if fig_name is None:
                fig = plt.figure()
            else:
                fig_name, fig_format, fig_width, fig_height, fig_dpi = self.__get_dev_params(
                    fig_name, fig_format, fig_width, fig_height, fig_dpi, int(0.8 * freq_dataframe.shape[0]), 6)
                fig = plt.figure(figsize = (fig_width, fig_height), dpi = fig_dpi)
        
            ax = fig.add_subplot(111)
            ax = self.__set_fig_labels(ax, main, xlab, ylab)
            freq_dataframe.plot.bar(ax = ax)
            
            if fig_name is None:
                plt.show()
            else:
                if fig_format != 'tiff':
                    if matplotlib.get_backend().lower() in ['agg', 'macosx']:
                        fig.set_tight_layout(True)
                    else:
                        fig.tight_layout()
                plt.savefig(fig_name, format = fig_format)
                plt.close()
    
    
    
    def hist_cdr3_len(self, xlim = None, prob = False,
                      main = '', xlab = None, ylab = None,
                      fig_name = None, fig_format = None,
                      fig_width = None, fig_height = None, fig_dpi = None):
        xlab = 'Length' if xlab is None else xlab
        if ylab is None:
            ylab = 'Frequency' if prob is False else ylab
            ylab = 'Probability' if prob is True else ylab
        
        sample_dists = []
        sample_names = []
        for bsample in self.__stats.samples:
            sample_dists.append(bsample.get_freq('cdr3_prot_len'))
            sample_names.append(bsample.name)
        
        dist_dataframe = pd.concat(sample_dists, axis = 1).fillna(0)
        dist_dataframe.columns = sample_names
        if prob:
            dist_dataframe = dist_dataframe / dist_dataframe.sum(axis = 0)
        dist_dataframe = dist_dataframe.set_index([[int(dx) for dx in dist_dataframe.index.values]])
        
        if xlim is None:
            xlim = [0, max([int(dx) for dx in dist_dataframe.index.values])]
        
        len_xlim = []
        for i in range(xlim[0], xlim[1] + 1):
            len_xlim.append(i)
        dist_dataframe = dist_dataframe.reindex(len_xlim).fillna(0)
        
        if self.__style == 'ggplot_dev':
            fig_name, fig_format, fig_width, fig_height, fig_dpi = self.__get_dev_params(
                    fig_name, fig_format, fig_width, fig_height, fig_dpi, int(0.8 * dist_dataframe.shape[0]), 6)
            dist_dataframe = pd.concat([pd.Series(dist_dataframe.index, name = 'Length',
                                                  index = dist_dataframe.index), dist_dataframe], axis = 1)
            datadf = pd.melt(dist_dataframe, id_vars = ['Length'], var_name  = 'Sample')
            p = ggplot(aes(x = 'Length', y = 'value', fill = 'factor(Sample)'), data = datadf)
            p = p + geom_bar(position = 'dodge', stat = 'identity')
            p = p + theme_bw()
            p = p + scale_fill_brewer(type = 'qual', palette = 'Set1')
            p = p + labs(title = main, x = xlab, y = ylab)
            p.save(fig_name)
        else:
            if fig_name is None:
                fig = plt.figure()
            else:
                fig_name, fig_format, fig_width, fig_height, fig_dpi = self.__get_dev_params(
                    fig_name, fig_format, fig_width, fig_height, fig_dpi, int(0.8 * dist_dataframe.shape[0]), 6)
                fig = plt.figure(figsize = (fig_width, fig_height), dpi = fig_dpi)
            ax = fig.add_subplot(111)
            ax = self.__set_fig_labels(ax, main, xlab, ylab)
            dist_dataframe.plot.bar(ax = ax)

            if fig_name is None:
                plt.show()
            else:
                if fig_format != 'tiff':
                    if matplotlib.get_backend().lower() in ['agg', 'macosx']:
                        fig.set_tight_layout(True)
                    else:
                        fig.tight_layout()
                plt.savefig(fig_name, format = fig_format)
                plt.close()
        
 
        








    
    
    def __scatter3d_freq(self, gene, prob, main, xlab, ylab, zlab, grid,
                         fig_name, fig_format, fig_width, fig_height, fig_dpi, sort):
        if gene != 'vdj':
            raise ValueError('Only vdj can be accepted.')
        if fig_name is not None:
            if len(fig_name) != len(self.__stats.samples):
                raise ValueError('The length of file names should be equal to the number of samples.')
        
        f = 0
        for bsample in self.__stats.samples:
            sample_name = bsample.name
            bfreq       = bsample.get_freq(gene)
            v_names = bfreq[0]
            d_names = bfreq[1]
            j_names = bfreq[2]
            freq    = bfreq[3]
            
            v_axis = []
            d_axis = []
            j_axis = []
            v_name_index = self.__get_gene_index(v_names)
            for v_name in v_names:
                v_axis.append(v_name_index[v_name])
            d_name_index = self.__get_gene_index(d_names)
            for d_name in d_names:
                d_axis.append(d_name_index[d_name])
            j_name_index = self.__get_gene_index(j_names)
            for j_name in j_names:
                j_axis.append(j_name_index[j_name])
        
            cm_heat = plt.cm.get_cmap('Blues')
        
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')
            ax.set_axis_bgcolor('white')
            ax.patch.set_facecolor('white')
            ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
            ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
            ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        
            xaxis = sorted(set(v_name_index.values()))
            yaxis = sorted(set(d_name_index.values()))
            zaxis = sorted(set(j_name_index.values()))
        
            if grid:
                ax.grid(False)
                for x in xaxis:
                    for y in yaxis:
                        ax.plot([x, x], [y, y], [min(zaxis), max(zaxis)], ls = 'solid', c = self.__col_grid, lw = 0.1)
                for x in xaxis:
                    for z in zaxis:
                        ax.plot([x, x], [min(yaxis), max(yaxis)], [z, z], ls = 'solid', c = self.__col_grid, lw = 0.1)
                for y in yaxis:
                    for z in zaxis:
                        ax.plot([min(xaxis), max(xaxis)], [y, y], [z, z], ls = 'solid', c = self.__col_grid, lw = 0.1)
        
            ax_scatter_obj = ax.scatter(v_axis, d_axis, j_axis, s = freq, c = freq, cmap = cm_heat)
            fig.colorbar(ax_scatter_obj)
            
            ax.set_xticks(xaxis)
            ax.set_yticks(yaxis)
            ax.set_zticks(zaxis)
            ax_xaxis = ax.set_xticklabels(v_name_index.keys())
            ax_yaxis = ax.set_yticklabels(d_name_index.keys())
            ax_zaxis = ax.set_zticklabels(j_name_index.keys())
            
            plt.setp(ax_xaxis, rotation = 90)
            plt.setp(ax_yaxis, rotation = 90)
            plt.setp(ax_zaxis, rotation = 0)
        
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.set_zlabel(zlab)
        
            ax.set_xlim(min(xaxis), max(xaxis))
            ax.set_ylim(min(yaxis), max(yaxis))
            ax.set_zlim(min(zaxis), max(zaxis))
        
            #ax.set_aspect('equal')
            if fig_name is None:
                plt.show()
            else:
                plt.savefig(fig_name[f], format = fig_format)
                plt.close()
            f += 1
    
    
    
    def plot_samplingresampling(self, plot_target = None, main = '', xlab = '', ylab = '',
                                fig_name = None, fig_format = None,
                                fig_width = None, fig_height = None, fig_dpi = None):
        pass   
        
        
        
        
    
    
    def plot_rarefaction(self, plot_target = None, main = '', xlab = '', ylab = '',
                         fig_name = None, fig_format = None,
                         fig_width = None, fig_height = None, fig_dpi = None):
        sample_diver = []
        sample_names = []
        max_samplingsize = 0
        for bsample in self.__stats.samples:
            sample_diver.append(bsample.div.rarefaction[plot_target].mean(axis = 1))
            sample_names.append(bsample.name)
            if max_samplingsize < bsample.div.rarefaction[plot_target].mean(axis = 1).max():
                max_samplingsize = bsample.div.rarefaction[plot_target].mean(axis = 1).max()
        
        if fig_name is None:
            fig = plt.figure()
        else:
            fig_name, fig_format, fig_width, fig_height, fig_dpi = self.__get_dev_params(
                fig_name, fig_format, fig_width, fig_height, fig_dpi, 6, 4)
            fig = plt.figure(figsize = (fig_width, fig_height), dpi = fig_dpi)
        
        ax = fig.add_subplot(111)
        ax = self.__set_fig_labels(ax, main, xlab, ylab)
        xlim = [0, max_samplingsize + 1000]
        for sample_idx in range(len(sample_diver)):
            sample_df = pd.DataFrame({'x': sample_diver[sample_idx].index.values,
                                      'y': sample_diver[sample_idx].values})
            sample_df.plot.line(x = 'x', y = 'y', ax = ax, label = sample_names[sample_idx], xlim = xlim)
        
        if fig_name is None:
            plt.show()
        else:
            if fig_format != 'tiff':
                if matplotlib.get_backend().lower() in ['agg', 'macosx']:
                    fig.set_tight_layout(True)
                else:
                    fig.tight_layout()
            plt.savefig(fig_name, format = fig_format)
            plt.close()
        
