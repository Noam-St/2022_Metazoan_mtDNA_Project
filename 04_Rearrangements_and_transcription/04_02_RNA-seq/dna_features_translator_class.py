from typing import Iterable
from dna_features_viewer import BiopythonTranslator
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from annotation import record_check
from Bio import SeqIO
from os import getcwd
from sys import platform
import numpy as np 
from mycolorpy import colorlist as mcp

PATH = getcwd()

if platform != 'win32':
    slash = '/'
else:
    slash = '\\'

plt.style.use('ggplot')

class MitoTranslator(BiopythonTranslator):
  """
  A class of custom coloration and filtering for usage by dna_features_viewer
  """
  def compute_feature_color(self, feature):
    """
    Colors CDS blue, misc_feature grey, rRNA green and introns red
    otherwise gold
    """
    if feature.type == "CDS":
      return 'blue'
    elif feature.type == 'misc_feature':
      return 'grey'
    elif feature.type == 'rRNA':
      return 'green'
    elif feature.type == 'intron':
      return 'red'
    elif feature.type == 'tRNA':
      return 'purple'
    else:
      return 'gold'

  def compute_feature_label(self, feature):
    """
    If feature is not CDS misc rRNA or intron, dont label it
    """
    wanted = ['CDS','rRNA','intron']

    if feature.type not in wanted:
      return None
    else:
      return BiopythonTranslator.compute_feature_label(self, feature)

  def compute_filtered_features(self, features):
    """
    If feature not in wanted, dont graph it aswell"""
    wanted = ['CDS','misc_feature','rRNA','intron', 'tRNA']
    return [feature for feature in features if feature.type in wanted]

class MitoTranslatorT(BiopythonTranslator):
  """
  A class of custom coloration and filtering for usage by dna_features_viewer
  """
  def compute_feature_color(self, feature):
    """
    Colors CDS blue, misc_feature grey, rRNA green and introns red
    otherwise gold
    """
    if feature.type == "CDS":
      return 'blue'
    elif feature.type == 'misc_feature':
      return 'grey'
    elif feature.type == 'rRNA':
      return 'green'
    elif feature.type == 'intron':
      return 'red'
    elif feature.type == 'tRNA':
      return 'purple'
    else:
      return 'gold'

  def compute_feature_label(self, feature):
    """
    If feature is not CDS misc rRNA or intron, dont label it
    """
    wanted = ['CDS','rRNA','intron','tRNA']

    if feature.type not in wanted:
      return None
    else:
      return BiopythonTranslator.compute_feature_label(self, feature)

  def compute_filtered_features(self, features):
    """
    If feature not in wanted, dont graph it aswell"""
    wanted = ['CDS','misc_feature','rRNA','intron','tRNA']
    return [feature for feature in features if feature.type in wanted]

def modify_y(y, plot, log, to_window, window):
  """
  Helper function to coverage_graph. Takes the y-axis and modifies it according to the parameters sent to coverage_graph.

  Parameters
  ----------
  y : list
    List of y-axis values
  plot : str
    Name of the df column to be plotted
  log : bool
    If True, logarithmic y-axis
  
  Returns
  -------
  y : list
    List of y-axis values (modified)
  """
  if any(i in plot for i in to_window):
    y = y.rolling(window = window, min_periods = 1).mean()
  if log:
    y = np.log10(y)
  return y

def coverage_graph(
  org, sample, strand, addit_samples = None, labels = [''], plotlist = ['coverage'], prediction = None, transcripts = None, window = 100, scale = 'linear', log = False, peaks = None, peaks_alpha = 1, peaks_annotate = False, tis_term_df = None, pos_color = 'red', neg_color = 'blue', neut_color = 'purple', to_window = ['z', 'coverage', 'RPM'], hline = '', title_org = True, custom_title = '', savefig = '', plot_range = None, alpha = 1, upscaling = 1, fill = False, y_up_lim = False, style = 'default', despine = True):
  """
  Create two graphs aligned, annotations and coverage.

  Parameters
  ----------
  org : str
    Organism name
  sample : pd.DataFrame
    sample df created by annotation.py's main
  strand : bool 
    strand of coverage graph
  addit_samples : pd.DataFrame, None or list
    additional sample df created by annotation.py's main
  labels : list
    list of labels for the samples
  plotlist : list
    A list of column names to plot on the y axis against position.
  prediction : pd.DataFrame
    A df of the prediction dataframe created by annotation.py's main
  transcripts : list
    nested list in this format [initiation, termination, strand] for the graph annotations.
  window : int
    rolling function value
  scale : str
    scale of the y-axis
  log : bool
    take the log values of the plotlist
  peaks : list
    A list of x axis values (positions) to insert a vertical line in
  peaks_alpha : float
    alpha value for the peaks
  peaks_annotate : bool
    If True, annotate the peaks
  tis_term_df : pd.DataFrame
    A df of the tis_term dataframe created by annotation.py's main
  pos_color : str
    Positive expression coloration
  neg_color : str
    Negative expression coloration
  neut_color : str
    Neut expression coloration
  to_window : list
    List of sample df columns to apply the window function to
  hline : str
    If "mean" - include a horizontal dashed line at the mean value acrosss the y axis
  title_org : bool
    If True - include the organism name in the title
  custom_title : str
    If not empty - include a custom title 
  savefig : str
    If not empty - save the figure to this path
  plot_range : list
    If not empty - plot a range of values from the list
  alpha : float
    Transparency of the plot
  upscaling : float
    Upscaling factor for the plot
  fill : bool
    If True - fill the area under the plot
  
  Notes
  -----
  The plotlist can be a list of strings or a list of lists. If a list of lists, the first element of the list is the name of the plot and the second element is the name of the column in the sample df.
  Does not return anything, create graph inline
  """
  RANDOM_COLORS = mcp.gen_color('hsv', n = 100)
  POS_COLORS = ['black', '#ffcc00', '#abff00', '#00ff00', '#00ffab', '#00ffff', '#00abff', '#0000ff', '#ab00ff', '#ff00ff', '#ff00ab'] #colors for the positive strand
  NEG_COLORS = ['green', '#00e7ff', '#00ffd4', '#00ff9c', '#00ff6c', '#00ff3c', '#00ff00', '#3cff00', '#6cff00', '#9cff00', '#d4ff00'] #colors for the negative strand
  NEUT_COLORS = ['#ff9c00', '#ffc600', '#ffd400', '#ffd400', '#ffd400', '#ffd400', '#ffd400', '#ffd400', '#ffd400', '#ffd400', '#ffd400'] #colors for the neutral strand
  #TODO(Noam) - This 'kind of' works but is not very good. Need to go over line by line and fix this for multiple samples.
  plt.style.use(style)
  if strand == 'both':
    coverage_graph(org = org, sample = sample,strand = True, addit_samples = addit_samples, labels = labels, plotlist = plotlist, prediction = prediction,transcripts = transcripts,window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peaks_annotate =  peaks_annotate, tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =  custom_title, savefig =  savefig[0: -4] + '_heavy' + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine)
    coverage_graph(org = org, sample = sample,strand = False,addit_samples =  addit_samples,labels =  labels,  plotlist = plotlist, prediction = prediction,transcripts = transcripts,window =  window,scale = scale,log = log, peaks = peaks,peaks_alpha =  peaks_alpha,peaks_annotate =  peaks_annotate, tis_term_df = tis_term_df, pos_color = pos_color,neg_color =  neg_color,neut_color = neut_color,to_window =  to_window,hline = hline, title_org = title_org,custom_title =   custom_title,savefig =  savefig[0: -4] + '_light' + savefig[-4:], plot_range = plot_range,alpha =  alpha, upscaling = upscaling, fill = fill, y_up_lim=y_up_lim, style = style, despine = despine)
    return
  try: prediction = list(prediction)
  except TypeError: pass
  if plot_range == None:
    plot_range = (1, len(sample))
  sample = sample[sample.Position.isin(range(plot_range[0], plot_range[1] + 1))] # Limit the sample df to the plot range
  record_check(org) # Check if the organism is in the saved and if not create it
  record = SeqIO.read(f'{PATH}{slash}genbank_DB{slash}{org}.gbk', format = 'genbank') # Read the genbank file
  if strand == False:
    plotlist = ['neg_' + i for i in plotlist]
    color = neg_color
    color_list = NEG_COLORS   
  elif strand == True:
    plotlist = ['pos_' + i for i in plotlist]
    color = pos_color
    color_list = POS_COLORS
  elif strand == 'both':
    plotlist_temp = ['pos_' + i for i in plotlist]
    plotlist = ['neg_' + i for i in plotlist]
    plotlist += plotlist_temp
  else:
    color = neut_color
    color_list = NEUT_COLORS

  fig, axes = plt.subplots(len(plotlist) + 1, 1, figsize = (15, upscaling * 1.5 * (len(plotlist) + 1)), sharex = True, gridspec_kw = {'height_ratios': [1/upscaling] + [1 for _ in range(len(plotlist))]})    
  x = sample.Position
  for i, plot in enumerate(plotlist):
    if despine:
      sns.despine(ax = axes[i + 1], top = True, right = True, trim = False)
    i+=1
    if strand == 'both':
      if 'pos_' in plot:
        color = pos_color
        y = sample[plot]
      elif 'neg_' in plot:
        color = neg_color
        y = sample[plot]
    else:
      try: y = sample[plot]
      except KeyError:
        y = sample[plot.replace('neg_','').replace('pos_','')]
        print(f'Parameter {y} does not have a strand-specific representation!\n')
    y = modify_y(y, plot, log, to_window, window)
    axes[i].plot(x, y, color = color, label = labels[0], alpha = alpha)
    if fill: axes[i].fill_between(x, y, color = color, alpha = alpha/2)
    try:
      if addit_samples:
        if not isinstance(addit_samples, list): raise ValueError('addit_samples must be a list')
    except ValueError:
      addit_samples = [addit_samples]
    if addit_samples:
      for j, addit_sample in enumerate(addit_samples):
        if len(addit_samples) + 1 > len(labels): raise ValueError('The number of labels must match the number of addit_samples')
        addit_sample = addit_sample[addit_sample.Position.isin(range(plot_range[0], plot_range[1] + 1))]
        try: y_cur = addit_sample[plot]
        except KeyError:
          y_cur = addit_sample[plot.replace('neg_','').replace('pos_','')]
        y_cur = modify_y(y_cur, plot, log, to_window, window)
        axes[i].plot(x, y_cur, color = color_list[j], label = labels[j+1], alpha = alpha)
        if fill: axes[i].fill_between(x, y_cur, color = color_list[j], alpha = alpha/2)

    if i == len(plotlist):
      axes[i].set_xlabel('Position')
    axes[i].set_ylabel(plot.replace('_', ' ').replace('neg ','').replace('pos ','').capitalize())

    #If hline is requested, plot a horizontal line of the mean value
    if hline != '':
      if hline == 'mean':
        if plot in to_window:
          y_pos = y[y > 1].mean()
        else:
          y_pos = y[y > 0].mean()
      elif hline == 'median':
        if plot in to_window:
          y_pos = y[y > 1].median()
        else:
          y_pos = y[y > 0].median()
      elif type(hline) == int:
        y_pos = hline
      else: y_pos = 1
      axes[i].axhline(y = y_pos, color = 'black', alpha = 0.7, linestyle = ':')
    
    #Plot terminations/initiations in the transcripts list which is in this format: [start,end,strand]
    if transcripts:
      c=0
      for t in transcripts:
        if t[2] == strand:
          if c == 0:
            axes[i].axvline(x = t[0], color = 'lawngreen', label = 'Initiation')
            axes[i].axvline(x = t[1] ,color = 'lightcoral', label = 'Termination')
            c+=1
          else:
            axes[i].axvline(x = t[0], color = 'lawngreen')
            axes[i].axvline(x = t[1] ,color = 'lightcoral')
        elif strand != True and strand != False:
          if c == 0:
            axes[i].axvline(x = t[0], color = 'lawngreen')
            axes[i].axvline(x = t[1] ,color = 'lightcoral')
            c+=1
          else:
            axes[i].axvline(x = t[0], color = 'lawngreen')
            axes[i].axvline(x = t[1] ,color = 'lightcoral')

    #Plot whatever predictions are in the prediction list (designed to represent ML model predictions)
    if prediction:
      for i,p in enumerate(prediction):
        if p == 1:
          axes[i].axvline(x = i+1, color = 'green', alpha = 0.9, ymax = 0.3)
        elif p == 2:
          axes[i].axvline(x = i+1, color = 'red', alpha = 0.9, ymax = 0.3)
    
    #Plot whatever positions are in the peaks list
    if peaks:
      if type(peaks) != list: raise TypeError('Peaks must be a list of indices')
      for j in peaks:
        if j < 20: j = 25
        axes[i].axvline(x = j, color = 'black', alpha = peaks_alpha, ymax = 1)
    
    #Plot tis_term from the tis_term dataframe created by the proseq peak caller program
    if type(tis_term_df) == pd.DataFrame:
      if strand in [True, False]:
        tis_term_df = tis_term_df.loc[tis_term_df.Strand == ('Heavy' if strand else 'Light'), :]
      for _, row in tis_term_df.iterrows():
        checker = False
        loc = row.Position
        if loc < 20:
          loc += 20
          checker = True
        colour = 'green' if row.Type == 'TIS' else 'red'
        alpha = row.Confidence_score
        axes[i].axvline(x = loc, color = colour, alpha = alpha, ymax = 1)
        if peaks_annotate: # If peak annotation option is turned on, adds position values (x-axis) to the TIS/TERM.
          if plot_range[0] <= loc <= plot_range[1]:
            axes[i].text(x = loc + ((abs(plot_range[1] - plot_range[0])) * .005), y = y.max() * .7, s = loc - 20 if checker else loc, fontsize = 10, color = 'black', alpha = .5, fontstyle = 'oblique')

  plotsize = plot_range[1] - plot_range[0]
  if plotsize >= 1000:
    graphic_record = MitoTranslator().translate_record(record)
  else:
    graphic_record = MitoTranslatorT().translate_record(record)
  graphic_record = graphic_record.crop(plot_range)
  graphic_record.plot(ax = axes[0], with_ruler = False, annotate_inline = True, strand_in_label_threshold = 7)
  if plotsize <= 120:
    graphic_record.plot_sequence(ax = axes[0], y_offset = 2)
  if log: axes[1].set_ylabel('log(Coverage)')
  plt.yscale(scale)
  if title_org:
    plt.suptitle(f'{org}{"_" if custom_title != "" else ""}{custom_title}', fontsize = 12)
  else:
    plt.suptitle(f'{custom_title}', fontsize = 12)
  axes[1].set_title(('Heavy' if strand else 'Light') + ' strand')
  if strand == None:
    axes[1].set_title('Both strands')
  #axes[1].set_ylim(bottom = 0) #TODO Make sure this is a good idea
  if y_up_lim:
    axes[1].set_ylim(top = y_up_lim)
  plt.tight_layout()
  axes[1].legend()
  
  if savefig != '':
    fig.savefig(fname = savefig, dpi = 300)

def plot_gorder(org):
  """
  Simple function to plot only the organism - without any data
  """
  _, ax = plt.subplots(figsize = (15, 2.5))
  record_check(org)
  record = SeqIO.read(f'{PATH}{slash}genbank_DB{slash}{org}.gbk', format = 'genbank')
  graphic_record = MitoTranslator().translate_record(record)
  graphic_record.plot(ax = ax, with_ruler = True, annotate_inline = True, strand_in_label_threshold = 10)
  plt.suptitle(org)
  plt.tight_layout()


      

