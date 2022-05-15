import pandas as pd
import numpy as np
import os, glob, random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
import flanking_regions as fr
import dna_features_translator_class as dftc
from importlib import reload
from statannot import add_stat_annotation
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from Bio import SeqIO
from scipy.stats import mannwhitneyu
import snoop
from Bio import Entrez
from ast import literal_eval
import plotly.express as px

PATH = os.getcwd()
reload(fr)
reload(dftc)

def retrieve_gb(ID) -> str:
    """
    Retrieve GenBank record from NCBI given an ID.

    Parameters
    ----------
    ID : str
        The ID of the record to retrieve
    
    Returns
    -------
    filename : str
        The path to the record retrieved from NCBI
    """
    if not os.path.isdir(os.path.join(PATH, 'genebank.DB')):
        os.mkdir(os.path.join(PATH, 'genebank.DB'))
    filename = os.path.join(PATH, 'genebank.DB', ID + '.gbk')
    if not os.path.isfile(filename):
        print(f'Downloading {ID}\n')
        net_handle = Entrez.efetch(db ='nucleotide', id = ID, rettype = 'gb', retmode = 'text')
        out_handle = open(filename, 'w')
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print(f'Saved {ID}\n')
    return filename

def deseq(sample) -> pd.DataFrame:
    """
    Perform DESeq normalization on a sample.

    Parameters
    ----------
    sample : pd.DataFrame   
        The sample to be normalized
    
    Returns
    -------
    sample : pd.DataFrame
        The normalized sample
    """
    deseq_log = np.log(sample)
    deseq_log = deseq_log.replace([np.inf, -np.inf], np.nan).dropna(how = 'all')
    deseq_log['row_mean'] = deseq_log.mean(axis = 1)
    deseq_log_ratio = deseq_log.loc[:, deseq_log.columns != 'row_mean'].sub(deseq_log.row_mean, axis = 'rows')
    deseq_log_ratio_median = np.exp(deseq_log_ratio.median(axis = 0))
    sample = sample.divide(deseq_log_ratio_median)
    return sample

def deseq_from_long(sample, gene_col, count_col, sample_col) -> pd.DataFrame:
    """
    Convert long format to wide format, perform DESeq normalization, and return to long format.

    Parameters
    ----------
    sample : pd.DataFrame
        The long-format dataframe
    gene_col : str
        The gene name column
    count_col : str
        The raw count column
    sample_col : str
        

    Returns
    -------
    sample : pd.DataFrame
        The DESeq-normalized dataframe

    """
    sample = sample.pivot_table(values = count_col, index = gene_col, columns = sample_col)
    sample = deseq(sample)
    sample = sample.reset_index().melt(id_vars = gene_col, value_name = 'junc_deseq')
    return sample


def mtdna_region(pos, window, total, left, inclusive = True) -> str:
    """
    Designed for circular DNA, return a list of positions (INDEX 1 BASED) to the left or right side of pos

    Parameters
    ----------
    pos : int
        The current position to return a window around
    window : int
        The size of the window
    total : int
        The total size of the DNA
    left : bool
        True if the window should be to the left of the
    inclusive :
        (Default value = True)

    Returns
    -------
    positions : list
        A list of positions (INDEX 1 BASED) to the left or right side of pos
    """

    if type(left) != bool: raise TypeError(f'Parameter left must be either True (left side window) or False (right side window! Given {left} instead')
    if total < window or total < pos:
        raise ValueError(f'The total size of the DNA must be smaller than both the window and the position!\nParameters given:\nTotal = {total}\nPosition = {pos}\nWindow = {window}')
    
    positions = []
    if left:
        if pos - window < 1:
            positions += list(range(total - abs(window - pos) + + (1 if inclusive else 0), total + 1))
            positions += list(range(1, pos + (1 if inclusive else 0)))
        else:
            positions += list(range(pos - window + (1 if inclusive else 0), pos + (1 if inclusive else 0)))
    else: #right
        if pos + window > total:
            positions += list(range(1, (pos + window) - total + (0 if inclusive else 1)))
            positions += list(range(pos + (0 if inclusive else 1), total + 1))
        else:
            positions += list(range(pos + (0 if inclusive else 1), pos + window + 1))
    return positions

def plotly_sample(sample_path, plot_col = 'coverage', window = 100, peaks = False):
    """
    Create a plotly graph of a single sample expression, graphs both strands with negative values for light strand and positive values for the heavy strand.

    Parameters
    ----------
    sample_path :
        
    plot_col :
        (Default value = 'coverage')

    window :
        (Default value = 100)
    peaks :
        (Default value = False)

    Returns
    -------

    """
    sample = pd.read_csv(sample_path, index_col = 0)
    sample_name = os.path.split(sample_path)[-1].replace('.csv','')
    pos_col = 'pos_' + plot_col
    neg_col = 'neg_' + plot_col
    sample[[pos_col, neg_col]] = sample[[pos_col, neg_col]].fillna(0).rolling(window = window, min_periods = 1).mean()

    sub_sample = sample[['Position', pos_col, neg_col]]
    sub_sample = sub_sample.rename({pos_col : 'Heavy strand', neg_col : 'Light strand'}, axis = 1)
    sub_sample['Light strand'] = sub_sample['Light strand'] * -1
    sub_sample = sub_sample.melt(id_vars = 'Position', value_vars = ['Heavy strand', 'Light strand'], var_name = 'Strand', value_name = 'Coverage')
    fig = px.line(sub_sample, x = 'Position', y = 'Coverage', color = 'Strand', title = f'Expression profile of {sample_name}')
    if type(peaks) == pd.DataFrame:
        for i, peak in peaks.iterrows():
            fig.add_vline(
                x = peak.Position, \
                line_color = 'green' if peak.Type == 'TIS' else 'red',
                line_dash = 'dash' if peak.Strand == 'Heavy' else 'dot',
                annotation_text = peak.Strand)
    
    return fig

def ratio_score(i, df, strand, up_window = 250, down_window = 50, normalized = False) -> float:
    """
    Calculate the ratio of downstream coverage to downstream + upstream coverage, with downstream and upstream defined as right and left or left and right for the heavy strand and light strand respectively.

    Parameters
    ----------
    i : int
        The current position
    df : pd.DataFrame
        The dataframe containing the coverage data
    strand : str
        The strand to calculate the ratio for
    up_window : int
        The upstream window size
    down_window : int
        The downstream window size
    normalized : bool
        (Default value = False) If True, the ratio will be normalized by the mean of the upstream and downstream windows
    
    Returns
    -------
    float
        The ratio of downstream coverage to downstream + upstream coverage
    

    """
    strand_name = 'pos' if strand else 'neg'
    cov_col = 'RPM' if normalized else 'coverage' 
    total = df.Position.max()
    if strand:
        down_window_pos = mtdna_region(i, down_window, total, False)
        down_mean_cov = df.loc[df.Position.isin(down_window_pos), f'{strand_name}_{cov_col}'].mean() + 1
        
        up_window_pos = mtdna_region(i, up_window, total, True)
        up_mean_cov = df.loc[df.Position.isin(up_window_pos), f'{strand_name}_{cov_col}'].mean() + 1
    else:
        down_window_pos = mtdna_region(i, down_window, total, True)
        down_mean_cov = df.loc[df.Position.isin(down_window_pos), f'{strand_name}_{cov_col}'].mean() + 1
        
        up_window_pos = mtdna_region(i, up_window, total, False)
        up_mean_cov = df.loc[df.Position.isin(up_window_pos), f'{strand_name}_{cov_col}'].mean() + 1
    
    # Return the ratio of downstream coverage to downstream + upstream coverages, this ratio is supposed to be low for TERM and high for TIS
    return down_mean_cov/(up_mean_cov + down_mean_cov)

def confidence_score(i, df, type, strand, normalized = True) -> float:
    """
    Defined as the inverse of downstream and upstream ratio for TERM sites and the downstream and upstream ratio for TIS sites. (TIS must have large ratio_scores and TERM must have low ratio scores.

    Parameters
    ----------
    i : int
        The current position
    df : pd.DataFrame
        The dataframe containing the coverage data
    type : str
        The type of site to calculate the confidence score for
    strand : str
        The strand to calculate the ratio for
    normalized : bool
        (Default value = False) If True, the ratio will be normalized by the mean of the upstream and downstream windows
    
    Returns
    -------
    float
        The confidence score of the site
    
    """
    ratio_sc = ratio_score(i, df, strand, normalized = normalized)
    if type == 'TERM':
        confidence_sc = 1 - ratio_sc
    elif type == 'TIS':
        confidence_sc = ratio_sc
    return confidence_sc

def dist_from_source(i, total_len) -> int:
    """
    Calculate the mtDNA positions in terms of distance relative to 0 position (negative for left side and positive for right side).

    Parameters
    ----------
    i : int
        The current position
    total_len : int
        The total length of the mtDNA region
    
    Returns
    -------
    int
        The distance from the source
    """
    if i >= total_len/2:
        return ((total_len - i) + 1) * -1
    elif i < total_len/2:
        return i

def z_score(n, mean, std) -> float:
    """
    Compute Z-score.

    Parameters
    ----------
    n : int 
        The number to compute the Z-score for
    mean : float
        The mean of the distribution
    std : float
        The standard deviation of the distribution
    
    Returns
    -------
    float
        The Z-score
    
    References
    ----------
    https://en.wikipedia.org/wiki/Standard_score
    """
    try: return (n - mean) / std
    except ZeroDivisionError: return 0

def scale_min_max(n, mini, maxi) -> float:
    """
    Compute scaled min_max value

    Parameters
    ----------
    n : int 
        The number to compute the scaled value for
    mini : int
        The minimum value of the scaled range
    maxi : int
        The maximum value of the scaled range
    
    Returns
    -------
    float
        The scaled value
    
    References
    ----------
    https://en.wikipedia.org/wiki/Normalization_(statistics)
    """
    try: return (n - mini) / (maxi - mini)
    except ZeroDivisionError: return 0

def get_medians(df, sep_col, value_col) -> pd.DataFrame:
    """
    Grab the median of value_col grouped by sep_col

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the coverage data
    sep_col : str
        The column to group by
    value_col : str
        The column to compute the median for
    
    Returns
    -------
    pd.Series
        The median of value_col grouped by sep_col
    """
    mean_list = []
    for i in df[sep_col].unique():
        mean_list.append(df.loc[df[sep_col] == i, value_col].median())
    return mean_list

def num_sim(n1 : int, n2 : int) -> int:
    """
    calculates a similarity score between 2 numbers

    Parameters
    ----------
    n1 : int
        The first number
    n2 : int
        The second number
    
    Returns
    -------
    int
        The similarity score
    """
    return 1 - abs(n1 - n2) / (n1 + n2)

def get_org_name_from_refseq(ID):
    """
    Get the organism name from a refseq ID

    Parameters
    ----------
    ID : str
        The refseq ID
    
    Returns
    -------
    str
        The organism name
    """
    return '_'.join(SeqIO.read(retrieve_gb(ID), 'gb').features[0].qualifiers.get('organism')[0].replace(' ', '_').split('_')[0:2])

def get_org_filtered_expr(org, plotlist = ['RPM'], hline = 'mean', filtered = True, window = 50, tis_term_df = '', dataset_title = False) -> pd.DataFrame:
    """
    Iterate over each file in annotated_pileup folder, create a gigantic dataframe of all files, and plot
    the means of values in plotlist into a coverage_graph plot.

    Parameters
    ----------
    org :
        The organism to plot     
    plotlist :
        (Default value = ['RPM']) The list of columns to plot
    hline :
        (Default value = 'mean') The value to plot as a horizontal line
    filtered :
        (Default value = True) If True, the dataframe will be filtered by the mean of the values in plotlist
    window :
        (Default value = 50) The window size to use for the rolling mean
    tis_term_df :
        (Default value = '') The dataframe containing the TIS and TERM sites
    dataset_title :
        (Default value = False) If True, the title of the plot will be the name of the dataset
    """
    title = ''
    org_initials = ''.join([i[0] for i in org.replace('_',' ').split(' ')]).upper()
    all_files = glob.glob(os.path.join('data', 'annotated_pileup', f"*{org_initials}.csv"))
    if all_files == []:
        return None
    li = []
    c=0
    real_id = ''
    for filename in all_files:
        if org_initials not in filename: continue
        #Count organisms
        try: df = pd.read_csv(filename, index_col=0, header=0)
        except pd.EmptyDataError: continue
        try: id = df.Chromosome.dropna().iloc[0]
        except IndexError:
            continue
        if real_id == '':
            cur_org = '_'.join(SeqIO.read(retrieve_gb(id), 'gb').features[0].qualifiers.get('organism')[0].replace(' ', '_').split('_')[0:2])
            if cur_org == org:
                real_id = id
            else:
                continue
        else:
            if real_id != id:
                #print(f'Wrong ID {id}')
                continue
        df['org'] = org
        li.append(df)
        c+=1
    if li == []:
        print('None')
        return None
    #Combine all sample junction dfs to a single one
    annotated_df = pd.concat(li, axis=0, ignore_index=True)
    agg_dict = {'coverage':'median', 'RPM':'median','RPM':'std', 'z':'median',
                'ends_ratio':'median', 'Gene':'first', 'Strand':'first'}
    for i in annotated_df.dataset.unique():
        cur_df = annotated_df.loc[(annotated_df.dataset == i) & (annotated_df.org == org), :]
        cur_df = cur_df.groupby('Position').agg(agg_dict).reset_index().sort_values(by = 'Position')
        if len(li) == 1:
            cur_df = annotated_df
        if filtered: cur_df = fr.main(sample = cur_df, filtered = False, org = org.replace('_', ' '), flank = window,
                                      direction = 'downstream')

        if dataset_title:
            title = i                              
        dftc.coverage_graph(org.replace('_', ' '), cur_df, None, plotlist = plotlist, hline = hline, title_org = True, tis_term_df=tis_term_df, custom_title=title)

def falloff_score(row) -> float:
    """
    Compute the falloff score for a row

    Parameters
    ----------
    row : pd.Series
        The row to compute the falloff score for
    
    Returns
    -------
    float
        The falloff score
    """
    st1 = row['left_strand']
    st2 = row['right_strand']
    left_exp = row['left_z']
    right_exp = row['right_z']
    junc_exp = row['junc_tpm_z']
    if (st1 == True and st2 == True) or (st1 == True and st2 == False):
        return (left_exp - junc_exp)
    if (st1 == False and st2 == False) or (st1 == False and st2 == True):
        return (right_exp - junc_exp)

def junc_type(row) -> str:
    """
    Compute the junction type for a row

    Parameters
    ----------
    row : pd.Series
        The row to compute the junction type for    
    
    Returns
    -------
    str
        The junction type
    """
    st1 = row['left_strand']
    st2 = row['right_strand']

    if st1 == True and st2 == True:
        return 'A' # Both sorrounding genes are in the Heavy strand
    elif st1 == False and st2 == False:
        return 'B' # Both sorrounding genes are in the Light strand
    elif st1 == True and st2 == False:
        return 'C' # The left gene is in the Heavy strand and the right gene is in the Light strand (head to head)
    elif st1 == False and st2 == True:
        return 'D' # The left gene is in the Light strand and the right gene is in the Heavy strand (tail to tail)
    else:
        print('Wrong mode! choose either type or falloff!\n')

def custom_boxplot(test_type = 'Mann-Whitney', loc = 'inside', box_pairs = [('DSJ','SSJ')],
                   t_format = 'star', xlabel = 'Strand switched', ylabel = 'log(expression)', savefig = False, style = 'default', despine = True, **kwargs):
    """
    Recieve a FacetPlot object, create a boxplot with **kwargs and use add_stat_annotation to perform a statistical test.

    Parameters
    ----------
    test_type :
        (Default value = 'Mann-Whitney') The statistical test to perform
    loc :
        (Default value = 'inside') The location of the annotation
    box_pairs :
        (Default value = [('DSJ') ('SSJ')]) The list of box pairs to perform the statistical test on
    'SSJ')] :
        
    t_format :
        (Default value = 'star') The format of the annotation
    xlabel :
        (Default value = 'Strand switched') The xlabel of the boxplot
    ylabel :
        (Default value = 'log(expression)') The ylabel of the boxplot
    **kwargs : 
        The kwargs to pass to FacetGrid.map_dataframe
    """
    if 'DSJ' not in kwargs['data']['strand_switch'].unique(): return 
    plt.style.use(style)
    _, ax = plt.subplots(figsize = (4,4))
    sns.boxplot(ax = ax, **kwargs)
    ax.set(xlabel = xlabel, ylabel = ylabel, title = kwargs['data']['org'].iloc[0].replace('_', ' '))
    add_stat_annotation(
        ax, data=kwargs['data'], x=kwargs['x'], y=kwargs['y'],
        box_pairs=box_pairs,
        test = test_type, loc = loc,
        text_format = t_format, verbose = 1
        )
    plt.tight_layout()
    if despine:
        sns.despine(ax = ax, offset = 10, trim = False)
    if savefig:
        plt.savefig(os.path.join(PATH, 'figures', f'{kwargs["data"]["org"].iloc[0].replace("_", " ")}_{kwargs["data"]["dataset"].iloc[0]}.svg'), dpi = 300)
    

def skree_plot(pca) -> None:
    """
    Recieve a PCA object, create a skree plot

    Parameters
    ----------
    pca : sklearn.decomposition.PCA
        The PCA object to plot
    
    """
    features = range(pca.n_components_)
    plt.plot(features, pca.explained_variance_ratio_ * 100, color = 'black', linestyle = '--', marker = 'o')
    plt.xlabel('PCA Features')
    plt.ylabel('Variance %')
    plt.xticks(features)
    plt.tight_layout()
    plt.show()

def elbow_plot(pca_reduced_df):
    """
    Recieve a dataframe of all PCA components - create elbow plot

    Parameters
    ----------
    pca_reduced_df : pd.DataFrame
        The dataframe of all PCA components
    """
    ks = range(1, 10)
    inertias = []
    for k in ks:
        model = KMeans(n_clusters = k)
        model.fit(pca_reduced_df)
        inertias.append(model.inertia_)
    plt.plot(ks, inertias, '-o', color = 'black')
    plt.xlabel('N clusters (k)')
    plt.ylabel('Inertia')
    plt.xticks(ks)
    plt.tight_layout()
    plt.show()

def sample_loader(folder):
    """
    Recieve a folder, load all samples in the folder

    Parameters
    ----------  
    folder : str
        The folder to load samples from
    
    Returns
    -------
    sample_path : str
        The path to the sample
    sample : pd.DataFrame
        The sample dataframe
    sample_name : str
        The name of the sample
    """
    sample_path = os.path.join(PATH, 'proseq', 'data', folder, random.choice([i for i in os.listdir(os.path.join(PATH, 'proseq', 'data', folder)) if os.path.split(i)[-1].replace(".csv","") if i.endswith('.csv')]))
    
    sample = pd.read_csv(sample_path, index_col = 0)
    sample_name = os.path.split(sample_path)[-1].replace(".csv","")
    print(f'The choice is {sample_name}')
    return sample_path, sample, sample_name

def peaks_loader(name):
    """
    Recieve a name, load all peaks in the folder

    Parameters
    ----------
    name : str
        The name of the peaks to load
    
    Returns
    -------
    peaks : pd.DataFrame
        The peaks dataframe
    
    Example
    -------
    >>> peaks_loader('DSJ')

    """
    peaks_path = os.path.join(PATH, 'proseq', 'data',name + '_refined.csv')
    try:
        peaks = pd.read_csv(peaks_path)
        return peaks
    except FileNotFoundError:
        raise FileNotFoundError(f'No peaks for this study name! {name}')


def add_label_band(ax, top, bottom, label, *, spine_pos=-0.05, tip_pos=-0.02, orientation = 'vertical'):
    """Helper function to add bracket around y-tick labels.

    Parameters
    ----------
    ax : matplotlib.Axes
        The axes to add the bracket to
    top, bottom :
        The positions in *data* space to bracket on the y-axis
    label : str
        The label to add to the bracket
    spine_pos, tip_pos :
        The position in *axes fraction* of the spine and tips of the bracket.
        These will typically be negative
    top :     
    bottom :        
    spine_pos :
        (Default value = -0.05)
    tip_pos :
        (Default value = -0.02)
    orientation :
        (Default value = 'vertical')

    Returns
    -------
        
    """
    # grab the yaxis blended transform
    transform = ax.get_yaxis_transform()

    # add the bracket
    bracket = mpatches.PathPatch(
        mpath.Path(
            [
                [tip_pos, top],
                [spine_pos, top],
                [spine_pos, bottom],
                [tip_pos, bottom],
            ]
        ),
        transform=transform,
        clip_on=False,
        facecolor="none",
        edgecolor="k",
        linewidth=2,
    )
    ax.add_artist(bracket)

    # add the label
    txt = ax.text(
        spine_pos,
        (top + bottom) / 2,
        label,
        ha="right",
        va="center",
        rotation= orientation,
        clip_on=False,
        transform=transform,
    )

    return bracket, txt

def alter_cluster_model(gorder):
    """
    Receive a gene order, return True if it behaves according to the alternating gene clusters model.
    Alternative clustering model - groups of more than 2+ protein coding genes that are alternating between the heavy and the light strand.

    Parameters
    ----------
    gorder : list
        List of genes in the exact order they appear on the mtDNA.
    
    Returns
    -------
    : bool
        True if the gorder is arranged in AGC and False if it isnt
    """
    # Remove tRNA
    gorder = [i for i in gorder if 'trn' not in i]
    gcluster = []
    cluster_count = 0
    for gene in gorder:
        # If the current cluster list is empty, add the gene to it
        if gcluster == []:
            gcluster.append(gene)
            continue 
        # Check the current gene's strand and the last gene's strand
        cur_strand = '-' not in gene
        last_strand = '-' not in gcluster[-1]
        # If the current gene is the same strand as the last gene, add to cluster and keep going, otherwise end cluster.
        if cur_strand == last_strand:
            gcluster.append(gene)
        else:
            if len(gcluster) > 1:
                cluster_count += 1
            else:
                return False
            gcluster = []
            gcluster.append(gene)
    if cluster_count >= 2:
        return True
    else:
        return False

def sort_x(value, df, cat_a, cat_b, hue, x, y):
    """
    Sort a list of categories (x) according to the pvalue generated by comparing y between hue.
    """
    switched = df.loc[(df[x] == value) & (df[hue] == cat_a), y]
    not_switched = df.loc[(df[x] == value) & (df[hue] == cat_b), y]
    return mannwhitneyu(switched, not_switched, alternative='less').pvalue        

def plot_y_per_x_by_hue(
    df, x, y, figpath = None, hue = 'strand_switch',
    cat_a = 'DSJ',
    cat_b = 'SSJ',
    xlab = 'Organism',
    ylab = 'log(expression)',
    legend_loc = 'lower right',
    add_stats = True,
    test = 'Mann-Whitney',
    multiple_test_correction = True,
    verbose = 1,
    stat_y = None,
    title = '',
    style = 'default',
    despine = True):
    """
    Create a nice looking box plot to compare the values of y between hue for each x.
    Default is comparisons of the junc_tpm expression between DSJ and SSJ for each organism.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to plot
    x : str
        The column to plot on the x-axis
    y : str
        The column to plot on the y-axis
    figpath : str
        The path to save the figure to
    hue : str
        The column to group by
    cat_a : str
        The first category to compare
    cat_b : str
        The second category to compare
    xlab : str
        The label for the x-axis
    ylab : str
        The label for the y-axis
    legend_loc : str
        The location of the legend
    add_stats : bool
        Whether to add the p-value and t-statistics to the plot
    test : str
        The statistical test to use
    multiple_test_correction : bool
        Whether to use the Bonferroni correction
    verbose : int
        The verbosity level
    stat_y : str
        The column to use for the statistics
    title : str
        The title of the plot
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    
    """
    plt.style.use(style)
    if not stat_y: stat_y = y
    to_plot = df.replace('_', ' ', regex = True)
    order = sorted([i for i in to_plot[x].unique()], key = lambda w : sort_x(value = w, df = to_plot, cat_a = cat_a, cat_b = cat_b, hue = hue, x = x, y = y), reverse = False)
    box_pairs = [((org, cat_a),(org, cat_b)) for org in to_plot[x].unique()]
    _, ax = plt.subplots(figsize = (8,6))
    sns.boxplot(ax = ax, data = to_plot, y = y, x = x, hue = hue, orient = 'v', order = order)
    if add_stats:
        ax, _ = add_stat_annotation(ax, data = to_plot, y = stat_y,
                        x = x, hue = hue, box_pairs = box_pairs, loc = 'inside', test = test, verbose = verbose,
                    text_format = 'star', comparisons_correction = 'bonferroni' if multiple_test_correction else None, order = order)

    plt.xticks(rotation = 90, fontsize = 11)
    plt.xlabel(xlab, fontsize = 16)
    plt.ylabel(ylab, fontsize = 16)
    plt.title(title, fontsize = 16)
    plt.legend(loc = legend_loc)
    plt.tight_layout()
    if despine:
        sns.despine(ax = ax, right = True, top = True, trim = True)
    if figpath:
        plt.savefig(figpath, dpi = 300)


def gene_junc_type(row, df, org_col = 'org', junctype_col = 'strand_switch', DSJ_name = 'DSJ'):
    """
    """
    org_df = df.loc[(df[org_col] == row[org_col]) & (df[junctype_col] == DSJ_name), :]
    left_genes = list(set(org_df['left_gene'].squeeze()))
    right_genes = list(set(org_df['right_gene'].squeeze()))
    all_genes = left_genes + right_genes
    return 'DSJ' if row.gene in all_genes else 'SSJ'

def generate_2_x_2_conting(
    df,
    row_cat = 'Alter_model',
    comparison_cat = 'strand_switch',
    col_cat_A = 'DSJ',
    col_cat_B = 'SSJ',
    group_factor = 'org',
    cont_factor = 'junc_tpm'):
    """
    Generate a 2x2 chi_squared table for the following factors: significant/non significant and the row_cat.
    The p-values are generated by comparing the cont_factor values between two categories of comparison_cat for each group_factor category.
    """
    AGC_list = []
    non_AGC_list = []
    alpha = 0.05/len(df[group_factor].unique())
    for i in df[group_factor].unique():
        grp_1 = df.loc[(df[comparison_cat] == col_cat_A) & (df[group_factor] == i), cont_factor]
        grp_2 = df.loc[(df[comparison_cat] == col_cat_B) & (df[group_factor] == i), cont_factor]
        AGC = df.loc[df[group_factor] == i, row_cat].iloc[0].squeeze()
        sign = mannwhitneyu(grp_1, grp_2, alternative = 'less').pvalue < alpha
        if AGC:
            AGC_list.append(sign)
        else:
            non_AGC_list.append(sign)
    return np.array([[sum(AGC_list), len(AGC_list) - sum(AGC_list)], [sum(non_AGC_list), len(non_AGC_list) - sum(non_AGC_list)]])

def log2fc(list1, list2):
    """
    Calculate the log2 fold change of the means of list1 to list2.
    """
    return np.log2(np.mean(list1)/np.mean(list2))

def per_pos_log2fc(df_1, df_2, pos_col = 'Position', expr_col = 'RPM', pos_expr_col = 'pos_RPM', neg_expr_col = 'neg_RPM'):
    """
    Calculate the per position log2fold change of df_1 to df_2.
    """
    merged = df_1[[pos_col, expr_col, pos_expr_col, neg_expr_col]].merge(df_2[[pos_col, expr_col, pos_expr_col, neg_expr_col]], on = [pos_col], how = 'inner', suffixes = ('_1', '_2'))
    merged['log2fc'] = np.log2((merged[expr_col + '_1'] + 1)/(merged[expr_col + '_2'] + 1))
    merged['pos_log2fc'] = np.log2((merged[pos_expr_col + '_1'] + 1)/(merged[pos_expr_col + '_2'] + 1))
    merged['neg_log2fc'] = np.log2((merged[neg_expr_col + '_1'] + 1)/(merged[neg_expr_col + '_2'] + 1))
    return merged

def log2fc_sd(list1, list2):
    """
    Calculate the fold change of the std of list1 to list2.
    """
    # Calculate the log2 fold change for all combinations of the two lists
    log2fc_list = []
    for i in list1:
        for j in list2:
            log2fc_list.append(np.log2(i/j))
    return np.std(log2fc_list, ddof = 1)

def get_pairs(x):
    """
    Returns a list of pairs of genes in a list
    """
    try:x = literal_eval(x)
    except ValueError:pass
    return [f'{i}_{j}' for i, j in zip(x, x[1:])]

def combine_annotated_pileup(org, saveloc = 'combined_annotated_pileups'):
    """
    """
    org_initials = ''.join([i[0] for i in org.replace('_',' ').split(' ')]).upper()
    saveloc = os.path.join(PATH, 'data', saveloc)
    if os.path.exists(os.path.join(saveloc, f'{org}_agg_counts.csv')):
        print(f'Combined file for {org} already exists!')
        annotated_df = pd.read_csv(os.path.join(saveloc, f'{org}_agg_counts.csv'), index_col = 0)
        return annotated_df
    
    if not os.path.exists(saveloc):
        os.mkdir(saveloc)
    all_files = glob.glob(os.path.join('data', 'annotated_pileup', f"*{org_initials}.csv"))
    if all_files == []:
        return None
    li = []
    c=0
    real_id = ''
    for filename in all_files:
        #Count organisms
        try: df = pd.read_csv(filename, index_col=0, header=0) # Read in the csv file
        except pd.EmptyDataError: continue # If the file is empty, skip it
        try: id = df.Chromosome.dropna().iloc[0] # Get the RefSeq ID
        except IndexError:
            continue
        if real_id == '':
            cur_org = '_'.join(SeqIO.read(retrieve_gb(id), 'gb').features[0].qualifiers.get('organism')[0].replace(' ', '_').split('_')[0:2]) # Get the organism name based on the RefSeq ID
            if cur_org == org: # If the organism is the same as the one we are looking for, save the RefSeq ID
                real_id = id
            else:
                continue
        else:
            if real_id != id: # If the RefSeq ID is not the same as the one we are looking for, skip the file
                #print(f'Wrong ID {id}')
                continue
        df['org'] = org 
        li.append(df)
        c+=1
    if li == []:
        print('None')
        return None
    #Combine all sample junction dfs to a single one
    annotated_df = pd.concat(li, axis=0)
    agg_dict = {
        'coverage' : 'sum',
        'pos_coverage' : 'sum',
        'neg_coverage' : 'sum',
        'end_start_counts' : 'sum',
        'pos_end_start_counts' : 'sum',
        'neg_end_start_counts' : 'sum',
        'ends_ratio' : 'mean',
        'Feature' : 'first',
        'Gene' : 'first',
        'Length' : 'first',
        'Strand' : 'first',
        'dataset' : 'first'}
    annotated_df = annotated_df.groupby('Position').agg(agg_dict).reset_index()
    annotated_df['end_start_log_ratio'] = np.log10((annotated_df['end_start_counts'] + 1)/(annotated_df['coverage'] + 1))
    annotated_df['pos_end_start_log_ratio'] = np.log10((annotated_df['pos_end_start_counts'] + 1)/(annotated_df['pos_coverage'] + 1))
    annotated_df['neg_end_start_log_ratio'] = np.log10((annotated_df['neg_end_start_counts'] + 1)/(annotated_df['neg_coverage'] + 1))
    annotated_df['organism'] = org
    # Save the annotated_df to a csv file
    annotated_df.to_csv(os.path.join(saveloc, f'{org}_agg_counts.csv'), index=True)
    return annotated_df

def convert_to_asterisk(pvalue):
    """
    Convert a given p-value float to asterisks based on the scientific convention
    """
    if pvalue <= 0.0001:
        return '***'
    elif pvalue <= 0.001:
        return '**'
    elif pvalue <= 0.05:
        return '*'
    else:
        return ' '

def unite_samples_from_folder(folder):
  annot = os.path.join(PATH, 'data', folder, '*.csv')
  annot_df_list = [pd.read_csv(f, index_col = 0) for f in glob.glob(annot)]
  annot_df = pd.concat(annot_df_list, axis = 0)
  test = annot_df.groupby('Position').agg({'pos_end_start_counts' : 'sum', 'neg_end_start_counts' : 'sum', 'end_start_counts' :'sum', 'coverage' : 'sum', 'neg_coverage' : 'sum', 'pos_coverage' :'sum', 'Feature' : 'first', 'Gene' : 'first','Length' : 'first', 'Strand' : 'first', 'dataset' : 'first', 'RPM' : 'mean', 'neg_RPM' : 'mean', 'pos_RPM' : 'mean'}).reset_index()
  test['end_to_coverage'] =np.log2( (test['end_start_counts'] + 1) / (test['coverage'] + 1))
  test['pos_end_to_coverage'] =np.log2( (test['pos_end_start_counts'] + 1) / (test['pos_coverage'] + 1))
  test['neg_end_to_coverage'] =np.log2( (test['neg_end_start_counts'] + 1) / (test['neg_coverage'] + 1))
  return test