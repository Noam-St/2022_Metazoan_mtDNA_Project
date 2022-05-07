
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import random, os
from matplotlib.colors import CSS4_COLORS as colorations, Normalize
from sklearn.manifold import TSNE
import umap
import importlib
import consts
importlib.reload(consts)
from consts import COLOR_DICT, REP_DICT

PATH = os.getcwd()

def convert_lower_tri(A): #This code turns an lower triangular matrix into a mirrored full distance matrix.
    out = A.T + A
    np.fill_diagonal(out,np.diag(A))
    return out


def generate_model(mode, X, savepath, metric = 'precomputed', umap_n_neighbors = 15, umap_min_dist = 0.1, random_state = 69, tsne_perp = 30, tsne_n_iter = 5000):
    mode = mode.upper()
    if os.path.exists(savepath):
        reduced = np.load(savepath)
    elif mode == 'TSNE':
        model = TSNE(perplexity=tsne_perp, n_iter = tsne_n_iter, metric = metric, random_state = random_state)
        reduced = model.fit_transform(X)
        np.save(file = savepath, arr = reduced)
    elif mode == 'UMAP':
        model = umap.UMAP(n_neighbors = umap_n_neighbors, min_dist = umap_min_dist, metric = metric)
        reduced = model.fit_transform(X)
        np.save(file = savepath, arr = reduced)
    return reduced

def mapper(ph, classes, a):
    if ph in classes.keys():
        return classes[ph]
    else:
        return a
    
def get_prop_from_embedding(a, x_lower, x_upper, y_lower, y_upper, df, tax_level, tax_name):
    """
    Get percentage of organisms from tax_name in df that exist in a rectangle (of any 2d reduction technique) 
    """
    selection = df[ (x_upper > a[: ,0]) & (a[: ,0] > x_lower) & (y_upper > a[: ,1]) & (a[: ,1] > y_lower)]
    total = len(df.loc[df[tax_level] == tax_name])
    prop = len(selection.loc[selection[tax_level] == tax_name])
    return round((prop/total) * 100, 1)


def generate_labels(df, level = 'phylum', specific_tax = '', above_level = '', n_org_thresh = 30):
    """
    Generate a label system to color any taxonomic level defined in specific_tax
    Paramters
    ---------
    df : pd.DataFrame
        The input df that contains all the organisms.
    level : str (default = 'phylum')
        The taxonomic level to color
    specific_tax : str (default = '')
        Unnecessary for phylum level, if below phylum, this defines the name of the taxonomic group to label one level below.
    above_level : str (default = '')
        The name of the taxonomic level that is one above the one provided in specific_tax
    n_org_thresh : int
        The minimum amount of organisms that must be in a taxonomic group for it to not be labled as others.
    """
    if level == 'phylum':
        class_count = df[level].value_counts(ascending = False)
    else: 
        class_count = df[level][df[above_level] == specific_tax].value_counts(ascending = False)
    #Isolate descending value counts
    tester = {k:v for k,v in zip(class_count.index, class_count.values) if v > n_org_thresh}
    classes = {k:v for v,k in enumerate(tester.keys())}
    colors = []
    a = len(tester.keys())

    labels = df[level].apply(mapper, args = [classes, a])
    labels = np.array(labels).reshape(-1,1)

    classes.update({'Other':len(classes)})
    for i in classes.keys():
        if i in COLOR_DICT.keys():
            colors.append(COLOR_DICT[i])
        else:
            colorchoice = random.choice(list([i for i in colorations.keys() if ('grey' not in i) and (i not in colors)]))
            colors.append(colorchoice)
            
    classes =  {v:k for k,v in classes.items()}

    return labels, classes, colors

def make_text_label(df, org):
    """
    """
    offset = 0.01*(df.x.max() - df.x.min())
    x = df.loc[df.organism == org, 'x'].iloc[0]
    y = df.loc[df.organism == org, 'y'].iloc[0]
    return x + offset, y, org

def plot_from_embedding(
    df,
    file,
    tax = 'Chordata',
    taxo_level = 'phylum',
    above_level = 'phylum',
    add_rect = False,
    rect_taxo_level = 'phylum',
    rect_taxo = 'Chordata',
    rect_x = [0, 0],
    rect_y = [0, 0],
    show_others = True,
    figsize = (6,6),
    alpha = .9,
    n_org_thresh = 70,
    style = 'default',
    despine = False):

    plt.style.use(style)
    reduced = df[['x', 'y']].to_numpy()
    labels, classes, colors = generate_labels(df,
                                          level = taxo_level,
                                          specific_tax = tax,
                                          above_level = above_level,
                                          n_org_thresh=n_org_thresh)
    try: colormap = [COLOR_DICT[classes[int(i)]] for i in labels]
    except KeyError: pass
    test = get_prop_from_embedding(a = reduced,
                               x_lower = rect_x[0],
                               x_upper = rect_x[1],
                               y_lower = rect_y[0],
                               y_upper = rect_y[1],
                               df = df,
                               tax_level = rect_taxo_level,
                               tax_name = rect_taxo)
    _, ax = plt.subplots(figsize = figsize)
    c=0
    total = []
    for g, c in enumerate(colors):
        if not show_others:
            if classes[g] == 'Other': continue
        i = np.where(labels == g)[0]
        total += list(i)
        ax.scatter(reduced[i,0], reduced[i,1], label=classes[g],c = c, alpha = alpha, s = 7, edgecolors='none')
    if add_rect:
        rect = patches.Rectangle((rect_x[0] - 0.5, rect_y[0] - 0.5), abs(rect_x[1] - rect_x[0]) + 0.5, abs(rect_y[1] - rect_y[0]), linewidth=1, edgecolor='black', facecolor='none', label = f'{test}% of {rect_taxo}', linestyle = ':')
        ax.add_patch(rect)
    #ax.text(x = txt_x, y = txt_y, s = txt_s, fontdict = {'fontsize' : 10})
    ax.legend(bbox_to_anchor = (1.02, 1), markerscale = 2, fontsize = 11)

    if despine:
        sns.despine(ax = ax, top = True, right = True)
    plt.savefig(file, bbox_inches = 'tight', dpi = 300)
    plt.show()

def rect_slice(df, rect_x, rect_y):
    """
    """
    return df.loc[(df.x < rect_x[1]) & (df.x > rect_x[0]) & (df.y < rect_y[1]) & (df.y > rect_y[0]), :]

def top_in_rect(df, rect_x, rect_y):
    """
    """
    sliced = rect_slice(df, rect_x, rect_y)
    return sliced.value_counts('Gene_order', normalize = True)[0:1].to_dict()

def gene_symbol_to_number(to_rep,number2symbol = False):
    """ recives a str and dict, replaces str with any matching values in dict and returns dict reversed to replace string annotations
    with numerical annotations
    
    Parameters
    ----------
    to_rep : str
        string to replace
    number2symbol : bool
        if true, returns dict reversed
    
    Returns
    -------
    str
        string with replaced values
    dict
        dict with replaced values
    """
    counter=0
    bools=False
    if number2symbol:
        temp_rep = {v:k for k,v in REP_DICT.items()}
    else:
        temp_rep = REP_DICT
    if '*' in to_rep: #checks if this gene is a duplicate marked by *
        counter=to_rep.count('*')
    to_rep=to_rep.replace('*','')
    if '-' in to_rep[0]: #checks if the gene is in the complementary strand
        bools=True
        to_rep=to_rep.replace('-','',1)
    for k,v in temp_rep.items(): #iterates over the replacement dictionary and looks for matching gene strings
        if to_rep==v:
            to_rep=k
            break
    to_rep=to_rep+'*'*counter
    if bools: to_rep='-'+to_rep #If gene is in complementry, if the answer is yes 
    return to_rep

def list_of_genes_to_number(genes, number2symbol = False):
    return [gene_symbol_to_number(i, number2symbol = number2symbol) for i in genes]
