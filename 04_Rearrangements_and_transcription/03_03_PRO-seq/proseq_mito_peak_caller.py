"""
Locate initiation and termination sites based on PRO-seq data, can run in one sample (single) or multi-sample (multi) modes.
Returns a csv file where each row is a position 
"""

#%%
import random
import pandas as pd
import os
import glob
import dna_features_translator_class as dftc
import numpy as np
from importlib import reload
import random
import datetime
from itertools import chain
import utils
PATH = os.getcwd()
reload(utils)

#%%
def pro_peaks(df, positions = [], up_window = 200, down_window = 500, window = 100, specificity = 0.05, logging = True):
    """
    Recieve a PRO-seq annotated_pileup, and locate TIS and TERM sites which are defined as sites with significant increase
    or decrease in the downstream direction (relative to strand) respectively.
    
    Parameters
    ----------
    df : pd.DataFrame
        The sample df in annotated_pileup format
    up_window : int, default = 500
        The amount of bases to average over on the left side.
    down_window : int, default = 500
        The amount of bases to average over on the right side.
    window : int, default = 100
        The window to calculate a rolling mean over
    specificity : float, default = 0.05
        Defines the percentage over or under the mean the values will have to be to define a TIS or TERM respectively
    logging : bool, default = True
        Whether to create a log file
    
    Returns
    -------
    heavy_peaks_tis : list
        A list of heavy strand TIS peak calls
    heavy_peaks_term : list
        A list of heavy strand TERM valley calls
    light_peaks_tis : list
        A list of light strand TIS peak calls
    light_peaks_term : list
        A list of light strand TERM valley calls
    """


    if positions == []:
        positions = df.loc[df['Position'] % 2 != 0, 'Position']
    local_func = locals()
    if logging:
        if not os.path.isdir(os.path.join(PATH, 'logs')):
            os.mkdir(os.path.join(PATH, 'logs'))
        logid = str(random.random())[-7:-1]
        log = open(os.path.join(PATH, 'logs', f'log.{logid}.csv'), 'w')

    # Use a rolling window to smoothen the expression values
    df.pos_RPM = df.pos_RPM.rolling(window, min_periods = 1).mean()
    pos_mean = df.pos_RPM.mean()
    df.neg_RPM = df.neg_RPM.rolling(window, min_periods = 1).mean()
    neg_mean = df.neg_RPM.mean()
    # For the filteration of coding regions within the relevant strand, if there is no gene in the region, define as np.nan to prevent things like D-loop from being included as a transcribed gene

    gene_feature_names = ['Gene', 'gene', 'tRNA', 'rRNA']
    df.loc[~df['Feature'].isin(gene_feature_names), 'Strand'] = np.nan
    df[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']] = df[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']].fillna(0)

    heavy_peaks_tis = []
    heavy_peaks_term = []
    
    light_peaks_tis = []
    light_peaks_term = []
    now = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    if logging:
        log.write(f'Run started on {now}\nCurrently running on: {df.dataset.iloc[0]}\n')
        log.write(f'Function parameters:\n')
        {log.write(f'{k} = {v}\n') for k,v in local_func.items() if type(v) != pd.DataFrame and type(v) != list}
        log.write(f'Positive mean coverage = {pos_mean}\nPositive {specificity*100}% of mean coverage = {pos_mean*specificity}\nNegative mean coverage = {neg_mean}\nNegative {specificity*100}% of mean coverage = {neg_mean*specificity}\n')
        log_dict = {k:[] for k in  ['i', 'up_left', 'up_right', 'up_left_circle','up_right_circle', 'down_left','down_right', 'down_left_circle', 'down_right_circle', 'up_heavy_mean', 'down_heavy_mean', 'down_light_mean', 'up_light_mean','pos_mean','neg_mean', 'accepted']}
    # Iterate over all the mtDNA positions
    for i in positions:
        accepted = []
        # Upstream and downstream window definitions
        # If a window sized position to the left of the current position is outside of the mtDNA, loop around the entirety of the mtDNA.
        if i - up_window < 0:
            up_left = 0
            up_right = i
            up_left_circle = len(df) - abs(up_window - i)
            up_right_circle = len(df)
        # If not, only define up_left and up_right as the end and start of the upstream window
        else:
            up_left = i - up_window
            up_right = i
            up_right_circle = 9999999
            up_left_circle = 9999999
        
        # If a window sized position to the right of the current position is outside of the mtDNA, loop around the entirity of the mtDNA
        if i + down_window > len(df):
            down_left = i
            down_right = len(df)
            down_left_circle = 0
            down_right_circle = (i + down_window) - len(df)
        # If not, only define the up_left and up_right as the end and start of the upstream window
        else:
            down_left = i
            down_right = i + down_window
            down_left_circle = 9999999
            down_right_circle = 9999999
        
        # HEAVY STRAND BLOCK
        # Calculate the mean values for the upstream and downstream windows with the borders defined above
        if df.loc[df.Position == i, 'Strand'].iloc[0] != True:
            up_heavy = df.loc[
                    ((df.Position >= up_left) & (df.Position < up_right)) | (
                    (df.Position >= up_left_circle) & (df.Position < up_right_circle)), 'pos_RPM']
            down_heavy = df.loc[
                    ((df.Position > down_left) & (df.Position <= down_right)) | (
                     (df.Position > down_left_circle) & (df.Position <= down_right_circle)), 'pos_RPM']
            up_heavy_mean = up_heavy.mean()
            down_heavy_mean = down_heavy.mean()

            pos_mean = np.mean(up_heavy.append(down_heavy))

            # If the upstream window mean to current position is smaller than *specificity* times overall mean and larger than down stream window mean of window size then define as peak (TIS)
            if (up_heavy_mean < pos_mean*specificity) and (down_heavy_mean > pos_mean):
                heavy_peaks_tis.append(i)
                accepted.append('heavy_tis')
            # If the upstream window mean to current position is larger than overall mean and the down stream window to current position is smaller than *specificity* times overall mean define as peak (TERMINATION) 
            elif (up_heavy_mean >  pos_mean) and (down_heavy_mean < pos_mean*specificity):
                heavy_peaks_term.append(i)
                accepted.append('heavy_term')
        else:
            up_heavy_mean = 0
            down_heavy_mean = 0
            pos_mean = 0
        
        # LIGHT STRAND BLOCK
        # Calculate the light strand windows mean values with borders as defined above
        if df.loc[df.Position == i, 'Strand'].iloc[0] != False:
            down_light = df.loc[
                    ((df.Position >= up_left) & (df.Position < up_right)) | (
                    (df.Position >= up_left_circle) & (df.Position < up_right_circle)), 'neg_RPM']
            up_light = df.loc[
                    ((df.Position > down_left) & (df.Position <= down_right) | (
                     (df.Position > down_left_circle) & (df.Position <= down_right_circle))), 'neg_RPM']
            down_light_mean = down_light.mean()
            up_light_mean = up_light.mean()

            neg_mean = np.mean(down_light.append(up_light))


            if (up_light_mean < neg_mean*specificity) and (down_light_mean > neg_mean):
                accepted.append('light_tis')
                light_peaks_tis.append(i)
            elif (up_light_mean >  neg_mean) and (down_light_mean < neg_mean*specificity):
                accepted.append('light_term')
                light_peaks_term.append(i)
        else:
            down_light_mean = 0
            up_light_mean = 0
            neg_mean = 0

        if logging:
            for j in log_dict.keys():    
                log_dict[j].append(eval(j,locals()))


    if logging:
        now = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        log_df = pd.DataFrame(log_dict)
        log_df.to_csv(log)
        log.write(f'\nRun finished on: {now}\n')
        log.close()
    return (heavy_peaks_tis, heavy_peaks_term, light_peaks_tis, light_peaks_term)

def pro_refine_temp(df, pos_list, tis, gap, strand, barren_range = 500, reduc = 0.7):
    """
    Refine the predictions made by pro_peaks function, this is done in two ways:
        For positions that are less than {gap} size apart, pick out the max or min for initiation and termination respectively.
        The last one is picked if light strand, otherwise first one.
        Also, if the maximum value in {barren_range} size area around the peak is lower than the overall mean, remove it.
    Parameters
    ----------
    df : pd.DataFrame
        The sample df
    pos_list : list
        The list of position outputted by pro_peaks
    tis : bool
        True if input pos_list is composed of TIS sites
    gap : int
        The range to aggregate peaks in
    strand : bool
        True if heavy strand, False otherwise
    barren_range : int default = 300
        The range to check for an empty area in
    Returns
    -------
    final : list
        A list of the TIS and TERM positions after refinement
    """
    temp = []
    final = []
    
    # Define strand names
    strand_name = "pos" if strand else "neg"
    # Iterate over positions identified as peaks by pro_peaks function
    for i in pos_list:

        # If first position
        if len(temp) == 0:
            temp.append(i)
        if (i - temp[-1] > gap) or (i == pos_list[-1]):
            final_temp = df[df.Position.isin(temp)]
            if tis:
                final_pos = final_temp[final_temp[f'{strand_name}_RPM'] == final_temp[f'{strand_name}_RPM'].min()].Position.iloc[-1 if strand else 0]
                #final_pos = final_temp.Position.iloc[-1 if strand else 0]
            else:
                final_pos = final_temp[final_temp[f'{strand_name}_RPM'] == final_temp[f'{strand_name}_RPM'].min()].Position.iloc[0 if strand else -1]
            
            if final_pos - barren_range < 0:
                barren_pos = list(range(0, final_pos))
                barren_pos += list(range(len(df) - (abs(barren_range - final_pos))))
                
                barren_pos+= list(range(final_pos, final_pos + barren_range))
                barren_pos = list(set(barren_pos))
                #print(f'{final_pos} left side is smaller than mtDNA start. Start point: {barren_pos[0]}, end point: {barren_pos[-1]}')
            elif final_pos + barren_range > len(df):
                barren_pos = list(range(0, (final_pos + barren_range) - len(df)))
                barren_pos += list(range(final_pos, len(df)))
                
                barren_pos += list(range(final_pos - barren_range, final_pos))
                #print(f'{final_pos} right side is bigger than mtDNA start. Start point: {barren_pos[0]}, end point: {barren_pos[-1]}')
            else:
                barren_pos = list(range(final_pos - barren_range, final_pos + barren_range))
            #print(f'Maximum coverage within range: {df.loc[(df.Position.isin(barren_pos)), f"{strand_name}_RPM"].max()}', f'Threshold to pass: {df[f"{strand_name}_RPM"].mean() * reduc}')
            if df.loc[(df.Position.isin(barren_pos)), f'{strand_name}_RPM'].max() >= df[f'{strand_name}_RPM'].mean() * reduc:
                #print(f'Added {final_pos}!')
                final.append(final_pos)
            else:
                print(f'Skipped {final_pos}!')
                #final.append(final_pos)
            temp = []
        else:
            temp.append(i)
    return final

def refine_all_peaks(df, heavy_tis, heavy_term, light_tis, light_term, gap = 100, barren_range = 500, reduc = .8):
    """
    """
    df = df.copy()
    final_heavy_tis = pro_refine_temp(df, heavy_tis, tis = True, strand = True, barren_range = barren_range, reduc = reduc, gap = gap)
    final_heavy_term = pro_refine_temp(df, heavy_term, tis = False, strand = True, barren_range = barren_range, reduc = reduc, gap = gap)
    final_light_tis = pro_refine_temp(df, light_tis, tis = True, strand = False, barren_range = barren_range, reduc = reduc, gap = gap)
    final_light_term = pro_refine_temp(df, light_term, tis = False, strand = False, barren_range = barren_range, reduc = reduc, gap = gap)
    return final_heavy_tis, final_heavy_term, final_light_tis, final_light_term
    

def pro_mean_peaks(agged_list):
    """
    Grab the mean position of all positions less than 300 bases apart.
    Parameters
    ----------
    agged_list : list
        A nested list of all the refined peak calls of all samples
    Returns 
    -------
    
    """
    unnested = sorted([i for i in list(chain(*agged_list)) if i != np.nan])
    print('Currently running pro_mean_peaks function!','#################')
    if unnested == [] or unnested == [[]]:
        print('Empty list!')
        return []
    final = []
    changeindex = 0
    for i,v in enumerate(unnested):
        if v - unnested[i - 1] > 300 or i == len(unnested) - 1:
            print(unnested[changeindex : i + 1])
            final.append(int(round(np.nanmean(unnested[changeindex : i+1]), 0)))
            changeindex = i
    return final

def pair_transcripts(heavy_tis, heavy_term, light_tis, light_term):
    """
    Recieve 4 lists of heavy and light strand tis and terms, assemble them into transcripts in the following format: [start_pos, end_pos, strand (True/False)]
    The transcripts contain ALL possible combinations of initations and terminations.
    """
    transcripts = []
    (i.sort() for i in (heavy_tis, heavy_term, light_tis, light_term))
    for i, _ in enumerate(heavy_tis):
        for j,_ in enumerate(heavy_term):
            transcripts.append([heavy_tis[i], heavy_term[j], True])
    for i, _ in enumerate(light_tis):
        for j, _ in enumerate(light_term):
            transcripts.append([light_tis[i], light_term[j], False])
    return transcripts

def check_peak_calls(sample, heavy_tis, heavy_term, light_tis, light_term, gap, reduc, barren_range, gene_feature_names, cover_prop = 1):
    """
    Recieve unrefined peak calls, return True if the peaks are coherent.
    Coherent peaks cover all the genes in the particular strand.
    """
    sample = sample.copy()
    heavy_tis, heavy_term, light_tis, light_term = refine_all_peaks(sample, heavy_tis, heavy_term, light_tis, light_term, gap = gap, reduc = reduc, barren_range = barren_range)
    # Generate lists of positions with genes in them in both the Heavy and the Light strand
    pos_genes_range = sample.loc[(sample.Feature.isin(gene_feature_names)) & (sample.Strand == True), 'Position'].to_list()
    neg_genes_range = sample.loc[(sample.Feature.isin(gene_feature_names)) & (sample.Strand == False), 'Position'].to_list()
    #print(f'{len(pos_genes_range)=}') #type:ignore
    #print(f'{len(neg_genes_range)=}') #type:ignore
    pos_genes_range = random.sample(pos_genes_range, k = int(round(len(pos_genes_range) * cover_prop, None)))
    neg_genes_range = random.sample(neg_genes_range, k = int(round(len(neg_genes_range) * cover_prop, None)))
    #print(f'{len(pos_genes_range)=}') #type:ignore
    #print(f'{len(neg_genes_range)=}') #type:ignore

    transcripts = pair_transcripts(heavy_tis, heavy_term, light_tis, light_term)
    pos_t_range = []
    neg_t_range = []
    for t in transcripts:
        if t[0] >= t[1] and t[2] == True:
            window = abs(t[0] - len(sample)) + t[1]
        elif t[0] <= t[1] and t[2] == False:
            window = abs(t[1] - len(sample)) + t[0]
        else:
            window = abs(t[1] - t[0])
        if t[2] == True:
            pos_t_range += utils.mtdna_region(t[0], window = window, left = False, total = len(sample))
        else:
            neg_t_range += utils.mtdna_region(t[0], window = window, left = True, total = len(sample))
    return all((set(pos_genes_range).issubset(set(pos_t_range)), set(neg_genes_range).issubset(set(neg_t_range))))

def main(sample, up_window, down_window, specificity, gap, logging = False, reduc = 0.8, barren_range = 500, ladder = 0.05):
    """
    Main function, optimize specificity value, defined to be the minimum specificity that still yields atleast one TIS and TERM in every strand, return the TIS and TERM locations on both strands with the optimal specificity value

    Parameters
    ----------
    sample : pd.DataFrame
        A dataframe that contains a single processed_pileup sample
    up_window : int
        The left side window size
    down_window : int
        The right side window size
    specificity : float
        The specificity value to send to propeaks
    gap : int
        The maximum gap allowed to aggregated samples (with pro_refine_temp)
    logging : bool, default = False
        Whether to log the work done by pro_peak
    reduc : float, default = 0.8
        The reduction value sent to pro_refine_temp
    barren_range : int, default = 500
        The barren_range value sent to pro_refine_temp
    
    Returns
    -------
    final_heavy_tis : list
        List of refined heavy strand TIS calls
    final_heavy_term : list
        List of refined heavy strand TERM calls
    final_light_tis : list
        List of refined light strand TIS calls
    final_light_tis : list
        List of refined light strand TERM calls
    """

    pro = pd.read_csv(sample, index_col = 0)
    # For the filteration of coding regions within the relevant strand, if there is no gene in the region, define as np.nan to prevent things like D-loop from being included as a transcribed gene
    gene_feature_names = ['Gene', 'gene', 'tRNA', 'rRNA']
    pro.loc[~pro['Feature'].isin(gene_feature_names), 'Strand'] = np.nan
    pro[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']] = pro[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']].fillna(0)

    spec_temp = specificity
    runs = 0
    positions = pro.Position.to_list()
    while True:
        heavy_tis, heavy_term, light_tis, light_term = pro_peaks(pro.copy(), positions = positions, up_window = up_window, down_window = down_window, specificity = spec_temp, logging = logging)

        if (any([i == [] for i in (heavy_tis, heavy_term, light_tis, light_term)])) or\
           (spec_temp < 0.049) or\
           (not (check_peak_calls(pro, heavy_tis, heavy_term, light_tis, light_term, gene_feature_names = gene_feature_names, gap = gap, reduc = reduc, barren_range = barren_range))):
           #NOTE - Three main options that cause the specification reduction iterations to halt:
           #    1. The specification value is so high that no TERM/TIS were found for one of the options.
           #    2. The specification value is less than 0.05.
           #    3. The specification value resulted in incoherent transcriptional units.
 
            if runs == 0:
                print('The initial specificity is being used!')
                final_heavy_tis_unrefined, final_heavy_term_unrefined, final_light_tis_unrefined, final_light_term_unrefined = heavy_tis, heavy_term, light_tis, light_term
            print(f'Done running after {runs} iterations!\nFinal specificity is {round(spec_temp + ladder ,1)}!\nHeavy tis : {heavy_tis[0:1]} to {heavy_tis[-2:-1]}\nHeavy term : {heavy_term[0:1]} to {heavy_term[-2:-1]}\nLight tis : {light_tis[0:1]} to {light_tis[-2:-1]}\nLight term : {light_term[0:1]} to {light_term[-2:-1]}')
            break
        positions = list(set(heavy_tis + heavy_term + light_tis + light_term))
        final_heavy_tis_unrefined, final_heavy_term_unrefined, final_light_tis_unrefined, final_light_term_unrefined = heavy_tis, heavy_term, light_tis, light_term
        spec_temp -= ladder
        runs += 1
    final_heavy_tis, final_heavy_term, final_light_tis, final_light_term = refine_all_peaks(pro.copy(),final_heavy_tis_unrefined, final_heavy_term_unrefined, final_light_tis_unrefined, final_light_term_unrefined, gap = gap, barren_range = barren_range, reduc = reduc)
    return final_heavy_tis, final_heavy_term, final_light_tis, final_light_term

def main_on_directory(folder, up_window, down_window, specificity, gap, logging = False, reduc = 0.8, barren_range = 500, save_on_folder = True):
    """
    Run main function on every .csv file in a directory specified, create a dataframe which contains the position, type (TIS/TERM), dataset of origin, sample of peaks, strand of peaks and confidence level on peaks that is calculted based on the ratio of downstream coverage to downstream + upstream coverage

    Parameters
    ----------
    folder : str
        The directory that contains processed_pileup sample .csv files
    up_window : int
        The left side window to send to pro_peaks
    down_window : int
        The right side window to send to pro_peaks
    specificity : float
        The specificity value to send to pro_peaks
    gap : int
        The minimal gap required to separate peak aggregation in pro_refine_temp function
    logging : bool, default = False
        Whether to log process done by proseq_peaks
    reduc : float, default = 0.8
        The reduction value sent to proseq_refine_temp
    barren_range : int, default = 500
        The range to check for low expression on in pro_refine_temp
    save_on_folder : bool, default = True
        Whether to save the df output into system
    
    Returns
    -------
    main_df : pd.DataFrame
        A dataframe that contains all peak calls for every single sample, each row is a peak call and the columns reveal the position, type (TIS/TERM), location of peak, dataset it originates from, sample it originates from and the confidence score in the call.
    """
    files = glob.glob(os.path.join(folder, '*.csv'))
    print(f'Currently running on {os.path.split(folder)[-1]} dataset!')
    tis_term = ['TIS','TERM','TIS','TERM']
    strand = ['Heavy', 'Heavy', 'Light', 'Light']
    test_dict = {'Position':[], 'Type':[], 'Strand':[], 'Dataset':[], 'Sample':[], 'Confidence_score' : [], 'Dist_from_source' : []}
    for pro_path in files:
        pro = pd.read_csv(pro_path, index_col = 0)
        pro[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']] = pro[['neg_coverage','pos_coverage','neg_RPM','pos_RPM']].fillna(0)
        sample = os.path.split(pro_path)[-1]
        print(f'Currently running on {sample.replace(".csv","")} file!')
        final_heavy_tis, final_heavy_term, final_light_tis, final_light_term = main(pro_path, up_window, down_window, specificity, gap, logging, reduc, barren_range)
        dataset = os.path.split(folder)[-1]
        for i,k in enumerate((final_heavy_tis, final_heavy_term, final_light_tis, final_light_term)):
            test_dict['Position'] += k
            test_dict['Type'] += [tis_term[i] for _ in k]
            test_dict['Dataset'] += [dataset for _ in k]
            test_dict['Sample'] += [sample.replace('.csv','') for _ in k]
            test_dict['Strand'] += [strand[i] for _ in k]
            test_dict['Confidence_score'] += [utils.confidence_score(pos, pro, tis_term[i], True if strand[i] == 'Heavy' else False) for pos in k]
            test_dict['Dist_from_source'] += [utils.dist_from_source(pos, len(pro)) for pos in k]
    main_df = pd.DataFrame(test_dict)
    if save_on_folder:
        main_df.to_csv(os.path.join(os.path.dirname(folder), f'{dataset}.csv'), index = False)
    return main_df

def refined_main_on_directory(folder, up_window, down_window, specificity, gap, logging = False, reduc = 0.8, barren_range = 500, save_on_folder = True, min_rep = 50, agg_gap = 500, replace = False):
    """
    Run the main_on_directory function and then refine the output across samples.
    """
    if os.path.exists(os.path.join(PATH,  'data', os.path.basename(folder) + '.csv')) and replace == False:
        combined_df = pd.read_csv(os.path.join(PATH,  'data', os.path.basename(folder) + '.csv'))
    else:
        combined_df = main_on_directory(folder, up_window, down_window, specificity, gap, logging, reduc, barren_range, save_on_folder = True)
    combined_df = combined_df.rename({'Distance from source' : 'Dist_from_source'}, axis = 1)
    combined_df = combined_df.sort_values(by = 'Dist_from_source')

    new_df = {k:[] for k in combined_df.columns if k != 'Dist_from_source'}
    rand_sample_path = os.path.join(folder, random.choice([i for i in os.listdir(folder) if i.endswith('.csv')]))
    sample = pd.read_csv(rand_sample_path, index_col = 0)
    mtdna_len = len(sample)
    del sample
    total_samples = len(combined_df.Sample.unique())
    dataset = combined_df.Dataset.iloc[0]
    
    for strand in ['Light', 'Heavy']:
        for termtis in ['TIS', 'TERM']:
            templist = []
            cur_df = combined_df.loc[(combined_df.Type == termtis) & (combined_df.Strand == strand), :].sort_values(by = 'Dist_from_source')
            cur_df_len = len(cur_df)
            cur_df_list = cur_df.Dist_from_source.to_list()
            print(f'Currently running: {strand} - {termtis}')
            for i, v in enumerate(cur_df_list):
                try:
                    gapped = abs(cur_df_list[i + 1] - v) > agg_gap
                except IndexError: # Last element
                    gapped = True
                templist.append(v)
                if gapped:
                    pos_df = cur_df.loc[cur_df['Dist_from_source'].isin(templist), :]
                    templist = []
                    sample_rep = (len(pos_df.Sample.unique())/total_samples) * 100
                    if  sample_rep > min_rep:
                        new_df['Strand'].append(strand)
                        new_df['Type'].append(termtis)
                        new_df['Confidence_score'].append(round(pos_df.Confidence_score.mean(), 2))
                        new_df['Position'].append(int(round(np.average(pos_df.Dist_from_source, weights = pos_df.Confidence_score), 0)))
                        new_df['Sample'].append(sample_rep)
                        new_df['Dataset'].append(dataset)
                        if new_df['Position'][-1] < 0: # If position is to the left of origin
                            new_df['Position'][-1] = mtdna_len + 1 - abs(new_df['Position'][-1])
                    else:
                        print(f'Not enough representation for position, removing it!\nPosition={round(np.mean(pos_df.Dist_from_source), 0)}\nType={termtis}\nStrand={strand}')                     
    new_df = pd.DataFrame(new_df)
    if save_on_folder:
        new_df.to_csv(os.path.join(os.path.dirname(folder), f'{dataset}_refined.csv'), index = False)
    return new_df
    
