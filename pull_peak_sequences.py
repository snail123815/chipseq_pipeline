from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import concurrent.futures 
import sys  # for error handling

from funcs import getSpanFetures 
import pandas as pd


def readPeak(file, thresh=150):
    ''' Three kind of dataframe generated
    name, start, end, abs_summit, fold_enrichment
    name, start, end, abs_summit, fold_enrichment_A, fold_enrichment_B
    name, start, end, log10_likely (of the peak presented in this file)
    thresh is the threshold of the -LOG10(pvalue)
    '''
    if file.endswith('.xls'):
        # direct peak calling result
        data = pd.read_csv(
            file, delimiter='\t',
            comment='#', index_col='name',
            usecols=[
                'name', 'start', 'end', 'abs_summit','-LOG10(pvalue)', 'fold_enrichment'
                ]
            )
        data = data[data['-LOG10(pvalue)']>=thresh]
    elif file.endswith('.bed'):  # different peaks called by Macs2
        data = pd.read_csv(file, delimiter='\t', skiprows=1,
                           header=None, usecols=[1, 2, 3, 4], index_col=2)
        if 'common' in file:
            data.columns = ['start', 'end', 'likely_difference']
        else:
            data.columns = ['start', 'end', 'log10_likely']
        data.index.name = 'name'
        # TODO Thresh
    elif file.endswith('.tsv'):
        # peaks by comparing peak calling result
        if 'common_peaks' in file:
            cols = ['name', 'start', 'end', 'abs_summit',
                    'fold_enrichment_A', 'fold_enrichment_B']
        else:
            cols = ['name', 'start', 'end', 'abs_summit','-LOG10(pvalue)', 'fold_enrichment']
        data = pd.read_csv(file, delimiter='\t',
                           usecols=cols,
                           index_col='name')
        # TODO Thresh
    else:
        raise NameError
    data = data[~data.index.duplicated(keep='first')]
    return data


def slice(sourceSeq, location, id=None):
    #start, end = peakDict[peak]
    start, end = location
    try:
        sliceFull = sourceSeq[start:end]
    except:
        print(sys.exc_info()[0])
        print(start, end)
        exit()
    # Expand features if the cut location is inside features
    sliceFull.features.extend(getSpanFetures(sourceSeq, start, end))
    descrip = []
    if len(sliceFull.features) > 0:
        for feat in sliceFull.features:
            if feat.type == 'gene':
                try:
                    descrip.append(feat.qualifiers['locus_tag'][0])
                except:
                    pass
    sliceSeq = SeqRecord(sourceSeq.seq[start:end])
    if isinstance(id, type(None)):
        sliceSeq.id = f'{sourceSeq.id}_{start}-{end}'
    else:
        sliceSeq.id = id
    sliceSeq.description = '-'.join(descrip).replace(' ', '_')
    return sliceSeq
# slice



def getSeq(sourceSeq, peakDict):
    sourceSeq = sourceSeq[:] # make a copy
    resultSeqs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for peak in peakDict:
            future.append(executor.submit(slice, sourceSeq, peakDict, peak))
        for future in concurrent.futures.as_completed(futures):
            sliceSeq = future.result()
            resultSeqs.append(sliceSeq)
    return resultSeqs
# getSeq


def singleFilter(filter, method):
    if filter not in ['fold', 'length', 'summit', 'likely', None]:
        raise Exception(
            f"Filter error, {filter} ['fold', 'length', 'summit', 'likely', None]")
    print(f'Filtering: {filter} | {method}')
    if filter == 'fold':
        # file ext in ['.tsv', '.xls']
        if type(method) == list:
            # make sure the method returned by singleFilter() still a list
            subMethod = method[0]
        else:
            subMethod = method
        if 'common' in file:
            if subMethod not in ['max', 'min']:
                raise Exception(
                    f"Method error. For {filter} in {file} you should use ['max', 'min'], \n while {subMethod} has been passed.")
        else:
            if subMethod != 'single':
                raise Exception(
                    f"Method error. For {filter} in {file} you should use ['single'], \n while {subMethod} has been passed.")
        from fucs import filterFoldEnrichment as filterFunction

    elif filter == 'length':
        if method == None:
            method = [300, 500]
        elif method not in ['dist', 'polyfit'] and type(method) != list:
            raise Exception(
                f"Method error. For {filter} in {file} you should use one of ['dist','polyfit', [min, max]], \n while {method} has been passed.")
        from fucs import filterLength as filterFunction

    elif filter == 'summit':
        if method == None:
            method = 150
        elif type(method) != int:
            raise Exception(
                f"Method error. For {filter} in {file} you should use integer (+- int around summit), \n while {method} has been passed.")
        from fucs import evenLengthAroundSummit as filterFunction

    elif filter == 'likely':
        if type(method) not in [int, float]:
            raise Exception(
                f"Method error. For {filter} in {file} you should use a number, not {method}")
        from fucs import filterLikely as filterFunction

    elif filter == None:
        def filterFunction(df, method=method):
            return df

    return filterFunction, method
# singleFilter


def addTitle(filtered, title, filter, method):
    methodIsList = False
    if type(method) == list:
        methodIsList = True
        method = '_'.join([str(i) for i in method])
    title = f"{title}_{filter}_{method}"
    if filter == 'length' and method in ['dist', 'polyfit']:
        minLenght = int(filtered.length.min())
        maxLength = int(filtered.length.max())
        title = f"{title}_{minLenght}_{maxLength}"
    if filter == 'fold' and not methodIsList:
        threshFold = int(filtered.fold_enrichment.min())
        title = f"{title}_{threshFold}"
    return title
# addTitle


def pullPeakSequences(genomeFile, file, filter=None, method=None, multiFilter=False, multiFilterMethods=None):
    title = os.path.splitext(file.split('/')[-1])[0]
    peakDF = readPeak(file)
    print('*' * 100)
    print(file)

    if multiFilter != False:
        filtered = peakDF.copy()
        for i, filter in enumerate(multiFilter):
            filterFunction, method = singleFilter(
                filter=filter, method=multiFilterMethods[i])
            filtered = filterFunction(filtered, method=method)
    else:
        filterFunction, method = singleFilter(filter=filter, method=method)
        filtered = filterFunction(peakDF, method=method)
    if len(filtered.index) <= 10:
        print(f'Only {len(filtered.index)} peaks left, skip.')
        return

    peakPositions = {}
    for index in filtered.index:
        peakPositions[index] = list(
            filtered.loc[index, ['start', 'end']].astype(int))
        if type(peakPositions[index][0]) != int:
            raise Exception(
                f'some thing wrong with \n{index}\n{filtered.loc[index,:]}')

    genome = SeqIO.read(genomeFile, 'genbank')

    if multiFilter == False:
        title = addTitle(filtered, title, filter, method)
    else:
        for i, filter in enumerate(multiFilter):
            title = addTitle(filtered, title, filter, multiFilterMethods[i])

    print(f'Output file {title}_seqs.fasta')
    outputFile = f'/Users/durand.dc/Desktop/ChIP1839/peak_calling/{title}_seqs.fasta'
    if os.path.isfile(outputFile):
        print(f'The result file exists for {filter} | {method}, skip.')
    else:
        peakSeqs = getSeq(genome, peakPositions)
        SeqIO.write(peakSeqs,
                    f'/Users/durand.dc/Desktop/ChIP1839/peak_calling/{title}_seqs.fasta',
                    'fasta')
    outputTable = f'/Users/durand.dc/Desktop/ChIP1839/peak_calling/{title}_seqs.xlsx'
    print(f'Output table {title}_seqs.xlsx')
    if os.path.isfile(outputTable):
        print(f'The result file exists for {filter} | {method}, skip.')
    else:
        filtered.to_excel(outputTable)


filterMethod = {'fold': ['single', 'min', 'max'],
                'length': [[300, 500], 'dist', 'polyfit'],
                'summit': [150],
                'likely': [1, 100],
                'none': [None]}
fileFilterMethod = {'.xls': {'fold': [0], 'length': [0, 1, 2], 'summit': [0], 'none': [0]},
                    'common_peaks.tsv': {'fold': [1, 2], 'length': [0], 'summit': [0], 'none': [0]},
                    'uniq_peaks.tsv': {'fold': [0], 'length': [0], 'summit': [0], 'none': [0]},
                    'common.bed': {'length': [0], 'likely': [0], 'none': [0]},
                    'cond1.bed': {'length': [0], 'likely': [1], 'none': [0]},
                    'cond2.bed': {'length': [0], 'likely': [1], 'none': [0]},
                    }
files = [
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/ChIP-25h_keepDup_model/ChIP-25h_peaks.xls',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/ChIP-48h_keepDup_model/ChIP-48h_peaks.xls',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/25-48_common_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/25-48_25_uniq_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/25-48_48_uniq_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/48-25_common_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/48-25_25_uniq_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks/48-25_48_uniq_peaks.tsv',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c3.0_common.bed',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c3.0_cond1.bed',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c3.0_cond2.bed',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c5.0_common.bed',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c5.0_cond1.bed',
    '/Users/durand.dc/Desktop/ChIP1839/peak_calling/diff_peaks_macs2/25-48_c5.0_cond2.bed',
]

if __name__ == '__main__':
    genomeFile = '/Users/durand.dc/Documents/works/Resources/Genomes_Info/Streptomyces_coelicolor/M145.gb'
    '''All defaults'''
    # for file in files:
    #     for ext in fileFilterMethod:
    #         if file.endswith(ext):
    #             filters = fileFilterMethod[ext]
    #             for filter in filters:
    #                 methodFilter = fileFilterMethod[ext][filter]
    #                 methods = [filterMethod[filter][i] for i in methodFilter]
    #                 if filter == 'none':
    #                     filter = None
    #                 for method in methods:
    #                     pullPeakSequences(genomeFile, file, filter, method)

    '''specifics'''
    # pullPeakSequences(genomeFile, files[0], 'length', [180, 250])
    # pullPeakSequences(genomeFile, files[1], 'length', [200, 270])
    # pullPeakSequences(genomeFile, files[0], multiFilter=[
    #                   'fold', 'summit'], multiFilterMethods=[['single', 20], 150])
    # pullPeakSequences(genomeFile, files[1], multiFilter=[
    #                   'fold', 'summit'], multiFilterMethods=[['single', 20], 170])
    # for file in files[3:5] + files[6:8]:
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'fold', 'summit'], multiFilterMethods=[['single', 20], 170])
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'fold', 'summit'], multiFilterMethods=[['single', 5], 170])
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'fold', 'summit'], multiFilterMethods=[['single', 10], 170])
    # for file in [files[8], files[11]]:
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'length', 'likely'], multiFilterMethods=[[200, 400], 1])
    #     pullPeakSequences(genomeFile, file, 'length', [200, 400])
    # for file in files[9:11] + files[12:]:
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'length', 'likely'], multiFilterMethods=[[200, 400], 100])
    # for file in [files[2], files[5]]:
    #     pullPeakSequences(genomeFile, file, multiFilter=[
    #                       'fold', 'summit'], multiFilterMethods=[['min', 20], 170])


    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', help='genome file, fastta or genbank')
    parser.add_argument('-f', '--files', nargs="+")

    args = parser.parse_args()
    genome = args.genome
    files = args.files