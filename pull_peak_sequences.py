from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import concurrent.futures 
import sys  # for error handling

from funcs import getSpanFetures 
import pandas as pd


def readPeak(file, thresh=0):
    ''' Three kind of dataframe generated
    name, start, end, abs_summit, fold_enrichment
    name, start, end, abs_summit, fold_enrichment_A, fold_enrichment_B
    name, start, end, log10_likely (of the peak presented in this file)
    thresh is the threshold of the -log10(pvalue)
    '''
    if file.endswith('.xls'):
        # direct peak calling result
        data = pd.read_csv(
            file, delimiter='\t',
            comment='#', index_col='name',
            usecols=[
                'name', 'start', 'end', 'abs_summit','-log10(pvalue)', 'fold_enrichment'
                ]
            )
        data = data[data['-log10(pvalue)']>=thresh]
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
            cols = ['name', 'start', 'end', 'abs_summit','-log10(pvalue)', 'fold_enrichment']
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
            futures.append(executor.submit(slice, sourceSeq, peakDict[peak], peak))
        for future in concurrent.futures.as_completed(futures):
            sliceSeq = future.result()
            resultSeqs.append(sliceSeq)
    return resultSeqs
# getSeq


def singleFilter(filter, method, filename):
    if filter not in ['fold', 'length', 'summit', 'likely', None]:
        raise Exception(
            f"Filter error, {filter} ['fold', 'length', 'summit', 'likely', None]")
    print(f'Filtering: {filter} | {method}')
    if filter == 'fold':
        # filename ext in ['.tsv', '.xls']
        if type(method) == list:
            # make sure the method returned by singleFilter() still a list
            subMethod = method[0]
        else:
            subMethod = method
        if 'common' in filename:
            if subMethod not in ['max', 'min']:
                raise Exception(
                    f"Method error. For {filter} in {filename} you should use ['max', 'min'], \n while {subMethod} has been passed.")
        else:
            if subMethod != 'single':
                raise Exception(
                    f"Method error. For {filter} in {filename} you should use ['single'], \n while {subMethod} has been passed.")
        from funcs import filterFoldEnrichment as filterFunction

    elif filter == 'length':
        if method == None:
            method = [300, 500]
        elif method not in ['dist', 'polyfit'] and type(method) != list:
            raise Exception(
                f"Method error. For {filter} in {filename} you should use one of ['dist','polyfit', [min, max]], \n while {method} has been passed.")
        from funcs import filterLength as filterFunction

    elif filter == 'summit':
        if method == None:
            method = 150
        elif type(method) != int:
            raise Exception(
                f"Method error. For {filter} in {filename} you should use integer (+- int around summit), \n while {method} has been passed.")
        from funcs import evenLengthAroundSummit as filterFunction

    elif filter == 'likely':
        try:
            method = float(method)
        except:
            raise Exception(
                f"Method error. For {filter} in {filename} you should use a number, not {method}")
        from funcs import filterLikely as filterFunction

    else:
        def filterFunction(df, method=None):
            return df

    return filterFunction, method
# singleFilter


def addTitle(filtered, title, filter, method):
    methodIsList = False
    if type(method) == list:
        methodIsList = True
        method = '_'.join([str(i) for i in method])
        title = f"{title}_{filter}_{method}"
    elif filter == 'length' and method in ['dist', 'polyfit']:
        minLenght = int(filtered.length.min())
        maxLength = int(filtered.length.max())
        title = f"{title}_{minLenght}_{maxLength}"
    elif filter == 'fold' and not methodIsList:
        threshFold = int(filtered.fold_enrichment.min())
        title = f"{title}_{threshFold}"
    elif filter == 'likely':
        title = f'{title}_{method}'
    else:
        raise NameError(f'Non-supported filter {filter}, method {method}')
    return title
# addTitle

# Main function:
def pullPeakSequences(genomeFile, file, output,
                      filter=None, method=None, 
                      multiFilter=False, multiFilterMethods=None):

    title = os.path.splitext(file.split('/')[-1])[0]
    peakDF = readPeak(file)
    print('*' * 100)
    print(file)

    if multiFilter:
        filtered = peakDF.copy()
        for i, filter in enumerate(multiFilter):
            filterFunction, method = singleFilter(
                filter=filter, method=multiFilterMethods[i], filename=file)
            filtered = filterFunction(filtered, method=method)
    else:
        filterFunction, method = singleFilter(filter=filter, method=method, filename=file)
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

    ext = os.path.splitext(genomeFile)[1]
    if ext[1:] in ['gb', 'genbank', 'gff']:
        genome = SeqIO.read(genomeFile, 'genbank')
    elif ext[1:] in ['fa', 'fasta', 'faa']:
        genome = SeqIO.read(genomeFile, 'fasta')
    else:
        raise NameError(genomeFile, ext)

    if multiFilter == False:
        title = addTitle(filtered, title, filter, method)
    else:
        for i, filter in enumerate(multiFilter):
            title = addTitle(filtered, title, filter, multiFilterMethods[i])

    print(f'Output file {title}_seqs.fasta')
    if os.path.isfile(output):
        print(f'The result file exists for {filter} | {method}, skip.')
    else:
        peakSeqs = getSeq(genome, peakPositions)
        SeqIO.write(peakSeqs,
                    output,
                    'fasta')
    outputTable = f'{output}.xlsx'
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
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', help='genome file, fastta or genbank')
    parser.add_argument('-f', '--files', nargs="+")
    parser.add_argument('--filter', help='''
    filter method, use one of the keys of dict
    filterMethod = {'fold': ['single', 'min', 'max'],
                    'length': [[300, 500], 'dist', 'polyfit'],
                    'summit': [150],
                    'likely': [1, 100],
                    'none': [None]}
     ''')
    parser.add_argument('-o', help='output path')
    parser.add_argument('--likely', help='threshold for filter "likely"')

    args = parser.parse_args()
    genome = args.genome
    files = args.files
    fm = args.filter
    fmlikely = args.likely
    outputPath = args.o

    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)
    
    for file in files:
        output = os.path.join(outputPath, f'{os.path.splitext(os.path.split(file)[1])[0]}.fa')
        pullPeakSequences(genome, file, output, filter=fm, method=fmlikely)
