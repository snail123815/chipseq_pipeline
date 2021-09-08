import numpy.polynomial.polynomial as poly
import numpy as np

def flatten(lst):
    from functools import reduce
    from operator import add
    flatten = lambda lst: [lst] if type(lst) != list else reduce(add, [flatten(ele) for ele in lst])
    flattened = flatten(lst)
    flattened = [str(ele) for ele in flattened]
    return flattened
# flatten

def findStartLine(csvFile):
    count = 0
    with open(csvFile, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'):
                if line != '':
                    break
            count += 1

    return count
# findStartLine


def turning_points(array):
    ''' turning_points(array) -> min_indices, max_indices
    Finds the turning points within an 1D array and returns the indices of the minimum and
    maximum turning points in two separate lists.
    https://stackoverflow.com/a/48360671/6823079
    answered Jan 20 2018 at 20:18
    by Timo
    '''
    idx_max, idx_min = [], []
    if (len(array) < 3):
        return idx_min, idx_max

    NEUTRAL, RISING, FALLING = range(3)

    def get_state(a, b):
        if a < b:
            return RISING
        if a > b:
            return FALLING
        return NEUTRAL

    ps = get_state(array[0], array[1])
    begin = 1
    for i in range(2, len(array)):
        s = get_state(array[i - 1], array[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING:
                    idx_max.append((begin + i - 1) // 2)
                else:
                    idx_min.append((begin + i - 1) // 2)
            begin = i
            ps = s
    return idx_min, idx_max
# turning_points


def readCoverage():
    import pickle
    from os.path import isfile

    coverageFile = '/Users/durand.dc/Desktop/ChIP1839/bamSorted/coverage.txt'
    coverageDataframeBz2 = '/Users/durand.dc/Desktop/ChIP1839/bamSorted/coverage.pickle.bz2'
    coverageDataframe = '/Users/durand.dc/Desktop/ChIP1839/bamSorted/coverage.pickle'

    if isfile(coverageDataframe):
        with open(coverageDataframe, 'rb') as dataFile:
            coverage = pickle.load(dataFile)
    elif isfile(coverageDataframeBz2):
        import bz2
        with bz2.open(coverageDataframe, 'rb') as dataFile:
            coverage = pickle.load(dataFile)
    else:
        import bz2
        import pandas as pd
        coverage = pd.read_csv(coverageFile, delimiter='\t', usecols=[
                               'pos', 'ChIP-48h', 'ChIP-25h', 'gDNA-25h', 'gDNA-48h'], index_col='pos')
        with bz2.open(coverageDataframeBz2, 'wb') as dataFile:
            pickle.dump(coverage, dataFile)
        with open(coverageDataframe, 'wb') as dataFile:
            pickle.dump(coverage, dataFile)

    return coverage
# readCoverage


def getSpanFetures(genome, startOri, endOri, expand=20000):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from math import ceil
    newStart = max(0, startOri-expand)
    newEnd = min(len(genome), endOri+expand)
    sourceSeq = genome[newStart:newEnd]

    start = startOri-newStart 
    end = start + (endOri - startOri)

    spanFeats = []
    for feat in sourceSeq.features:
        spanStart = start in feat
        spanEnd = end in feat
        # feat.__contains__(self, value)
            #Check if an integer position is within the feature.

        spanFeat = SeqFeature(type=feat.type)

        if spanStart and spanEnd:
            # Target position is inside feature
            if feat.type == 'CDS':
                # calculate correct start and end location to make it inframe
                newStart = (3 - abs(start - feat.location.start)) % 3
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newStart = 0
                newEnd = end - start
            spanFeat.location = FeatureLocation(newStart,
                                                newEnd,
                                                strand=feat.location.strand)
            for key in feat.qualifiers:
                if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
                    spanFeat.qualifiers[key] = [
                        f'{keyStr} (slice)' for keyStr in feat.qualifiers[key]]
                elif key == 'translation':
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPointA = ceil((start - feat.location.start) / 3)
                        cutPointB = (end - feat.location.start) // 3
                    else:
                        len(spanFeat.qualifiers[key][0])
                        cutPointA = len(
                            spanFeat.qualifiers[key][0]) + 1 - (end - feat.location.start) // 3
                        cutPointB = len(
                            spanFeat.qualifiers[key][0]) + 1 - ceil((start - feat.location.start) / 3)
                    spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPointA:cutPointB]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanStart:
            # Start position inside feature, feature ends in this range
            if feat.type == 'CDS':
                newStart = (3 - abs(start - feat.location.start)) % 3
            else:
                newStart = 0
            newEnd = feat.location.end - start
            spanFeat.location = FeatureLocation(newStart,
                                                newEnd,
                                                strand=feat.location.strand)
            for key in feat.qualifiers:
                if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
                    spanFeat.qualifiers[key] = [
                        f'{keyStr} (right part)' for keyStr in feat.qualifiers[key]]
                elif key == 'translation':
                    spanFeat.qualifiers[key] = [
                        keyStr for keyStr in feat.qualifiers[key]]
                    if feat.location.strand == 1:
                        cutPoint = len(
                            spanFeat.qualifiers[key][0]) - ceil(len(spanFeat) / 3) + 1
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPoint:]
                    else:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][:cutPoint]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanEnd:
            # End position inside feature, feature ends in this range
            if feat.type == 'CDS':
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newEnd = end - start
            newStart = feat.location.start - start
            spanFeat.location = FeatureLocation(newStart,
                                                newEnd,
                                                strand=feat.location.strand)
            for key in feat.qualifiers:
                if key in ['gene_synonym', 'locus_tag', 'product', 'protein_id', 'db_xref', 'mol_type']:
                    spanFeat.qualifiers[key] = [
                        f'{keyStr} (left part)' for keyStr in feat.qualifiers[key]]
                elif key == 'translation':
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][:cutPoint]
                    else:
                        cutPoint = ceil((len(feat) - len(spanFeat)) / 3)
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][cutPoint:]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        else:
            # Not in range, ignore
            continue

    return spanFeats
# getSpanFetures


def peakLengthAtLoc(data, series):
    if not hasattr(series, "__iter__"):
        series = [series, ]
    lengthNums = []
    for x in series:
        if 0 <= x <= len(data.index):
            lengthNums.append(data.length.iloc[np.int(x)])
        else:
            lengthNums.append('')
    if len(lengthNums) == 1:
        return lengthNums[0]
    return lengthNums
# peakLengthAtLoc


def minMaxPolyfit(x, y, degree=10):
    """x,y,degree = 10
    x,y are list of values"""
    coefs = poly.polyfit(x, y, degree)
    Fit = poly.Polynomial(coefs)
    yTrend = [np.NaN, ]
    yTrend_trend = [np.NaN, ]
    for index in range(len(y) - 1):
        if index < 500:
            yTrend.append(np.NaN)
            yTrend_trend.append(np.NaN)
        else:
            yTrend.append(Fit(index + 1) - Fit(index))
            if np.isnan(yTrend[index - 1]):
                yTrend_trend.append(np.NaN)
            else:
                # this is to find real turning point - tangent angle change
                yTrend_trend.append(yTrend[index] - yTrend[index - 1])
    min_indices, max_indices = turning_points(yTrend_trend)
    return min_indices, max_indices, Fit
# minMaxPolyfit


def filterFoldEnrichment(peakDF, method=['single', 10]):
    methodIsList = False
    if type(method) == list:
        methodIsList = True
        thresh = method[1]
        method = method[0]
    else:
        thresh = 10

    if method not in ['single', 'max', 'min']:
        raise Exception(f"Wrong method '{method}' ('single','max','min')")
    if method == 'single':
        if 'fold_enrichment' not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}")
        y = peakDF.fold_enrichment
        if len(y) <= 2000 or methodIsList:
            print(
                f'No polyfit, setting thresh for fold_enrichment = {thresh}')
        else:
            bondarieMin = 1300
            x = range(1, len(y) + 1)
            min_indices, max_indices, Fit = minMaxPolyfit(x, y)
            thresh = min([Fit(lim) for lim in [
                         indice for indice in max_indices if indice >= bondarieMin]])
            print(f'Done polyfit, thresh for fold_enrichment = {thresh}')
    if method == 'max':
        if 'fold_enrichment_A' not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}")
        maxFold = peakDF.loc[:, ['fold_enrichment_A',
                                 'fold_enrichment_B']].max(axis=1)
        peakDF = peakDF.assign(fold_enrichment=maxFold)
    if method == 'min':
        if 'fold_enrichment_A' not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}")
        minFold = peakDF.loc[:, ['fold_enrichment_A',
                                 'fold_enrichment_B']].min(axis=1)
        peakDF = peakDF.assign(fold_enrichment=minFold)
    filtered = peakDF.query('fold_enrichment >= @thresh')
    return filtered
# filterFoldEnrichment method = ['single','max','min']


def filterLength(peakDF, method='dist'):
    peakDF = peakDF.assign(length=peakDF.end - peakDF.start)
    peakDF.sort_values('length', inplace=True)
    if method == 'dist':
        medians = peakDF.length.median()
        minLength = int(peakDF.length.quantile(q=[0.25]))
        maxLength = int(peakDF.length.quantile(q=[0.75]))
    elif method == 'polyfit':
        # data has to be origional peak files
        # this has to be checked by plotting in peak_enrichment_modeling
        y = peakDF.fold_enrichment
        if len(y) <= 2000:
            raise Exception(
                'You are doing polyfit in a data set less then 2000 entries.')
        bondarieMin = 1300
        x = range(1, len(y) + 1)
        min_indices, max_indices, Fit = minMaxPolyfit(x, y)
        foldThresh = min([Fit(lim) for lim in [
                         indice for indice in max_indices if indice >= bondarieMin]])
        max_indices = [
            indice for indice in max_indices if indice >= bondarieMin]
        if len(max_indices) != 2:
            raise Exception(
                f'Please check the ployfit graph, there are {len(max_indices)} turning points for now.')
        # gether data for output
        minLength, maxLength = [peakLengthAtLoc(peakDF, indice)
                                for indice in max_indices if indice >= bondarieMin]
    elif type(method) == list:
        try:
            minLength, maxLength = method
        except:
            raise Exception(
                'Method should be [minLength, maxLength] if you do not fillter')
    else:
        raise Exception(f'Method {method} wrong for length filter')
    filtered = peakDF.query('@minLength <= length <= @maxLength')
    filtered.sort_index()
    return filtered
# filterLength method = ['dist','polyfit']


def evenLengthAroundSummit(peakDF, method):
    rangeSummit = method
    if 'abs_summit' not in peakDF.columns:
        raise Exception('No summit info in dataframe')
    peakDF.start = peakDF.abs_summit - rangeSummit
    peakDF.end = peakDF.abs_summit + rangeSummit
    return peakDF
# evenLengthAroundSummit


def filterLikely(peakDF, method):
    thresh = float(method)
    if 'likely_difference' in peakDF.columns:
        if thresh >= 10:
            raise Exception(
                'Threshold error, for common peaks the value should be the lesser the better, set up a small value ~1')
        filtered = peakDF.query('likely_difference <= @thresh')
    elif 'log10_likely' in peakDF.columns:
        if thresh < 10:
            raise Exception(
                'Threshold error, for unique peaks, the value should be the bigger the better, set up a big value >=10')
        filtered = peakDF.query('log10_likely >= @thresh')
    elif '-log10(pvalue)' in peakDF.columns:
        if thresh < 10:
            raise Exception(
                'Threshold error, for peaks, the value should be the bigger the better, set up a big value >=10')
        filtered = peakDF.query('`-log10(pvalue)` >= @thresh')
        
    else:
        raise Exception(
            f'Likely filter failed for this data frame with header\n{peakDF.columns}')
    return filtered
# filterLikely

