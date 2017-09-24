######################################################################################################################################################################
#       Copyright 2017, Sebastian Pechmann, Universite de Montreal (The MIT License)
#
#       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
#       to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#       and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#       FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#       WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################################################################################################################

import sys, io
import numpy as np 
import scipy.stats

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects.numpy2ri import numpy2ri



tables = np.load('data/reference/tables.npz')
synonymous = tables['synonymous']
conservative = tables['conservative']
neighbours = tables['neighbours']

aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

codons = [ 'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
           'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
           'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
           'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT' ]

gencode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
           "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P", "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

stopidx = [48, 50, 56]

aa_hydrophobicity = {"I": 4.5,  "V":  4.2, "L":  3.8, "F":  2.8, "C":  2.5, "M":  1.9, "A":  1.8, "G": -0.4, "T": -0.7, "S": -0.8,
                     "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5 }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def loadBatch(replica):
    """
    load all passages for wt and ga
    """

    cc_wt = np.zeros(( 2218, 64, 7 ))
    cc_ga = np.zeros(( 2218, 64, 7 ))


    for i in range(7):
        no = int(i) + 2

        name_wt = 'data/counts/' + replica + '.wt.p' + str(no) + '.txt'
        name_ga = 'data/counts/' + replica + '.ga.p' + str(no) + '.txt'

        cc_wt[:,:,i] = np.genfromtxt(name_wt)
        cc_ga[:,:,i] = np.genfromtxt(name_ga)


    return cc_wt, cc_ga



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_average_coverage(wt, ga):

    wt = np.sum(wt, 1)
    ga = np.sum(ga, 1)

    data_wt = numpy2ri(wt)
    data_ga = numpy2ri(ga)

    r.assign('wt', data_wt)
    r.assign('ga', data_ga)

    r(' wt <- as.matrix(wt) ')
    r(' ga <- as.matrix(ga) ')
    r(' source("src/R/figure_coverage.R") ')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mutation_rate(data):
    """
    1. simple mutations sum(mutations)/sum(coverage)
    2. stop codon mutation rate: mu = 3N / M
    N = number of stop codons
    M = number of NSMT (non-sense mutation targets)
    """

    seq = loadRefSeq('data/reference/polio.codonseq')   

    NSMT = np.array([ 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                      1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                      0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0 ])
 
    rate = np.zeros((2200))
    stoprate = [] 

    N = 0
    M = 0

    for pos in range(2200): #np.shape(data[0]):

        coverage  = 0
        mutations = 1

        current_codon = seq[pos]
        current_idx = codons.index(current_codon)

        current_mutations = sum(data[pos,stopidx] ) 
        current_neighbours = neighbours[current_idx,:]
        current_stop = sum(current_neighbours[stopidx])
        current_stoptarget = sum(data[pos,:]) * float(current_stop) / 9.

        coverage += current_stoptarget
        mutations += current_mutations

        current_N = current_stop
        current_M = sum(data[pos,:] *  NSMT )

        if current_N > 0:
            N += current_N
            M += current_M

        if current_stop > 0 and coverage > 0:
            rate[pos] = mutations / float(coverage)

        if current_N > 0 and ~np.isnan(current_N) and current_M > 0 and ~np.isnan(current_M):
            current_stoprate = 3*float(current_N)/current_M
            stoprate.append(current_stoprate)

    rate = np.median(rate)

    stoprate = np.array(stoprate)
    stoprate = np.mean(stoprate)

    return stoprate 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_mutrate(wt, ga):

    ratemat = np.zeros(( 7, 2 ))

    for i in range(7):
        ratemat[i,0] = mutation_rate(wt[:,:,i])
        ratemat[i,1] = mutation_rate(ga[:,:,i])

    d = numpy2ri(ratemat)
    r.assign('data', d)

    r(' source("src/R/figure_mutrate.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dnds_count(data):
    """
    compute evolutionary rates with Nei-Gojobori model 
    """

    sites = getsites_virus(data)
    allsynonymous = np.zeros((64))
    allnonsynonymous = np.zeros((64))

    for i in range(64):
        current_neighbours = neighbours[:,i] == 1
        current_synonymous = synonymous[:,i] == 1
        current_nonsynonymous = synonymous[:,i] == 0

        allsynonymous[i] = np.sum( current_neighbours * current_synonymous) / 3.
        allnonsynonymous[i] = np.sum( current_neighbours * current_nonsynonymous) / 3.

    allsites = ( np.sum(sites, axis=1) + np.sum(sites, axis=0) ) / 2.
    all_s = allsynonymous * allsites
    all_n = allnonsynonymous * allsites

    S   = np.sum(all_s)
    N   = np.sum(all_n)

    # get observed mutations
    all_sd = 0
    all_nd = 0

    for a in range(64):
        for b in range(64):
            if a not in stopidx and b not in stopidx:
                mutcount = sites[a,b]
                [diff_s, diff_n] = pathwaysdiff(codons[int(a)], codons[int(b)], codons, synonymous)
                all_sd += (diff_s * mutcount)
                all_nd += (diff_n * mutcount)
    Sd  = np.sum(all_sd)
    Nd  = np.sum(all_nd)

    # calculate rates
    if S == 0:
        S = 1
    if N == 0:
        N = 1

    pS  = Sd / float(S)
    pN  = Nd / float(N)

    if (pS < (3./4.)):
        ds = -(3./4.) * np.log(1-((4./3.)*pS))
    else:
        ds = NaN

    if (pN < (3./4.)):
        dn = -(3./4.) * np.log(1-((4./3.)*pN))
    else:
        dn = NaN

    if ds != 0:
        omega = dn/ds
    else:
        omega = NaN

    return round(ds, 6), round(dn, 6)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getsites_virus(codon_counts):
    stopindx = [48, 50, 56]
    sites = np.zeros(( 64, 64 ))
    for i in range(codon_counts.shape[0]):
        w = np.argmax(codon_counts[i,:])
        if sum(codon_counts[i,:]) > 0:
            sites[w,:] +=  (codon_counts[i,:] /  sum(codon_counts[i,:]))
    sites[stopindx,:] = 0
    sites[:,stopindx] = 0

    return np.around(sites, 6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pathwaysdiff(codon1, codon2, codons, synonymous):
    import itertools
    sd, nd = (0,)*2

    # Get all possible evolutionary paths between two codons
    codon1 = np.array(list(codon1))
    codon2 = np.array(list(codon2))
    c1 = ''.join(codon1)

    diffs = np.where(codon1 != codon2)                     # nt positiosn at which codons differ
    diffs = diffs[0].tolist()
    order = np.array(list( itertools.permutations(diffs) ) )  # all possible permutations
    (paths, changes) = order.shape
    pathways = np.zeros((paths, changes+1))
    pathways[:,0] = codons.index(c1)

    for i in range(paths):
        codon = np.copy(codon1)
        for j in range(changes):
            mutation = order[i,j]
            codon[mutation] = codon2[mutation]
            c =  ''.join(codon)
            ci = codons.index(c)
            pathways[i,j+1] = ci

    # Remove paths with stop codons
    if paths > 1:
        pathways = pathways[~(pathways==48).any(1), :]
        pathways = pathways[~(pathways==50).any(1), :]
        pathways = pathways[~(pathways==56).any(1), :]
    paths2 = pathways.shape[0]

    # Classify all mutations for all paths
    syn_p, nsyn_p = (float(0),)*2

    for i in range(paths2):
        for j in range(changes):
            first  = int(pathways[i,j])
            second = int(pathways[i,j+1])
            syn = synonymous[first, second]
            if (syn == 1):
                syn_p += 1
            else:
                nsyn_p += 1
    dff = np.array([syn_p, nsyn_p])
    dff = dff /float(paths2)
    dff = np.around(dff, 5)
    return dff


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_dnds(data):
    """
    compute dnds for all passages
    """

    dnds = np.zeros(( 7, 2 ))
    for i in range(7):
        dnds[i,:] = dnds_count(data[:,:,i])

    return dnds

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_dnds(wt, ga):

    dnds_wt = numpy2ri(wt)
    dnds_ga = numpy2ri(ga)

    r.assign('wt', dnds_wt)
    r.assign('ga', dnds_ga)

    r(' source("src/R/figure_dnds.R") ')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mutationalEntropy(datacc, pos, seq, AA):

    H = 0
    theta = 1

    for i in pos:
        current_codon = seq[i]
        current_aa = gencode[current_codon]

        current_idx = codons.index(current_codon)
        current_counts = np.copy(datacc[i,:])
        current_counts[current_idx] = 0

        for aa in aminoacids:
            current_aa = np.array(AA) == aa

            sel = np.sum( current_counts[current_aa] )
            all = np.sum(current_counts)

            if sel > 0 and all > theta:

                P = sel / all
                H += P * np.log(P)

    return -H



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def loadRefSeq(fname):
    """
    load polio reference sequence
    """

    fileIN = io.open(fname, 'r')
    line = fileIN.readline()
    seq = []
    while line:
        current_line = line.rstrip('\n\r')
        seq.append(str( current_line) )
        line = fileIN.readline()

    return seq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_entropy(data):
    """
    compute entropy of mutation space for all passages
    """

    seq = loadRefSeq('data/reference/polio.codonseq')

    AAcode = []
    for i in range(64):
        AAcode.append(gencode[codons[i]])

    entropy = np.zeros(( 2200, 7 ))

    for i in range(7):
        current_data = data[:,:,i] 
        for j in range(2200): 
            entropy[j, i] = mutationalEntropy(current_data, np.array([j]), seq, AAcode)

    return entropy


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_entropy(wt, ga):

    ent_wt = numpy2ri(wt)
    ent_ga = numpy2ri(ga)

    r.assign('wt', ent_wt)
    r.assign('ga', ent_ga)

    r(' source("src/R/figure_entropy.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mutationalDiversity(data, pos, seq, AA):

    D = 0
    theta = 1

    for i in pos:
        current_codon = seq[i]
        current_aa = gencode[current_codon]

        current_idx = codons.index(current_codon)
        current_counts = np.copy(data[i,:])
        current_counts[current_idx] = 0

        for aa in aminoacids:
            current_aa = np.array(AA) == aa

            n = np.sum( current_counts[current_aa] )
            N = np.sum(current_counts)

            if n > 0 and N > theta:

                D += np.square(float(n) / float(N))

    return D

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_diversity(data):
    """
    compute entropy of mutation space for all passages
    """

    seq = loadRefSeq('data/reference/polio.codonseq')

    AAcode = []
    for i in range(64):
        AAcode.append(gencode[codons[i]])

    diversity = np.zeros(( 2200, 7 ))

    for i in range(7):
        for j in range(2200):
            diversity[j, i] = mutationalDiversity(data[:,:,i], np.array([j]), seq, AAcode)

    return diversity


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_diversity(wt, ga):

    div_wt = numpy2ri(wt)
    div_ga = numpy2ri(ga)

    r.assign('wt', div_wt)
    r.assign('ga', div_ga)

    r(' source("src/R/figure_diversity.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hela_ws = {"AAA":0.4138, "AAC":0.8276, "AAG":0.5462, "AAT":0.3633, "ACA":0.2759, "ACC":0.2731, "ACG":0.1572, "ACT":0.3793, 
           "AGA":0.2069, "AGC":0.2759, "AGG":0.2386, "AGT":0.1211, "ATA":0.1725, "ATC":0.3917, "ATG":0.3448, "ATT":0.4786,
           "CAA":0.3448, "CAC":0.3448, "CAG":0.4552, "CAT":0.1514, "CCA":0.2414, "CCC":0.2483, "CCG":0.2152, "CCT":0.3448, 
           "CGA":0.2069, "CGC":0.1738, "CGG":0.2386, "CGT":0.2414, "CTA":0.069, "CTC":0.2483, "CTG":0.3324, "CTT":0.3448,
           "GAA":0.1379, "GAC":0.5172, "GAG":0.5959, "GAT":0.2271, "GCA":0.3794, "GCC":0.4469, "GCG":0.2248, "GCT":0.6207, 
           "GGA":0.3103, "GGC":0.4138, "GGG":0.3062, "GGT":0.1817, "GTA":0.1725, "GTC":0.2731, "GTG":0.6069, "GTT":0.3793,
           "TAC":0.4483, "TAT":0.1968, "TCA":0.0345, "TCC":0.3476, "TCG":0.149, "TCT":0.4828, "TGC":1, "TGG":0.2414, 
           "TGT":0.439, "TTA":0.1724, "TTC":0.3793, "TTG":0.2621, "TTT":0.1665}

human_ws = {"AAA": 0.6077, "AAC": 0.9867, "AAG": 0.7746, "AAT": 0.4521, "ACA": 0.1658, "ACC": 0.1989, "ACG": 0.2188, "ACT": 0.2762, 
            "AGA": 0.1657, "AGC": 0.2409, "AGG": 0.3293, "AGT": 0.1246, "ATA": 0.1382, "ATC": 0.5591, "ATG": 0.5525, "ATT": 0.5666,
            "CAA": 0.7459, "CAC": 0.3039, "CAG": 0.9845, "CAT": 0.1334, "CCA": 0.2763, "CCC": 0.2862, "CCG": 0.1989, "CCT": 0.3712, 
            "CGA": 0.1658, "CGC": 0.1392, "CGG": 0.1912, "CGT": 0.1934, "CTA": 0.1105, "CTC": 0.2586, "CTG": 0.3116, "CTT": 0.3591,
            "GAA": 0.5249, "GAC": 0.5525, "GAG": 0.8309, "GAT": 0.2425, "GCA": 0.304, "GCC": 0.5967, "GCG": 0.2354, "GCT": 0.8287, 
            "GGA": 0.3039, "GGC": 0.4144, "GGG": 0.3459, "GGT": 0.1819, "GTA": 0.1658, "GTC": 0.2188, "GTG": 0.6055, "GTT": 0.3039,
            "TAC": 0.4619, "TAT": 0.2217, "TCA": 0.1934, "TCC": 0.2387, "TCG": 0.1724, "TCT": 0.3315, "TGC": 0.8762, "TGG": 1, 
            "TGT": 0.4036, "TTA": 0.3039, "TTC": 0.4696, "TTG": 0.3182, "TTT": 0.2062}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tai_mutations(site, refseq, te, cc):
    """
    compute tAI for mutational space
    """

    current_codon = refseq[site]
    current_idx = codons.index(current_codon)
    current_aa = gencode[current_codon]
    current_neighbours = neighbours[current_idx,:]
    current_synonymous = synonymous[current_idx,:]
    current_te = te.get(current_codon, np.nan)
    current_count = int(cc[site,current_idx])

    tai = 0
    counter = 0

    syn = True 
    if syn:
        neighbourhood = current_neighbours * current_synonymous
    else:
        neighbourhood = np.copy(current_neighbours)

    neighbourhood = np.where( neighbourhood == 1 )[0]

    rnei = []
    for i in neighbourhood:
        if i != current_idx:
            mut_codon = codons[i]
            mut_te = float( te.get(mut_codon, np.nan) )
            mut_count = cc[site,i]

            if mut_te != np.nan:
                tai += (mut_te) * mut_count
                counter += mut_count
                rnei.append(mut_te)

    if (len(rnei) > 1) and counter > 0:
        tai /= counter
        tai = round(float(tai), 4)
    else:
        tai = np.nan

    if len(rnei) > 1:
        rtai = np.mean(np.array(rnei))
    else:
        rtai = np.nan

    return current_te, rtai, tai, current_count, int(counter)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_mutationtAI(wt, ga, scale):
    """
    delta-tAI of mutation space
    """

    if scale == "hela":
        te = hela_ws
    elif scale == "human":
        te = human_ws

    seq = loadRefSeq('data/reference/polio.codonseq')
    muttai = np.zeros(( 2200, 8, 7))

    for passage in range( 7 ):
        for i in range( 2200 ):
            te_seq, te_r, te_wt, count_seq_wt, cn_wt = tai_mutations(i, seq, te, wt[:,:,passage])
            te_seq, te_r, te_ga, count_seq_ga, cn_ga = tai_mutations(i, seq, te, ga[:,:,passage])
            muttai[i,:,passage] = np.array([te_seq, te_r, te_wt, te_ga, count_seq_wt, cn_wt, count_seq_ga, cn_ga])

    return muttai


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_muttai(data):

    d = numpy2ri(data)
    r.assign('data', d)

    r(' source("src/R/figure_tAI.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def MHtest(data):
    """
    perform MH test in R (faster)
    """ 

    data = numpy2ri(data)
    r.assign('D', data)
    r(' result <- mantelhaen.test(D) ')
    r(' pval <- result$p.value ')
    p = np.array(r.pval)
    r(' odds <- result$estimate ')
    OR = np.array(r.odds)

    return round(p[0], 10), OR[0]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Ftest(data):
    """
    perfor Fisher's exact test in R (faster)
    """

    p = 1
    OR = 1

    if amin( data ) > 0:
        data = numpy2ri(data)
        r.assign('D', data)
        r(' result <- fisher.test(D) ')
        r(' pval <- result$p.value ')
        p = np.array(r.pval)
        r(' odds <- result$estimate ')
        OR = np.array(r.odds)

        p  = p[0]
        OR = OR[0]

    return round(p, 10), round(OR, 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def test_allpassages(site, refseq, cc_wt, cc_ga, mutmat, syn):

    pssgs = range(7)
    n = len(pssgs)
 
    mhmutations = np.zeros(( 2, 2, n ))
    mhcounts = np.zeros(( 2, 2, n ))
    s1 = np.zeros(( n ))
    s2 = np.zeros(( n ))

    for ix, i in enumerate(pssgs):
        if np.size(site) > 1:
            for j in site:
                current_mutations = get_mutclass(j, i, refseq, cc_wt, cc_ga, mutmat, syn)
                if np.amin(current_mutations) > -1:
                    mhmutations[:,:,ix] += current_mutations
        else:
            current_mutations = get_mutclass(int(site), i, refseq, cc_wt, cc_ga, mutmat, syn)
            if np.amin(current_mutations) > 0:
                mhmutations[:,:,ix] = current_mutations

        s1[ix] = np.amin(mhmutations[:,:,ix] )

    nonzero1 = s1 != 0
    mhmutations = mhmutations[:,:,nonzero1]

    if sum(nonzero1) > 4: 	# require >0 mutations in at least 5 out of 7 passages
        out1   = MHtest(mhmutations)
        pval_mut = out1[0]
        OR_mut = out1[1]
    else:
        pval_mut = 1
        OR_mut = 1

    return OR_mut, pval_mut 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_mutclass(pos, passage, refseq, cc_wt, cc_ga, mutmat, syn):

    current_wt_codon = refseq[pos]
    current_wt_idx = codons.index(current_wt_codon)

    current_neighbours = neighbours[current_wt_idx,:]
    current_synonymous = synonymous[current_wt_idx,:]

    sel_class = mutmat[pos,:] == 1 
    sel_nclass = mutmat[pos,:] == 0 

    if syn:
        sc = np.array(sel_class, dtype=bool) * np.array(current_neighbours, dtype=bool) * np.array(current_synonymous, dtype=bool)
        sn = np.array(sel_nclass, dtype=bool) * np.array(current_neighbours, dtype=bool) * np.array(current_synonymous, dtype=bool)
    else:
        sc = np.array(sel_class, dtype=bool) * np.array(current_neighbours, dtype=bool) 
        sn = np.array(sel_nclass, dtype=bool) * np.array(current_neighbours, dtype=bool) 

    wt_class    = np.sum( cc_wt[pos, sc, passage] )
    wt_nonclass = np.sum( cc_wt[pos, sn, passage] )
    ga_class    = np.sum( cc_ga[pos, sc, passage] )
    ga_nonclass = np.sum( cc_ga[pos, sn, passage] )

    a = np.array([[wt_class, ga_class], [wt_nonclass, ga_nonclass]])
    a = np.around(a, 0)

    minmut = True
    if minmut:
        if np.amin(a) < 1:
            a = np.zeros(( 2,2 ))

    return a




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_te(wt, ga, type):
    """
    test for differences in mutations to slow codons
    """
    
    diff = np.zeros(( 2200, 2 )) 
    seq = loadRefSeq('data/reference/polio.codonseq')

    hela_tAI =   ['TCA','CTA','AGT','GAA','TCG','CAT','ACG','TTT','TTA','ATA','GTA','CGC'] 	# bottom 20%
    human_tAI =  ['CTA','AGT','CAT','ATA','CGC','AGA','ACA','CGA','GTA','TCG','GGT','CGG'] 	# bottom 20% 

    if type == "hela":
        slowcodons = hela_tAI[:]									
    elif type == "human":
        slowcodons = human_tAI[:]


    mutmat = np.zeros(( 64 ))
    for ix, i in enumerate(codons):
        if i in slowcodons:
            mutmat[ix] = 1
    mutmat = np.tile(mutmat, (2220, 1) )

    for i in range(2200):
        OR, p = test_allpassages(i, seq, wt, ga, mutmat, True)
        diff[i,] = np.array([p, OR])

    return diff




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_tediff(data):

    d = numpy2ri(data)
    r.assign('data', d)

    r(' source("src/R/figure_nopt.R") ')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_agg(wt, ga):
    """
    test for sites that differ in aggregation prone mutations
    """

    # initiate
    diff = np.zeros(( 2200, 2 ))
    seq = loadRefSeq('data/reference/polio.codonseq')
    theta = 2	

    # prepare mutation matrix
    mutmat = np.zeros(( 2220, 64 ))

    fileIN = io.open('data/precomputed/polio.AP.txt', 'r')
    line = fileIN.readline()
    while line:
        current_line = line.split()
        pos = int(current_line[0])
        codon = current_line[1]

        if seq[int(pos)-1] == codon and len(current_line) > 4:
            mut = current_line[3]
            tango = round(float(current_line[4]), 2)
            if tango < -1000:
                tango = 0

            for ix, i in enumerate(codons):
                current_aa = gencode[i]
                if current_aa == mut:

                    mutmat[pos-1, ix] = tango

        line = fileIN.readline()

    mm2 = mutmat > theta
    mutmat = np.array(mm2, dtype=bool)

    # compute differences
    for i in range(2200):
        OR, p = test_allpassages(i, seq, wt, ga, mutmat, False)
        diff[i,] = np.array([p, OR])

    return diff


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_aggdiff(data):

    d = numpy2ri(data)
    r.assign('data', d)

    r(' source("src/R/figure_sites_agg.R") ')
    r(' source("src/R/figure_P1_agg.R") ')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_stab(wt, ga):
    """
    test for sites that differ in destabilizing mutations
    """

    # initiate
    diff = np.zeros(( 2200, 2 ))
    seq = loadRefSeq('data/reference/polio.codonseq')
    theta = 1           

    # prepare mutation matrix
    mutmat = np.zeros(( 2220, 64 ))

    fname='data/precomputed/polio.ddg.txt'
    fileIN = io.open(fname, 'r')
    line = fileIN.readline()
    while line:
        current_line = line.split()
        pos = int(current_line[0])
        codon = current_line[1]

        if seq[int(pos)-1] ==  codon and len(current_line) > 4:
            mut = current_line[3]
            ddg = current_line[4]

            for ix, i in enumerate(codons):
                current_aa = gencode[i]

                if current_aa == mut:
                    mutmat[pos-1, ix] = ddg

        line = fileIN.readline()

    mm2 = mutmat > theta
    mutmat = np.array(mm2, dtype=bool)

    # compute differences
    for i in range(2200):
        OR, p = test_allpassages(i, seq, wt, ga, mutmat, False)
        diff[i,] = np.array([p, OR])

    return diff


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_stabdiff(data):

    d = numpy2ri(data)
    r.assign('data', d)

    r(' source("src/R/figure_sites_stab.R") ')
    r(' source("src/R/figure_P1_stab.R") ')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_agg(seq):
    """
    load aggregation propensity of mutations into dict 
    """

    tango_dict = {}

    fname='data/precomputed/polio.AP.txt'
    fileIN = io.open(fname, 'r')

    line = fileIN.readline()
    while line:
        current_line = line.split()
        pos = int(current_line[0])
        codon = current_line[1]
        aa1 = current_line[2]
        aa2 = current_line[3]

        if seq[int(pos)-1] == codon and len(current_line) > 4:
            id = str(pos) + str(aa1) + str(aa2)
            tango = round(float(current_line[4]), 2)

            if tango < -1000:
                tango = 0

            tango_dict[id] = tango

        line = fileIN.readline()

    return tango_dict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_stab(seq):
    """
    get dict for stability predictions
    """

    stab_dict = {}

    fname='data/precomputed/polio.ddg.txt'
    fileIN = io.open(fname, 'r')

    line = fileIN.readline()
    while line:
        current_line = line.split()
        pos = int(current_line[0])
        codon = current_line[1]
        aa1 = current_line[2]
        aa2 = current_line[3]

        if seq[int(pos)-1] ==  codon and len(current_line) > 4:

            id = str(pos) + str(aa1) + str(aa2)
            ddg = current_line[4]

            stab_dict[id] = ddg

        line = fileIN.readline()

    return stab_dict





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_hyd(wt, ga):
    """
    test for sites that differ in hydrophobicity of the  mutations
    """

    # initiate
    diff = np.zeros(( 2200, 2 ))
    seq = loadRefSeq('data/reference/polio.codonseq')
    theta = 4           # threshold of mutations

    # prepare mutation matrix
    hydrophobics = ['V', 'L', 'I', 'M', 'F', 'Y', 'W']
    mutmat = np.zeros(( 64 ))
    for ix, i in enumerate(codons):
        if gencode[i] in hydrophobics:
            mutmat[ix] = 1
    mutmat = np.tile(mutmat, (2220, 1) )

    mutmat = np.zeros(( 2220, 64 ))
    for i in range(2200):
        aa_wt = gencode[ seq[i] ]
        for j in range(64):
            aa2 = gencode[ codons[j] ]
            if aa_wt != "*" and aa2 != "*":
                mutmat[i,j] = aa_hydrophobicity[aa2] - aa_hydrophobicity[aa_wt]

    mm2 = mutmat > theta
    mutmat = np.array(mm2, dtype=bool)


    # compute differences
    for i in range(2200):
        OR, p = test_allpassages(i, seq, wt, ga, mutmat, False)
        diff[i,] = np.array([p, OR])

    return diff


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_hyddiff(data):

    d = numpy2ri(data)
    r.assign('data', d)

    r(' source("src/R/figure_sites_hyd.R") ')
    r(' source("src/R/figure_P1_hyd.R") ')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_stabhyddiff(stab, hyd):

    d1 = numpy2ri(stab)
    r.assign('stab', d1)

    d2 = numpy2ri(hyd)
    r.assign('hyd', d2)

    r(' source("src/R/figure_sites_stab_nothyd.R") ')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_sigsites(data):
    """
    print list of sig sites
    """

    # adjust p-value
    data[:,0] = pval_adjust(data[:,0])

    # only P1
    data = data[0:880,:] 

    sel_sig = data[:,0] < 0.05
    sel_wt = data[:,1] > 1
    sel_ga = data[:,1] < 1

    sel_sigwt = sel_sig * sel_wt
    sel_sigga = sel_sig * sel_ga

    pos_sigwt = np.where(sel_sigwt)[0]
    pos_sigga = np.where(sel_sigga)[0]

    print pos_sigwt
    print pos_sigga

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pval_adjust(pval):
    """ 
    correct for multiple testing
    Benjamini & Hochberg FDR method (R funct p.adjust)
    """ 
    pval = np.array(pval)
    if np.size(pval) > 1:
        padj = np.nan * np.ones( (np.size(pval) ))
        nn = np.isnan(pval)
        pval = pval[~nn]
        n  = len(pval)
        i  = np.array( range(n)[::-1]) + 1
        o  = np.array( sorted(range(n), key=lambda k: pval[k])[::-1] )
        ro = np.array( sorted(o, key=lambda k: o[k]) )
        adj = np.minimum.accumulate( float(n)/np.array(i) * pval[o] ) #, array(pval)[o]
        for i in range(len(adj)):
            adj[i] = min(1, adj[i])
        padj[~nn] = adj[ro]
    else:
        padj = pval
    
    return padj



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_allmutations(pos, cc_wt, cc_ga):
    """
    test all possible mutations of a site individually
    across all passages
    """

    refseq = loadRefSeq('data/reference/polio.codonseq')
    allstab = load_stab(refseq)
    alltango = load_agg(refseq)

    current_wt_codon = refseq[pos]
    current_wt_idx = codons.index(current_wt_codon)
    current_neighbours = neighbours[current_wt_idx,:]
    current_neighbours_idx = np.where(current_neighbours)[0]

    out = np.zeros(( sum(current_neighbours), 8 ))

    counter_neighbours = 0 
    for i in current_neighbours_idx:
        mhcounts = np.zeros(( 2, 2, 7 ))
        s = np.zeros(( 7 ))

        aa1 = gencode[current_wt_codon]
        aa2 = gencode[codons[i]]
        id = str(pos+1) + str(aa1) + str(aa2)

        for j in range(7):
            mutation = codons[i]
            mut_wt = int(cc_wt[pos, i, j])
            nonmut_wt = int( cc_wt[pos,current_wt_idx,j]  )
            mut_ga = int(cc_ga[pos, i, j])
            nonmut_ga = int(cc_ga[pos, current_wt_idx, j]  )

            mhcounts[:,:,j] = np.array( [[mut_wt, mut_ga], [nonmut_wt, nonmut_ga]] )
            s[j] = np.amin( mhcounts[:,:,j] )

        nonzero = s != 0
        mhcounts = mhcounts[:,:,nonzero]

        if sum(nonzero) > 4:
            out1   = MHtest(mhcounts)
            pval_mut = out1[0]
            OR_mut = out1[1]
        else:
            pval_mut = 1
            OR_mut = 1

        if gencode[codons[i]] != "*":
            out[counter_neighbours,:] = np.array([pos, pval_mut, OR_mut, aa_hydrophobicity[ aa2 ] - aa_hydrophobicity[ aa1 ], allstab.get(id, np.nan), alltango.get(id, np.nan), aminoacids.index(aa1), aminoacids.index(aa2)])

        counter_neighbours += 1

    return out 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def random_properties(n):
    """
    geneates distribution of random mutations
    hydrophobicity; stability; aggregation
    """
    import random

    refseq = loadRefSeq('data/reference/polio.codonseq')
    allstab = load_stab(refseq)
    alltango = load_agg(refseq)

    stab = allstab.values()
    agg = alltango.values()

    out = np.zeros((n, 3))

    for i in range(n):
        out[i,0] = aa_hydrophobicity[random.choice(aminoacids)] - aa_hydrophobicity[random.choice(aminoacids)]
        out[i,1] = random.choice(stab)
        out[i,2] = random.choice(agg)

    return out

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_mutations_sigsites(data, wt, ga):
    """
    explore individual mutations that give rise to significant sites
    """

    data[:,0] = pval_adjust( data[:,0] ) # correct pval for mult. testing (BH)
    theta = 0.05

    sig = data[:,0] < theta
    pos = data[:,1] > 1
    neg = data[:,1] < 1

    sig_wt = np.where( sig * pos )[0]
    sig_ga = np.where( sig * neg )[0]

    out_wt = np.zeros(( len(sig_wt)*9, 8 ))
    out_ga = np.zeros(( len(sig_ga)*9, 8 ))

    counter_pos = 0
    for i in sig_wt:
        out_pos = test_allmutations(i, wt, ga)
        L_pos = np.shape(out_pos)[0]
        out_wt[counter_pos:counter_pos+L_pos,:] = out_pos
        counter_pos += L_pos

    counter_neg = 0
    for i in sig_ga:
        out_neg = test_allmutations(i, wt, ga)
        L_neg = np.shape(out_neg)[0]
        out_ga[counter_neg:counter_neg+L_neg,:] = out_neg
        counter_neg += L_neg

    out_wt = out_wt[:counter_pos,:]
    out_ga = out_ga[:counter_neg,:]

    return out_wt, out_ga




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_sigpos(stab_wt, stab_ga, agg_wt, agg_ga, hyd_wt, hyd_ga, rando):
    """
    filter sites of interest and plot in R
    """

    # only P1 + no AAs
    stab_wt = stab_wt[ stab_wt[:,0] < 880,:6]
    stab_ga = stab_ga[ stab_ga[:,0] < 880,:6]
    agg_wt  = agg_wt[ agg_wt[:,0] < 880,:6]
    agg_ga  = agg_ga[ agg_ga[:,0] < 880,:6]
    hyd_wt  = hyd_wt[ hyd_wt[:,0] < 880,:6]
    hyd_ga  = hyd_ga[ hyd_ga[:,0] < 880,:6]
   

    # adjust p-value
    stab_wt[:,1] = pval_adjust(stab_wt[:,1])
    stab_ga[:,1] = pval_adjust(stab_ga[:,1]) 
    agg_wt[:,1] = pval_adjust(agg_wt[:,1])
    agg_ga[:,1] = pval_adjust(agg_ga[:,1]) 
    hyd_wt[:,1] = pval_adjust(hyd_wt[:,1])
    hyd_ga[:,1] = pval_adjust(hyd_ga[:,1]) 


    stabwt = numpy2ri(stab_wt)
    r.assign('stab_wt', stabwt)
    stabga = numpy2ri(stab_ga)
    r.assign('stab_ga', stabga)

    aggwt = numpy2ri(agg_wt)
    r.assign('agg_wt', aggwt)
    aggga = numpy2ri(agg_ga)
    r.assign('agg_ga', aggga)

    hydwt = numpy2ri(hyd_wt)
    r.assign('hyd_wt', hydwt)
    hydga = numpy2ri(hyd_ga)
    r.assign('hyd_ga', hydga)

    rand = numpy2ri(rando)
    r.assign('rand', rand)

    r(' source("src/R/figure_sites_properties.R") ')
    r(' source("src/R/figure_sites_sel_properties.R") ')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mutation_network(stab_wt, stab_ga, agg_wt, agg_ga):
    """
    parse significant mutations in network format
    """

    theta = 0.05 	# significance level    

    # only P1
    stab_wt = stab_wt[ stab_wt[:,0] < 880,:]
    stab_ga = stab_ga[ stab_ga[:,0] < 880,:]
    agg_wt  = agg_wt[ agg_wt[:,0] < 880,:]
    agg_ga  = agg_ga[ agg_ga[:,0] < 880,:]

    # adjust p-value
    stab_wt[:,1] = pval_adjust(stab_wt[:,1])
    stab_ga[:,1] = pval_adjust(stab_ga[:,1])
    agg_wt[:,1] = pval_adjust(agg_wt[:,1])
    agg_ga[:,1] = pval_adjust(agg_ga[:,1])


    # only significant mutations
    stab_wt = stab_wt[ stab_wt[:,1] < theta,:]
    stab_ga = stab_ga[ stab_ga[:,1] < theta,:]
    agg_wt = agg_wt[ agg_wt[:,1] < theta,:]
    agg_ga = agg_ga[ agg_ga[:,1] < theta,:]

    # individual mutation with same trend as significant site
    stab_wt = stab_wt[ stab_wt[:,2] > 1,:]
    stab_ga = stab_ga[ stab_ga[:,2] < 1,:]
    agg_wt = agg_wt[ agg_wt[:,2] > 1,:]
    agg_ga = agg_ga[ agg_ga[:,2] < 1,:]

    # individual mutation of same kind as significant site
    stab_wt = stab_wt[ stab_wt[:,4] > 1,6:8]
    stab_ga = stab_ga[ stab_ga[:,4] > 1,6:8]
    agg_wt = agg_wt[ agg_wt[:,5] > 2,6:8]
    agg_ga = agg_ga[ agg_ga[:,5] > 2,6:8]


    def condense_mutations(data, anno):
        """
        get dictionary of counts of each mutation
        """
        
        out = {}
        for i in range( np.shape(data)[0] ):

            aa1 = aminoacids[ int(data[i,0]) ]
            aa2 = aminoacids[ int(data[i,1]) ] 

            mutation = str(aa1) + "_" + str(aa2) + "_" + anno

            if mutation in out.keys():
                out[mutation] += 1
            else:
                out[mutation] = 1

        return out

    mut_stab_wt = condense_mutations(stab_wt, "stab_wt")
    mut_stab_ga = condense_mutations(stab_ga, "stab_ga")
    mut_agg_wt  = condense_mutations(agg_wt, "agg_wt")
    mut_agg_ga  = condense_mutations(agg_ga, "agg_ga")

    mutnet = mut_stab_wt.copy()
    mutnet.update(mut_stab_ga)
    mutnet.update(mut_agg_wt)
    mutnet.update(mut_agg_ga)

    with open('mutnetwork.tab', 'w') as f:
        [f.write('{0}\t{1}\n'.format(key, value)) for key, value in mutnet.items()]







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUPPLEMENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_average_coverage_supplement(wt, ga):

    wt = np.sum(wt, 1)
    ga = np.sum(ga, 1)

    data_wt = numpy2ri(wt)
    data_ga = numpy2ri(ga)

    r.assign('wt', data_wt)
    r.assign('ga', data_ga)

    r(' wt <- as.matrix(wt) ')
    r(' ga <- as.matrix(ga) ')
    r(' source("src/R/figure_supplement_coverage.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_entropy_supplement(ent_wt, ent_ga, div_wt, div_ga):

    ent_wt = numpy2ri(ent_wt)
    ent_ga = numpy2ri(ent_ga)
    div_wt = numpy2ri(div_wt)
    div_ga = numpy2ri(div_ga)

    r.assign('entwt', ent_wt)
    r.assign('entga', ent_ga)
    r.assign('divwt', div_wt)
    r.assign('divga', div_ga)

    r(' source("src/R/figure_supplement_entropy.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_sites_supplement(data_stab, data_agg, data_hyd):

    stab = numpy2ri(data_stab)
    r.assign('stab', stab)
    agg = numpy2ri(data_agg)
    r.assign('agg', agg)
    hyd = numpy2ri(data_hyd)
    r.assign('hyd', hyd)

    r(' source("src/R/figure_supplement_sites.R") ')
    r(' source("src/R/figure_supplement_P1_stab.R") ')
    r(' source("src/R/figure_supplement_P1_agg.R") ')
    r(' source("src/R/figure_supplement_P1_hyd.R") ')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_dnds_supplement(wt1, ga1, wt2_P1, ga2_P1, wt2, ga2):

    wt1R = numpy2ri(wt1)
    ga1R = numpy2ri(ga1)
    r.assign('wt1', wt1R)
    r.assign('ga1', ga1R)

    wt2_P1R = numpy2ri(wt2_P1)
    ga2_P1R = numpy2ri(ga2_P1)
    r.assign('wt2_P1', wt2_P1R)
    r.assign('ga2_P1', ga2_P1R)

    wt2R = numpy2ri(wt2)
    ga2R = numpy2ri(ga2)
    r.assign('wt2', wt2R)
    r.assign('ga2', ga2R)

    r(' source("src/R/figure_supplement_dnds.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_mutrate_supplement(wt, ga):

    ratemat = np.zeros(( 7, 2 ))

    for i in range(7):
        ratemat[i,0] = mutation_rate(wt[:,:,i])
        ratemat[i,1] = mutation_rate(ga[:,:,i])

    d = numpy2ri(ratemat)
    r.assign('data', d)

    r(' source("src/R/figure_supplement_mutrate.R") ')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_muttai_supplement(r1_hela, r1_human, r2_hela, r2_human):

    r1hela = numpy2ri(r1_hela)
    r.assign('r1hela', r1hela)
    r1human = numpy2ri(r1_human)
    r.assign('r1human', r1human)
    r2hela = numpy2ri(r2_hela)
    r.assign('r2hela', r2hela)
    r2human = numpy2ri(r2_human)
    r.assign('r2human', r2human)

    r(' source("src/R/figure_supplement_tAI.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_tediff_supplement(r1_hela, r1_human, r2_hela, r2_human):

    r1hela = numpy2ri(r1_hela)
    r.assign('r1hela', r1hela)
    r1human = numpy2ri(r1_human)
    r.assign('r1human', r1human)
    r2hela = numpy2ri(r2_hela)
    r.assign('r2hela', r2hela)
    r2human = numpy2ri(r2_human)
    r.assign('r2human', r2human)

    r(' source("src/R/figure_supplement_nopt.R") ')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_sigpos_supplement(stab_wt, stab_ga, agg_wt, agg_ga, hyd_wt, hyd_ga, rando):
    """
    filter sites of interest and plot in R
    """

    # only P1 + no AAs
    stab_wt = stab_wt[ stab_wt[:,0] < 880,:6]
    stab_ga = stab_ga[ stab_ga[:,0] < 880,:6]
    agg_wt  = agg_wt[ agg_wt[:,0] < 880,:6]
    agg_ga  = agg_ga[ agg_ga[:,0] < 880,:6]
    hyd_wt  = hyd_wt[ hyd_wt[:,0] < 880,:6]
    hyd_ga  = hyd_ga[ hyd_ga[:,0] < 880,:6]

    # adjust p-value
    stab_wt[:,1] = pval_adjust(stab_wt[:,1])
    stab_ga[:,1] = pval_adjust(stab_ga[:,1])
    agg_wt[:,1] = pval_adjust(agg_wt[:,1])
    agg_ga[:,1] = pval_adjust(agg_ga[:,1])
    hyd_wt[:,1] = pval_adjust(hyd_wt[:,1])
    hyd_ga[:,1] = pval_adjust(hyd_ga[:,1])

    stabwt = numpy2ri(stab_wt)
    r.assign('stab_wt', stabwt)
    stabga = numpy2ri(stab_ga)
    r.assign('stab_ga', stabga)

    aggwt = numpy2ri(agg_wt)
    r.assign('agg_wt', aggwt)
    aggga = numpy2ri(agg_ga)
    r.assign('agg_ga', aggga)

    hydwt = numpy2ri(hyd_wt)
    r.assign('hyd_wt', hydwt)
    hydga = numpy2ri(hyd_ga)
    r.assign('hyd_ga', hydga)

    rand = numpy2ri(rando)
    r.assign('rand', rand)

    r(' source("src/R/figure_supplement_sites_sel_properties.R") ')


