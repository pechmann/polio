######################################################################################################################################################################
#	Copyright 2017, Sebastian Pechmann, Universite de Montreal (The MIT License)
#
#	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
#	to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#	and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
#	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################################################################################################################

import sys
sys.path.insert(0, 'src')
from utils import *


include_figure2 = True
include_figure3 = True
include_figure4 = True
include_supplement = True


def figure2(wt, ga):
    """
    reproduce panels from Figure 2
    """

    #cirseq coverage
    plot_average_coverage(wt, ga)

    #virus mutation rate
    plot_mutrate(wt, ga)

    #evolutionary rate of P1 
    dnds_wt = compute_dnds(wt[0:880,:,:])
    dnds_ga = compute_dnds(ga[0:880,:,:])
    plot_dnds(dnds_wt, dnds_ga)

    # mutational diversity
    entropy_wt = compute_entropy(wt)
    entropy_ga = compute_entropy(ga)
    plot_entropy(entropy_wt, entropy_ga)

    diversity_wt = compute_diversity(wt)
    diversity_ga = compute_diversity(ga)
    plot_diversity(diversity_wt, diversity_ga)


def figure3(wt, ga):
    """
    reproduce panels from Figure 3
    """

    aggdiff = test_agg(wt, ga)
    plot_aggdiff(aggdiff)

    stabdiff = test_stab(wt, ga)
    plot_stabdiff(stabdiff)

    hyddiff = test_hyd(wt, ga)
    plot_hyddiff(hyddiff)

    stab_wt, stab_ga = test_mutations_sigsites(stabdiff, wt, ga)
    agg_wt, agg_ga = test_mutations_sigsites(aggdiff, wt, ga)
    hyd_wt, hyd_ga = test_mutations_sigsites(hyddiff, wt, ga)
    rand = random_properties(1000)
    plot_sigpos(stab_wt, stab_ga, agg_wt, agg_ga, hyd_wt, hyd_ga, rand)

    #mutation_network(stab_wt, stab_ga, agg_wt, agg_ga)
    
    plot_stabhyddiff(stabdiff, hyddiff)


def figure4(wt, ga):
    """
    reproduce panels from Figure4
    """

    mutation_tAI = compute_mutationtAI(wt, ga, "hela")
    plot_muttai(mutation_tAI)

    tediff = test_te(wt, ga, "hela")
    plot_tediff(tediff)


def supplement(wt1, ga1, wt2, ga2):
    """
    reproduce panels from Supplement
    """

    # coverage
    plot_average_coverage_supplement(wt2, ga2)

    # mutation rate
    plot_mutrate_supplement(wt2, ga2)


    #dnds replica 1, full virus
    dnds_wt1 = compute_dnds(wt1)
    dnds_ga1 = compute_dnds(ga1)

    #dnds replica 2, P1
    dnds_wt2_P1 = compute_dnds(wt2[0:880,:,:])
    dnds_ga2_P1 = compute_dnds(ga2[0:880,:,:])

    #dnds replica2, full virus
    dnds_wt2 = compute_dnds(wt2)
    dnds_ga2 = compute_dnds(ga2)
    plot_dnds_supplement(dnds_wt1, dnds_ga1, dnds_wt2_P1, dnds_ga2_P1, dnds_wt2, dnds_ga2)  

    # mutational diversity
    entropy_wt2 = compute_entropy(wt2)
    entropy_ga2 = compute_entropy(ga2)
    diversity_wt2 = compute_diversity(wt2)
    diversity_ga2 = compute_diversity(ga2)
    plot_entropy_supplement(entropy_wt2, entropy_ga2, diversity_wt2, diversity_ga2)

    aggdiff2 = test_agg(wt2, ga2)
    stabdiff2 = test_stab(wt2, ga2)
    hyddiff2 = test_hyd(wt2, ga2)
    plot_sites_supplement(stabdiff2, aggdiff2, hyddiff2)


    stab_wt2, stab_ga2 = test_mutations_sigsites(stabdiff2, wt2, ga2)
    agg_wt2, agg_ga2 = test_mutations_sigsites(aggdiff2, wt2, ga2)
    hyd_wt2, hyd_ga2 = test_mutations_sigsites(hyddiff2, wt2, ga2)
    rand = random_properties(1000)
    plot_sigpos_supplement(stab_wt2, stab_ga2, agg_wt2, agg_ga2, hyd_wt2, hyd_ga2, rand)

    #mutation_network(stab_wt2, stab_ga2, agg_wt2, agg_ga2)

    mutation_tAI_hela_r1 = compute_mutationtAI(wt1, ga1, "hela")
    mutation_tAI_human_r1 = compute_mutationtAI(wt1, ga1, "human")
    mutation_tAI_hela_r2 = compute_mutationtAI(wt2, ga2, "hela")
    mutation_tAI_human_r2 = compute_mutationtAI(wt2, ga2, "human")
    plot_muttai_supplement(mutation_tAI_hela_r1, mutation_tAI_human_r1, mutation_tAI_hela_r2, mutation_tAI_human_r2)   

    tediff_r1_hela = test_te(wt1, ga1, "hela")
    tediff_r1_human = test_te(wt1, ga1, "human")
    tediff_r2_hela = test_te(wt2, ga2, "hela")
    tediff_r2_human = test_te(wt2, ga2, "human")
    plot_tediff_supplement(tediff_r1_hela, tediff_r1_human, tediff_r2_hela, tediff_r2_human)



if __name__ == '__main__':

    # load data
    rep1_wt, rep1_ga = loadBatch('replica1')
    rep2_wt, rep2_ga = loadBatch('replica2')

    if include_figure2:
        figure2(rep1_wt, rep1_ga)

    if include_figure3:
        figure3(rep1_wt, rep1_ga)

    if include_figure4:
        figure4(rep1_wt, rep1_ga)

    if include_supplement:
        supplement(rep1_wt, rep1_ga, rep2_wt, rep2_ga)



