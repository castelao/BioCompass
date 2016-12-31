# -*- coding: utf-8 -*-

""" Biocompass criteria to compare two genes
"""

import re


def score_match(gene1, gene2, criteria=None):
    score = 0
    score += criteria_category(gene1, gene2)
    score += criteria_best_git_BGC(gene1, gene2)

    return score


def criteria_category(gene1, gene2):
    """ Compare the category of the genes
    """
    if (gene1.category == 'hypothetical') and \
            (gene2.category == 'hypothetical'):
                return 1
    elif (gene1.category == 'hypothetical') or \
            (gene2.category == 'hypothetical'):
                return 2
    elif (gene1.category == gene2.category):
            return 5


def criteria_best_git_BGC(gene1, gene2):
    """ Compare ???
    """
    score = 0
    if gene1.best_hit_BGC != 'None' and gene2.best_hit_BGC != 'None':
        if gene1.best_hit_BGC == gene2.best_hit_BGC:
            score += 2
            gene1_best_hit_pos = re.search(
                    r'^\D*([0-9]*)', gene1.best_hit_gene_loc)
            gene2_best_hit_pos = re.search(
                    r'^\D*([0-9]*)', gene2.best_hit_gene_loc)
            dif_best_hit_pos = abs(
                    abs(int(gene2_best_hit_pos.group(1)) \
                           - int(gene1_best_hit_pos.group(1))) \
                    - abs(int(gene2.n) - int(gene1.n)))
            if dif_best_hit_pos == 0:
                score += 3
            elif dif_best_hit_pos == 1:
                score += 2
            elif dif_best_hit_pos == 2:
                score += 1
    else:
        score += 1

    return score
