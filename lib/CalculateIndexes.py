#!/usr/bin/env python3

# calculate competition as a normalized weighted sum
def MetabolicCompetitionIdx(SeedSetAConfidence, SeedSetBConfidence):
    ''' 
    Calculate metabolic compeition index between SetA and Set B
    Input: SeedSet dictionary with confidence scores
    Return: Normalized weight sum metabolic competition index
    '''

    # get seed set A & B
    SeedA = set(SeedSetAConfidence.keys())
    SeedB = set(SeedSetBConfidence.keys())

    # get intersects
    intersectAB = SeedA.intersection(SeedB)

    # store weighted sum
    normIntersect = 0.0
    
    # get weighted total
    sumA = sum(SeedSetAConfidence.values())

    # count intersect of confidence scores
    for seed in intersectAB:
        normIntersect += SeedSetAConfidence[seed]

    # calculate normalized weighted sum
    MetabolicCompetitionIdx = (normIntersect/sumA)

    return MetabolicCompetitionIdx

# calculated cooperation with seedsA vs nonseedsB
def MetabolicCooperationIdx(SeedSetAConfidence, SeedSetBConfidence, nonSeedB):
    '''
    Calculate metabolic cooperation index between SeedSetA and nonSeedSetB
    Normalized weighted sum
    '''
    
    SeedA = set(SeedSetAConfidence.keys())
    SeedB = set(SeedSetBConfidence.keys())
    nonSeedB = set(nonSeedB)
    SetB = (SeedB|nonSeedB)

    # get intersects intersect (A n nonB) &! B
    intersect_seedA_nonseedB = SeedA.intersection(nonSeedB)
    intersect_seedA_setB = SeedA.intersection(SetB)

    # calculate normalized weighted sum
    MetabolicCooperationIdx = (len(intersect_seedA_nonseedB)/len(intersect_seedA_setB))

    return MetabolicCooperationIdx

