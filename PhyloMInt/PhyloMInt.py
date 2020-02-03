#!/usr/bin/python3
import sys
import os
import subprocess
from pathlib import Path
from libsbml import *
from argparse import ArgumentParser
from lib import BuildGraphNetX, CalculateIndexes

def singleFile(infile1, infile2):
    '''
    Processes competition and cooperation index if two sbml files are provided
    Input: SBML A and B
    Returns: tabulated competition and cooperation index
    '''

    A = os.path.basename(infile1).rstrip('.xml')
    B = os.path.basename(infile2).rstrip('.xml')

    DG_A=BuildGraphNetX.buildDG(args.infile1)
    DG_B=BuildGraphNetX.buildDG(args.infile2)

    # Calculate the seed sets
    SeedSetAConfidence, SeedSetA, nonSeedSetA = BuildGraphNetX.getSeedSet(DG_A, maxComponentSize=args.maxcc)
    SeedSetBConfidence, SeedSetB, nonSeedSetB = BuildGraphNetX.getSeedSet(DG_B, maxComponentSize=args.maxcc)

    # Get Comeptition
    MetabolicCompetitionIdxAB = CalculateIndexes.MetabolicCompetitionIdx(SeedSetAConfidence, SeedSetBConfidence)
    MetabolicCooperationIdxAB = CalculateIndexes.MetabolicCooperationIdx(SeedSetAConfidence, SeedSetBConfidence, nonSeedSetB)

    return (f'{A}\t{B}\tCompetition:{MetabolicCompetitionIdxAB}\tCooperation:{MetabolicCooperationIdxAB}')

def runCarveMe(fna_path, outdir, FGS):
    ''' 
    Takes input directory of genomes, predicts CDS, runs CarveMe
    Input: directory path to genomes
    Output: directory path to SBML files
    '''
    
    # make temporary file for FGS 
    Path(f'{outdir}/tmp/FragGeneScan/').mkdir(parents=True, exist_ok=True)
    
    # acceptable suffixes
    my_suffixes = ('fna','fa','fq','fastq')

    # run FGS
    fna_files = [f for f in os.listdir(fna_path) if f.endswith(my_suffixes)]

    print('FragGeneScan: Predicting CDS from genome sequences')
    for fna in fna_files:
        #run FGS
        fna_prefix = os.path.splitext(os.path.basename(fna))[0]
        subprocess.run([FGS, '-s', f'{fna_path}/{fna}', '-o', f'{outdir}/tmp/FragGeneScan/{fna_prefix}', '-w', '1', '-t', 'complete',])

    # make dir for CarveMe output
    Path(f'{outdir}/tmp/CarveMe/').mkdir(parents=True, exist_ok=True)

    # run CarveMe
    print('CarveMe: Predicting GENRE')
    for faa in fna_files:
        #carve ${file}/${f}.faa -o ${OUTDIR}/${f}.xml
        faa_prefix = os.path.splitext(os.path.basename(faa))[0]
        print('faa_prefix')
        subprocess.run(['carve', f'{outdir}/tmp/FragGeneScan/{faa_prefix}.faa', '-o', f'{outdir}/tmp/CarveMe/{faa_prefix}.xml'])

def directoryALL(dirr_path, outdir, outfile):
    '''
    Processess all SBML(XML) files in the directory path.
    Input: dirrectory path, outpath
    Output: file with competition and cooperation
    '''

    # build initial dictionaries
    temp_set = set()
    SeedSetDic = dict()
    nonSeedSetDic = dict()
    ConfidenceDic = dict()

    # get all XML files in directory        
    sbml_files = [ f for f in os.listdir(dirr_path) if f.endswith('.xml') ]

    count = 0
    total = len(sbml_files)

    # fill dictionary
    for sbml in sbml_files:

        count += 1

        sbml_path = ('%s/%s' % (dirr_path, sbml))
        sbml_base = sbml.rstrip('.xml')

        print(f'Calculating Seeds: {sbml_base} - {count}/{total}')

        if sbml_base not in temp_set:

            # calculate SeedSets
            DG_sbml = BuildGraphNetX.buildDG(sbml_path)
            SeedSetConfidence, SeedSet, nonSeedSet = BuildGraphNetX.getSeedSet(DG_sbml, maxComponentSize=args.maxcc)

            # Save dicts
            SeedSetDic[sbml_base] = SeedSet
            nonSeedSetDic[sbml_base] = nonSeedSet
            ConfidenceDic[sbml_base] = SeedSetConfidence
            temp_set.add(sbml_base)

    # calculate pairwise competition & cooperation index
    count = 0
    pairwise_total = total*total

    with open(f'{outdir}/{outfile}', 'w') as out_file:
        out_file.write(f'\tA\tB\tCompetition\tComplementarity\n')
        for A in sbml_files:
            for B in sbml_files:
                
                A = A.rstrip('.xml')
                B = B.rstrip('.xml')
                
                count += 1
                print(f'Calculating Indexes: {A} vs {B} - {count}/{pairwise_total}')
    
                MetabolicCompetitionIdxAB = CalculateIndexes.MetabolicCompetitionIdx(ConfidenceDic[A], ConfidenceDic[B])
                MetabolicCooperationIdxAB = CalculateIndexes.MetabolicCooperationIdx(ConfidenceDic[A], ConfidenceDic[B], nonSeedSetDic[B])

                output = (f'{A}\t{B}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n')
                out_file.write(output)

if __name__ == '__main__' :
    parser = ArgumentParser(description = 'PhyloMInt: Takes input genomes or genome scale metabolic models in SBML format. Calculates metabolic complementarity and cooperation indexes. Outputs to output file, tsv format.')
    parser.add_argument('-i', '--infile1', help = 'Input SBML file (XML format) organism A.')
    parser.add_argument('-j', '--infile2', help = 'Input SBML file (XML format) organism B.')
    parser.add_argument('-c', '--carveme', help = '(Optional) Path to fasta files. Will run FragGeneScan, predict genes to be used as CarveMe output.', default=None)
    parser.add_argument('-d', '--dirr', help = '(Optional) Path to directory with SBML files. Will compute all pairwise comparisons & use dynamic programming for computational speed-up.', default=None)
    parser.add_argument('-o', '--outfile', help = 'Output file name.', required=True)
    parser.add_argument('--outdir', help = 'Outdir path. default = cwd', default = os.getcwd()) 
    parser.add_argument('--maxcc', help = 'Maximum number of nodes in a strongly connected component (SCC) to consider in SeedSet', default=5)
    parser.add_argument('--fraggenescan', help = 'path to FragGeneScan. (default: lib/FragGeneScan1.30/FragGeneScan)', default = f'{os.path.dirname(os.path.abspath(__file__))}/lib/FragGeneScan1.30/FragGeneScan')
    args = parser.parse_args()


    # if dirr not used
    if args.dirr == None and args.carveme == None:
        
        # get indexes
        output = singleFile(args.infile1, args.infile2)

        # write to outfile
        with open(f'{args.outdir}/{args.outfile}', 'a+') as outfile:
            outfile.write(output)

    # raise exception, only allow start from Carveme or SBML.
    elif args.dirr != None and args.carveme != None:
        raise ValueError('Error: Both CarveMe and SBML files have been selected. Please choose to execute starting from genome files (CarveMe), or SBML files.')

    # run pipeline starting with SBML files in path
    elif args.dirr != None and args.carveme == None:
        # process all files in directory path
        directoryALL(args.dirr, args.outdir, args.outfile)

    # start pipeline with genome files
    elif args.dirr == None and args.carveme != None:
        runCarveMe(args.carveme, args.outdir, args.fraggenescan)
        directoryALL(f'{args.outdir}/tmp/CarveMe', args.outdir, args.outfile)


