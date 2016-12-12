#!/usr/bin/env python
'''
Summary: This script will calculate CSS of pair of reverse repeat sequence.

Version: 0.1.0

Author: Rui Dong; Xiao-Ou Zhang; Xu-Kai Ma

Usage: cs_xm.py [options] -g GENOME circ_file

Options:
    -h --help                       Show help message.
    -v --version                    Show version.
    -g GENOME --genome=GENOME       Genome fasta file.
    -l LENGTH --length=LENGTH       Minimum element length. [default: 50]
    -p THREAD --thread=THREAD       Running threads. [default: 10]
    -o OUT --output=OUT             Output file. [default: circ_cs]
    --tmp                           Keep temporary BLAST results.
'''

import time
import re
import pysam
import subprocess
import os
import argparse
from multiprocessing import Pool
from collections import defaultdict
from file_parse import check_fasta
from dir_func import create_temp, delete_temp

def fetchFa(fa, fPath, intron):
    chrom, start, end = re.split(':|-', intron)
    start = int(start)
    end = int(end)
    fasta = fa.fetch(chrom, start, end)
    with open(fPath, 'w') as f:
        f.write('>' + intron + '\n')
        f.write(fasta + '\n')
    return (start, end)

def runBlast(file1, file2):
    #with open(file1) as f:
        #for line in f:
            #print(line)
    p = subprocess.Popen('blastn -query %s -subject %s ' % (file1, file2) +
                         '-word_size 11 -gapopen 5 -gapextend 2 -penalty -3 ' +
                         '-reward 2 -strand minus -outfmt 6',
                         shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, close_fds=True)
    out, err = p.communicate()
    return out.decode('utf-8')

def calPairScore(left, right, score):
    distance = ((right - left) * 1.0 / 1000) ** 2
    return score / distance

def filterBlast(blastResult, offset1, offset2, leftEnd, length, WITHIN_FLAG=False):
    blastInfo = []
    for line in blastResult.split('\n')[:-1]:
        loc1, loc2, loc4, loc3 = [ int(x) for x in line.split()[6:10] ]
        e, blastScore = [ float(x) for x in line.split()[10:] ]
        if WITHIN_FLAG and loc3 < loc1: continue
        if loc2 - loc1 < length or loc4 - loc3 < length: continue
        if e > 1e-5: continue

        region1Start = offset1 + loc1
        region1End = offset1 + loc2
        region2Start = offset2 + loc3
        region2End = offset2 + loc4

        #pairScore = calPairScore(region1End, region2Start, blastScore)
        if WITHIN_FLAG:
            pairScore = calPairScore(region1End, region2Start, blastScore)
        else:
            #pairScore = calPairScore(region1End, region2Start,
                                     #blastScore)
            pairScore = calPairScore(region1End, region2Start-offset2+leftEnd,
                                     blastScore)
            #symmetryScore = calSymmetryScore(region1End, leftEnd, offset2,
                                            #region2Start)
            #pairScore = pairScore * symmetryScore # for get effect pair ability
        blastInfo.append([region1Start, region1End,
                          region2Start, region2End,
                          pairScore])

    return blastInfo

def reflect(blastInfo, length):
    reflection = {}
    score = defaultdict(float)
    if not blastInfo:
        return (score, reflection)
    boundary = []
    scoreInfo = defaultdict(float)
    for info in blastInfo:
        (region1Start, region1End,
         region2Start, region2End,
         pairScore) = info
        
        region1 = '%d\t%d' % (region1Start, region1End)
        region2 = '%d\t%d' % (region2Start, region2End)

        if region1 not in scoreInfo:
            boundary.append([region1Start, region1End])
        scoreInfo[region1] += pairScore
        if region2 not in scoreInfo:
            boundary.append([region2Start, region2End])
        scoreInfo[region2] += pairScore

    boundary.sort()

    b1 = boundary[0]
    b1Info = '\t'.join(str(x) for x in b1)
    reflection[b1Info] = b1Info
    score[b1Info] = scoreInfo[b1Info]

    for b2 in boundary[1:]:
        b2Info = '\t'.join(str(x) for x in b2)
        if b1[1] - b2[0] >= length: # has overlap
            reflection[b2Info] = b1Info
            score[b1Info] += scoreInfo[b2Info]
        else: # no overlap
            reflection[b2Info] = b2Info
            score[b2Info] = scoreInfo[b2Info]
            # refresh element may can put outside it!
            b1 = b2
            b1Info = b2Info

    return (score, reflection)

def calSymmetryScore(left1, left2, right1, right2):
    leftLength = left2 - left1
    rightLength = right2 - right1
    if leftLength <= rightLength:
        return leftLength * 1.0 / rightLength
    else:
        return rightLength * 1.0 / leftLength

def calCompeteScore(start, end, scoreInfo, reflection, blastInfo):
    competeScore = 0
    for info in blastInfo:
        (region1Start, region1End,
         region2Start, region2End,
         pairScore) = info

        overlap1 = overlap(start, end, region1Start, region1End)
        overlap2 = overlap(start, end, region2Start, region2End)
        withinScore = pairScore

        #if overlap1 >= (end - start) * 1.0 / 2:
            ## I think here may use region2
            #region1 = '%d\t%d' % (region1Start, region1End)
            #competePotential = withinScore / scoreInfo[reflection[region1]]
            #competeScore += withinScore * competePotential
        #if overlap2 >= (end - start) * 1.0 / 2:
            ## I think here may use region1
            #region2 = '%d\t%d' % (region2Start, region2End)
            #competePotential = withinScore / scoreInfo[reflection[region2]]
            #competeScore += withinScore * competePotential

        if overlap1 >= (end - start) * 1.0 / 2:
            # I think here may use region2
            region2 = '%d\t%d' % (region2Start, region2End)
            competePotential = withinScore / scoreInfo[reflection[region2]]
            competeScore += withinScore * competePotential
        if overlap2 >= (end - start) * 1.0 / 2:
            # I think here may use region1
            region1 = '%d\t%d' % (region1Start, region1End)
            competePotential = withinScore / scoreInfo[reflection[region1]]
            competeScore += withinScore * competePotential

    return competeScore

def overlap(start1, end1, start2, end2):
    region = [[start1, end1],[start2, end2] ]
    region.sort()
    if region[1][0] < region[0][1]:
        return region[0][1] - region[1][0]
    else:
        return 0


def calCS(faF, circId, leftIntron, rightIntron, length, TMP_FLAG, tmpDir):
    '''
    1. Fetch intron sequence and run BLAST
    2. Filter BLASt results and create reflection between elements
    3. Calculate complementary score
    '''

    fa = pysam.FastaFile(faF)
    chrom, start, end = circId.split()
    start = int(start)
    end = int(end)

    # creat temporary files
    if TMP_FLAG:
        folder = '%s_%d_%d' % (chrom, start, end)
        os.mkdir('%s/%s' % (tmpDir, folder))
        leftF = '%s/%s/left_intron.fa' % (tmpDir, folder)
        rightF = '%s/%s/right_intron.fa' % (tmpDir, folder)
        acrossF = '%s/%s/across_blast.txt' % (tmpDir, folder)
    else:
        tempDir, leftF, rightF = create_temp()

    # fetch intron sequences
    leftStart, leftEnd = fetchFa(fa, leftF, leftIntron)
    rightStart, rightEnd = fetchFa(fa, rightF, rightIntron)

    # blast alignments
    acrossBlast = runBlast(leftF, rightF)
    leftBlast = runBlast(leftF, leftF)
    rightBlast = runBlast(rightF, rightF)

    #print(leftBlast)

    if TMP_FLAG:
        with open(acrossF, 'w') as f:
            f.write(acrossBlast)

    # filter blast results
    acrossBlastInfo = filterBlast(acrossBlast, leftStart, rightStart,
                                  leftEnd, length)
    leftBlastInfo = filterBlast(leftBlast, leftStart, leftStart,
                                leftEnd, length, WITHIN_FLAG=True)
    rightBlastInfo = filterBlast(rightBlast, rightStart, rightStart,
                                 leftEnd, length, WITHIN_FLAG=True)

    #print(leftBlastInfo)
    # create reflection between elements for within pairings
    leftScore, leftReflection = reflect(leftBlastInfo, length)
    rightScore, rightReflection = reflect(rightBlastInfo, length)

    # calculate complementary score
    maxScore, score1, score2, score3 = 0.0, 0.0, 0.0, 0.0
    leftRegion, rightRegion = 'NULL', 'NULL'
    #print(acrossBlast)
    #print(acrossBlastInfo)
    #print(leftScore)
    for info in acrossBlastInfo:
        (region1Start, region1End,
         region2Start, region2End,
         pairScore) = info

        symmetryScore = calSymmetryScore(region1End, start, end,
                                         region2Start)
        leftCompeteScore = calCompeteScore(region1Start, region1End,
                                           leftScore, leftReflection,
                                           leftBlastInfo)
        rightCompeteScore = calCompeteScore(region2Start, region2End,
                                            rightScore, rightReflection,
                                            rightBlastInfo)

        acrossScore = pairScore
        pairingPotential = acrossScore / (acrossScore + leftCompeteScore +
                                          rightCompeteScore)
        complementaryScore = symmetryScore * pairingPotential * acrossScore

        #print(symmetryScore)

        if complementaryScore > maxScore:
            maxScore = complementaryScore
            score1 = symmetryScore
            score2 = pairingPotential
            score3 = acrossScore

            leftRegion = '%s:%d-%d' % (chrom, region1Start, region1End)
            rightRegion = '%s:%d-%d' % (chrom, region2Start, region2End)

    if not TMP_FLAG:
        delete_temp(tempDir)

    return '\t'.join([circId, str(maxScore), str(score1), str(score2),
                      str(score3), leftRegion, rightRegion])


def cs(options):
    localTime = time.strftime('%H:%M:%S', time.localtime(time.time()))
    print('Start cs_xm.py at %s' % localTime)
    # check fasta file
    faF = check_fasta(options.genome, False)

    # check whether keep tmp files
    if options.tmp:
        TMP_FLAG = True
    else:
        TMP_FLAG = False

    os.mkdir(options.output)
    tmpDir = options.output
    with open(options.circ_file, 'r') as circ:
        p = Pool(options.thread)
        results = []
        for line in circ:
            circType = line.split()[13]
            if circType == 'ciRNA': continue

            leftIntron, rightIntron = line.split()[16].split('|')
            # not first/last exons
            if leftIntron == 'None' or rightIntron == 'None':
                continue

            circId = '\t'.join(line.split()[:3])
            length = options.length

            results.append(p.apply_async(calCS, args=(faF, circId,
                                                      leftIntron,
                                                      rightIntron,
                                                      length, TMP_FLAG,
                                                      tmpDir)))
        p.close()
        p.join()
        with open(options.output+'.txt', 'w') as outF:
            for r in results:
                outF.write(r.get() + '\n')

            #results.append(calCS(faF, circId,
                                 #leftIntron,
                                 #rightIntron,
                                 #length,
                                 #TMP_FLAG,
                                 #tmpDir))

        #with open(options.output+'.txt', 'w') as outF:
            #for r in results:
                #outF.write(r + '\n')

    if not TMP_FLAG:
        delete_temp(tmpDir)

        localTime = time.strftime('%H:%M:%S', time.localtime(time.time()))
        print('End cs_xm.py at %s' % localTime)

def main():
    #local_time = time.strftime('%H:%M:%S', time.localtime(time.time()))
    #print('Start CIRCexplorer2 cs at %s' % local_time)

    parser = argparse.ArgumentParser()

    parser.add_argument('circ_file',
                        help='The circular RNA file (CIRCexplorer format)')

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 0.1.0')
    parser.add_argument('-g', '--genome', dest='genome',
                        help='Genome fasta file.', required=True)
    parser.add_argument('-l', '--length', dest='length', type=int, default=50,
                        help='Minimum element length. [default: 50]')
    parser.add_argument('-p', '--thread', dest='thread', type=int, default=10,
                        help='Running threads. [default: 10]')
    parser.add_argument('-o', '--ouput', dest='output', default='circ_cs',
                        help='Output file. [default: circ_cs]')
    parser.add_argument('--tmp', default=False, action='store_true',
                        help='Keep temporary BLAST results.')

    options = parser.parse_args()
    print('###Parameters:')
    for key,val in options.__dict__.items():
        print('%s:\t%s'%(key,val))
    print('###Parameters:\n')

    cs(options)

if __name__ == '__main__':
    main()
