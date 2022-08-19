import os
import random,time
import string
import sys
import re
import argparse
from Bio import SeqIO
import subprocess
import shutil
from collections import OrderedDict
from threading import Thread
from Bio.Seq import Seq
import math
import time
#from ISSearch import ISfinder
from PdifFinder.angularPlasmid import angularPlasmid
#from echarts import echartsSequence
from PdifFinder.pdifmodulecharts import pdifmoduleechart

# 命令行设置



def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inFile", "-i", help="input nucleotide file")
    parser.add_argument("--genbankFile", "-g", help="input genbank file")
    parser.add_argument("--indir", "-d", help="input file dir")
    parser.add_argument("--outdir", "-o", help="output directory")
    parser.add_argument("--circle", "-c", default="true", help="output graph format,default is circle")
    args = parser.parse_args()
    if args.outdir:
        if args.inFile or args.genbankFile or args.indir:
            inFile = args.inFile
            outdir = args.outdir.rstrip('/')
            genbankFile = args.genbankFile
            indir = args.indir
            circle = args.circle
            return inFile, genbankFile, indir, outdir, circle
        else:
            print("Please check your input file")
            sys.exit(0)
    else:
        print("""
        *****     *  *  *****  *****  *     *  *****  ******
        *   *     *  *  *      *      *     *  *      *    *
        *****     *  *  *****  *****  *     *  *      ******
        *      ****  *  *      *      *  ****  *****  * *
        *      *  *  *  *      *      *  *  *  *      *  *
        *      ****  *  *      *      *  ****  *****  *    *
    Author:     Shao Mengjie
    Email:      1437819081@qq.com
    This is a tool for pdif site and pdif-ARGs module search.This script accepts FASTA file or GENBANK file(plasmid sequence less than 1M).
    Output includes:
           AMRgene.txt            resistance gene
           pdif_site.txt          pdif site
           pdifmodule_list.txt    pdif-ARGs module
           pdifmoduleseq.fasta      pdif-ARGs sequence
           pdifmodule.svg         pdif-ARGs figure
           plasmid.html           circular graph for above features
        """)
        sys.exit(0)


def checkInfileFormat(inFile):
    baseCharacterList = ['A', 'T', 'C', 'G', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N', 'a', 't', 'c', 'g',
                         'r', 'y', 'm', 'k', 's', 'w', 'h', 'b', 'v', 'd', 'n']
    seq = ""
    breakFlag = 0
    format = True
    f = open(inFile)
    for lines in f:
        lines = lines.strip()
        if lines.startswith('>') and seq == "":
            id = lines
        elif not lines.startswith('>'):
            seq = seq + lines
        elif lines.startswith('>') and seq != "":
            for base in seq:
                if base in baseCharacterList:
                    continue
                else:
                    format = False
                    breakFlag = 1
                    break
            if breakFlag:
                break
            id = lines
            seq = ""
    f.close()
    return format


def check_dependencies():
    blastnPath = shutil.which("blastn")
    if not blastnPath:
        print("blastn not found")
    return blastnPath



def getBMBC(pattern):
    # 预生成坏字符表
    BMBC = dict()
    for i in range(len(pattern) - 1):
        char = pattern[i]
        # 记录坏字符最右位置（不包括模式串最右侧字符）
        BMBC[char] = i + 1
    return BMBC


def getBMGS(pattern):
    # 预生成好后缀表
    BMGS = dict()

    # 无后缀仅根据坏字移位符规则
    BMGS[''] = 0

    for i in range(len(pattern)):

        # 好后缀
        GS = pattern[len(pattern) - i - 1:]

        for j in range(len(pattern) - i - 1):

            # 匹配部分
            NGS = pattern[j:j + i + 1]

            # 记录模式串中好后缀最靠右位置（除结尾处）
            if GS == NGS:
                BMGS[GS] = len(pattern) - j - i - 1
    return BMGS


def BM(string, pattern):
    """
    Boyer-Moore算法实现字符串查找
    """
    m = len(pattern)
    n = len(string)
    i = 0
    j = m
    indies = {}
    BMBC = getBMBC(pattern=pattern)  # 坏字符表
    BMGS = getBMGS(pattern=pattern)  # 好后缀表
    while i < n:
        while (j > 0):
            if i + j - 1 >= n:  # 当无法继续向下搜索就返回值
                return indies

            # 主串判断匹配部分
            a = string[i + j - 1:i + m]

            # 模式串判断匹配部分
            b = pattern[j - 1:]

            # 当前位匹配成功则继续匹配
            if a == b:
                j = j - 1

            # 当前位匹配失败根据规则移位
            else:
                i = i + max(BMGS.setdefault(b[1:], m), j - BMBC.setdefault(string[i + j - 1], 0))
                j = m

            # 匹配成功返回匹配位置
            if j == 0:
                # indies.append(i)
                indies[i] = "good"
                i += 1
                j = len(pattern)


def findFeatureEndSeq(inFile, outdir):
    seqList = []
    pos0List = []
    featureEnd = ['TAA', 'TGT', 'CAT']
    for rec in SeqIO.parse(inFile, 'fasta'):
        id = rec.id
        seq = str(rec.seq).upper()
        break
    seqTempList = {}
    for featureSeq in featureEnd:
        tempList = BM(seq, featureSeq)
        for temp in tempList.keys():
            seqTemp = seq[temp - 8:temp + 20]
            if seqTemp != "" and seqTemp not in seqTempList:
                seqTempList[seqTemp] = "good"
                pos0List.append(temp - 8)
                seqList.append(seqTemp)
    return seq, pos0List, seqList


def findFeatureEndSeqThread():
    pass


def findPdif(inFile, outdir, blastnPath, pdifDB, resistanceGenePosList):
    possiblePdifSiteOutfile = outdir + '/tmp/possible_pdif_site.txt'
    possiblePdifPairOutfileTemp = outdir + '/tmp/possible_pdif_pair.temp.txt'
    possiblePdifPairOutfileTemp2 = outdir + '/tmp/possible_pdif_pair.temp2.txt'
    possiblePdifPairOutfile = outdir + '/pdif_pair.txt'
    pdifSiteOufile = outdir + '/pdif_site1.txt'
    with open(pdifSiteOufile, 'w') as w:
        w.write('')

    seedList = []
    for rec in SeqIO.parse(pdifDB, 'fasta'):
        seq1 = str(rec.seq)
        XerC = seq1[0:11]
        XerD = seq1[17:28]
        seedList.append(XerC)
        seedList.append(XerD)
        seedList.append(XerD)
        seedList.append(XerC)
    posList = []
    possiblePdifSiteDict = OrderedDict()
    possiblePdifSiteSeqList = []
    seq, pos0List, seqList = findFeatureEndSeq(inFile, outdir)
    print("%s featured sequences" % len(seqList))
    if len(seqList) != 0:
        for num4, seq2 in enumerate(seqList):
            threadList = []
            pos0 = pos0List[num4]
            if not os.path.exists(outdir + '/tmp/seedSearch/' + str(num4)):
                os.mkdir(outdir + '/tmp/seedSearch/' + str(num4))
            batch = 1
            seedListPairLength = int(len(seedList) / 2)
            if seedListPairLength < batch:
                batch = seedListPairLength
            batchNumber = math.floor(seedListPairLength / batch)
            restNumber = seedListPairLength - batchNumber * batch
            for i in range(1, batch + 1, 1):
                start = (i - 1) * batchNumber + 1
                end = start + batchNumber - 1
                p = Thread(target=findMatchFragmentThread, args=(num4, pos0, outdir, i, seedList, seq2, start, end,))
                threadList.append(p)
            if restNumber != 0:
                p1 = Thread(target=findMatchFragmentThread, args=(
                    num4, pos0, outdir, batch + 1, seedList, seq2, seedListPairLength - restNumber + 1,
                    seedListPairLength,))
                threadList.append(p1)

            for thread in threadList:
                thread.start()
            for thread in threadList:
                thread.join()
        for file in os.listdir(outdir + '/tmp/seedSearch/'):
            indir = outdir + '/tmp/seedSearch/' + file
            for file1 in os.listdir(indir):
                infile = indir + '/' + file1
                f = open(infile)
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if line:
                        posList.append(int(line))
                f.close()
        if posList != []:
            with open(possiblePdifSiteOutfile, 'w') as w:
                w.write('')
            i = 0
            newPosList = sorted(list(set(posList)))
            newPosList1 = []
            for pos in newPosList:
                seqC = seq[pos:pos + 11]
                seqD = seq[pos + 17:pos + 28]
                seqSpacer = seq[pos + 11:pos + 17]
                seqCD = seqC + '|' + seqD
                possiblePdifSiteSeqList.append(seqCD)
                possiblePdifSiteDict[i] = str(pos + 1) + '|' + str(pos + 28) + '|' + seqC + '|' + seqD + '|' + seqSpacer
                i = i + 1
                with open(possiblePdifSiteOutfile, 'a') as w:
                    w.write(str(pos + 1) + ' ' + str(pos + 28) + '\t' + seqC + '\t' + seqSpacer + '\t' + seqD + '\n')
                newPosList1.append(pos)
            print('seed search completed! %s match results' % (len(newPosList1)))
        else:
            newPosList1 = []
        if len(newPosList1) > 1:
            pairList = findPossiblePair(possiblePdifSiteSeqList, outdir)
            print('match results pair completed!')
            if pairList != []:
                pairList1 = []
                for i in range(0, len(pairList), 2):
                    left = pairList[i]
                    right = pairList[i + 1]
                    data = str(left) + '|' + str(right)
                    pairList1.append(data)
                pairListSet1 = set(pairList1)

                pairList2 = list(pairListSet1)
                pairList = []
                for pos in pairList2:
                    left = pos.split('|')[0]
                    right = pos.split('|')[1]
                    pairList.append(int(left))
                    pairList.append(int(right))

                with open(possiblePdifPairOutfileTemp, 'w') as w:
                    w.write('')
                index = 0
                pdifSiteList = []
                for i in range(0, len(pairList), 2):
                    left = pairList[i]
                    right = pairList[i + 1]
                    if left >= 0 and right >= 0:
                        start1 = int(possiblePdifSiteDict[left].split('|')[0])
                        end1 = int(possiblePdifSiteDict[left].split('|')[1])
                        start2 = int(possiblePdifSiteDict[right].split('|')[0])
                        end2 = int(possiblePdifSiteDict[right].split('|')[1])
                        if end1 < end2:
                            start = end1
                            end = start2
                        else:
                            start = end2
                            end = start1

                        # checkResult = checkResistanceGenePos(start1,end1,start2,end2,start,end,resistanceGenePosList)
                        checkResult = True
                        if checkResult:
                            index = index + 1
                            if left not in pdifSiteList:
                                pdifSiteList.append(left)
                            if right not in pdifSiteList:
                                pdifSiteList.append(right)
                            with open(possiblePdifPairOutfileTemp, 'a') as w:
                                s11 = possiblePdifSiteDict[left].split('|')[2]
                                s12 = possiblePdifSiteDict[left].split('|')[3]
                                s13 = possiblePdifSiteDict[left].split('|')[4]
                                s21 = possiblePdifSiteDict[right].split('|')[2]
                                s22 = possiblePdifSiteDict[right].split('|')[3]
                                s23 = possiblePdifSiteDict[right].split('|')[4]
                                w.write(' '.join(possiblePdifSiteDict[left].split('|')[:2]) + ' ' + ' '.join(
                                    possiblePdifSiteDict[right].split('|')[
                                    :2]) + ' ' + s11 + s12 + s13 + ' ' + s21 + s22 + s23 + '\n')
                        else:
                            continue
                tagDict = {}
                for l, data in enumerate(sorted(pdifSiteList)):
                    content = possiblePdifSiteDict[data]
                    contentText = '\t'.join(content.split('|')[:3]) + '\t' + content.split('|')[4] + '\t' + \
                                  content.split('|')[3]

                    pdifName = 'pdif%s' % (l + 1)
                    with open(pdifSiteOufile, 'a') as w:
                        w.write(pdifName + '\t' + contentText + '\n')
                    key = ' '.join(content.split('|')[:2])
                    tagDict[key] = pdifName
                f = open(possiblePdifPairOutfileTemp)
                with open(possiblePdifPairOutfileTemp2, 'w') as w:
                    w.write('')
                for lines in f:
                    lines = lines.strip()
                    if lines:
                        line = lines.split(' ')
                        key1 = ' '.join(line[:2])
                        key2 = ' '.join(line[2:4])
                        value1 = tagDict[key1]
                        value2 = tagDict[key2]
                        with open(possiblePdifPairOutfileTemp2, 'a') as w:
                            w.write(value1 + ' ' + value2 + ' ' + lines + '\n')
                f.close()
                if os.path.exists(possiblePdifPairOutfileTemp2):
                    f = open(possiblePdifPairOutfileTemp2)
                    contents = f.readlines()
                    f.close()
                    sort_contents = sorted(contents,
                                           key=lambda data: (int(data.split(' ')[2]), int(data.split(' ')[4])))
                    for cont in sort_contents:
                        with open(possiblePdifPairOutfile, 'a') as w:
                            w.write(cont)

                print('%s pdif sites found' % len(pdifSiteList))
            else:
                print('No pdif sites found')
        return True
    else:
        print("No featured sequences")
        return False


def findMatchFragmentThread(num4, pos0, outdir, name, seedList, seq, start1, end1):
    outfile = outdir + '/tmp/seedSearch/' + str(num4) + '/%s' % name
    start = (start1 - 1) * 2
    end = (end1 - 1) * 2 + 2
    for k in range(start, end, 2):
        seqc = seedList[k]
        seqd = seedList[k + 1]
        posList = []
        maxMistachXerC = 3
        maxMistachXerD = 2
        for i in range(1):
            initPos = i
            endPos = i + 11
            tempSeqC = seq[initPos:endPos]

            tempMistachNumberC = 0
            tempMistachNumberD = 0
            tempMistachNumberCLeft = 0
            for l in range(5):
                if tempSeqC[l] != seqc[l]:
                    tempMistachNumberCLeft = tempMistachNumberCLeft + 1
            if tempMistachNumberCLeft > maxMistachXerC:
                tempMistachNumberC = -1
                continue
            else:
                for j in range(11):
                    if tempSeqC[j] != seqc[j]:
                        tempMistachNumberC = tempMistachNumberC + 1

            if tempMistachNumberC > maxMistachXerC:
                continue
            else:
                if (endPos + 17) <= len(seq):
                    tempSeqD = seq[initPos + 17:endPos + 17]
                    for k in range(11):
                        if tempSeqD[k] != seqd[k]:
                            tempMistachNumberD = tempMistachNumberD + 1
                    if tempMistachNumberD > maxMistachXerD:
                        continue
                    else:
                        with open(outfile, 'a') as w:
                            w.write(str(initPos + pos0) + '\n')


def checkResistanceGenePos(start1, end1, start2, end2, start, end, resistanceGenePosList):
    containFlag = False
    intersectionFlag = False
    for pos in resistanceGenePosList:
        start0 = int(pos.split('-')[0])
        end0 = int(pos.split('-')[1])
        if start <= start0 and end >= end0:
            containFlag = True
            break
    for pos in resistanceGenePosList:
        start3 = int(pos.split('-')[0])
        end3 = int(pos.split('-')[1])
        if start3 <= end1 and end3 >= start1:
            intersectionFlag = True
            break
        if start3 <= end2 and end3 >= start2:
            intersectionFlag = True
            break
    if containFlag == True and intersectionFlag == False:
        return True
    else:
        return False


def findPossiblePair(possiblePdifSiteSeqList, outdir):
    pairList = []
    length = len(possiblePdifSiteSeqList)

    batchList = []
    length = len(possiblePdifSiteSeqList)
    numList = []
    for i in range(length):
        numList.append(i)
    while 1:
        if len(numList) > 1:
            num = numList.pop()
            for number in numList:
                combine = str(num) + '-' + str(number)
                batchList.append(combine)
        else:
            break

    batchLength = len(batchList)
    batch = 5
    if batchLength < batch:
        batch = batchLength
    batchNumber = math.floor(batchLength / batch)
    restNumber = batchLength - 50 * batchNumber
    threadList = []
    for i in range(1, batch + 1, 1):
        start = (i - 1) * batchNumber
        end = start + batchNumber - 1
        p = Thread(target=findPossiblePairThread, args=(outdir, i, possiblePdifSiteSeqList, batchList, start, end,))
        threadList.append(p)
    p1 = Thread(target=findPossiblePairThread, args=(
        outdir, batch + 1, possiblePdifSiteSeqList, batchList, batchLength - restNumber, batchLength - 1,))
    threadList.append(p1)
    for thread in threadList:
        thread.start()
    for thread in threadList:
        thread.join()
    for file in os.listdir(outdir + '/tmp/pairSearch'):
        f = open(outdir + '/tmp/pairSearch/' + file)
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            start = int(line.split(' ')[0])
            end = int(line.split(' ')[1])
            pairList.append(start)
            pairList.append(end)
        f.close()
    return pairList


def findPossiblePairThread(outdir, name, possiblePdifSiteSeqList, batchList, start, end):
    pairList = []
    maxMistachCD = 4
    maxMistachDC = 6
    for i in range(start, end + 1):

        pos1 = int(batchList[i].split('-')[0])
        pos2 = int(batchList[i].split('-')[1])
        seqC1 = possiblePdifSiteSeqList[pos1].split('|')[0]
        seqD1 = possiblePdifSiteSeqList[pos1].split('|')[1]
        reverCompleteSeqC1 = reverComplement(seqC1)
        reverCompleteSeqD1 = reverComplement(seqD1)
        seqC2 = possiblePdifSiteSeqList[pos2].split('|')[0]
        seqD2 = possiblePdifSiteSeqList[pos2].split('|')[1]
        misMatch1 = compareTwoSeq(reverCompleteSeqC1, seqD2)
        misMatch2 = compareTwoSeq(reverCompleteSeqD1, seqC2)
        misMatch3 = compareTwoSeq(reverCompleteSeqD1, seqC2)
        misMatch4 = compareTwoSeq(reverCompleteSeqC1, seqD2)
        if misMatch1 <= maxMistachCD and misMatch2 <= maxMistachDC:
            pairList.append(pos1)
            pairList.append(pos2)
        if misMatch3 <= maxMistachDC and misMatch4 <= maxMistachCD:
            pairList.append(pos1)
            pairList.append(pos2)

    for i in range(0, len(pairList), 2):
        with open(outdir + '/tmp/pairSearch/%s' % name, 'a') as w:
            w.write(str(pairList[i]) + ' ' + str(pairList[i + 1]) + '\n')


def reverComplement(seq):
    mySeq = Seq(seq)
    retSeq = str(mySeq.reverse_complement())
    return retSeq


def compareTwoSeq(seq1, seq2):
    misMatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            misMatch = misMatch + 1
    return misMatch


def findResistanceGene(inFile, outdir, blastnPath, scriptsDir):
    resistanceGenePosList = []
    database = "ncbi"
    default_ident = 90
    default_cov = 60
    cmd = blastnPath + ' -task blastn -dust no -evalue 1E-20 -culling_limit 1 -query ' + inFile + ' -db ' + scriptsDir + '/AMRDB/sequences -outfmt \"6 sseqid pident length  qstart qend slen sstrand\" -out ' + outdir + '/tmp/blast.out'
    subprocess.run(cmd, shell=True)
    blastResultDict = {}
    f = open(outdir + '/tmp/blast.out')
    lines = f.readlines()
    if len(lines) == 0:
        pass
    else:
        for line in lines:
            data = line.strip().split('\t')
            name = re.split('~~~', data[0])[1]
            ident = float(data[1])
            cov = float(data[2]) * 100 / int(data[5])
            leftPos = data[3]
            rightPos = data[4]
            strand = data[6]
            data = str(ident) + '|' + str(cov) + '|' + leftPos + '|' + rightPos + '|' + name + '|' + strand
            if name not in blastResultDict:
                blastResultDict[name] = []
                blastResultDict[name].append(data)
            else:
                blastResultDict[name].append(data)
    f.close()
    for name, dataList in blastResultDict.items():
        maxIdent = default_ident
        maxCov = default_cov
        pos = ''
        for data in dataList:
            ident = float(data.split('|')[0])
            cov = float(data.split('|')[1])
            leftPos = data.split('|')[2]
            rightPos = data.split('|')[3]
            name = data.split('|')[4]
            strand = data.split('|')[5]
            if ident >= maxIdent and cov >= maxCov:
                bestLeftPos = leftPos
                bestRightPos = rightPos
                pos = bestLeftPos + '-' + bestRightPos
                bestStand = strand
                bestIdent = format(ident, '0.2f')
                bestCov = format(cov, '0.2f')

        if pos != '':
            left1 = int(pos.split('-')[0])
            right1 = int(pos.split('-')[1])
            flag = True
            if resistanceGenePosList != []:
                for po in resistanceGenePosList:
                    left2 = int(po.split('-')[0])
                    right2 = int(po.split('-')[1])
                    if left1 < right2 and right1 > left2:
                        flag = False
                        break
                    else:
                        continue
            if flag:
                with open(outdir + '/AMRgene1.txt', 'a') as w:
                    w.write(name + '\t' + bestLeftPos + '\t' + bestRightPos + '\t' + bestStand + '\t' + str(
                        bestIdent) + '\t' + str(bestCov) + '\n')
                resistanceGenePosList.append(pos)
    return resistanceGenePosList


def getSeqFromGenbankFile(genbankFile, outdir):
    outfile = outdir + '/inputFile.fasta'
    with open(outfile, 'w') as w:
        w.write('')
    gb_seqs = SeqIO.parse(genbankFile, "gb")
    for gb_seq in gb_seqs:
        seq = str(gb_seq.seq)
        for key in gb_seq.features:
            if key.type == "source":
                if 'plasmid' in key.qualifiers:
                    id = key.qualifiers['plasmid'][0]
                    if len(id) != 0:
                        with open(outfile, 'a') as w:
                            w.write('>' + id + '\n' + seq + '\n')
                    else:
                        id = gb_seq.id
                        with open(outfile, 'a') as w:
                            w.write('>' + id + '\n' + seq + '\n')
                else:
                    id = gb_seq.id
                    with open(outfile, 'a') as w:
                        w.write('>' + id + '\n' + seq + '\n')
            break
    return outfile


def getSeqFromFastaFile(inFile, outdir):
    outfile = outdir + '/inputFile.fasta'
    for rec in SeqIO.parse(inFile, 'fasta'):
        print("Sequence length:", len(str(rec.seq)))
        SeqIO.write(rec, outfile, "fasta")
        break
    return outfile


def finalFilter(outdir):
    pdifSiteFile = outdir + '/pdif_site1.txt'
    pdifDict = {}
    f = open(pdifSiteFile)
    for lines in f:
        lines = lines.strip()  # 删掉换行符
        if lines:
            data = lines.split('\t')  # 遇到Tab隔开
            pdifDict[data[0]] = lines
    f.close()
    os.remove(pdifSiteFile)

    pairFile = outdir + '/pdif_pair.txt'
    pairList = []
    numList = []
    f = open(pairFile)
    for lines in f:
        lines = lines.strip()
        if lines:
            data = lines.split(' ')
            name1 = data[0]
            name2 = data[1]
            num1 = int(name1.replace('pdif', ''))
            num2 = int(name2.replace('pdif', ''))
            if num1 > num2:
                num1, num2 = num2, num1
            if num1 not in numList:
                numList.append(num1)
            if num2 not in numList:
                numList.append(num2)
            info = str(num1) + '-' + str(num2)
            pairList.append(info)
    f.close()
    os.remove(pairFile)
    index = 1
    completedList = []
    numList = sorted(numList)
    for num in numList:
        if num not in completedList:
            num1 = num + 1
            info1 = str(num) + '-' + str(num1)
            if info1 in pairList:
                completedList.append(num)
                completedList.append(num1)
                key1 = 'pdif' + str(num)
                key2 = 'pdif' + str(num1)

                newName1 = 'pdif' + str(index)
                newName2 = 'pdif' + str(index + 1)
                index = index + 2
                with open(pdifSiteFile, 'a') as w:
                    w.write(newName1 + '\t' + '\t'.join(pdifDict[key1].split('\t')[1:]) + '\tC|D\n')
                    w.write(newName2 + '\t' + '\t'.join(pdifDict[key2].split('\t')[1:]) + '\tD|C\n')



def singleThread(inFile, outdir, blastnPath, scriptsDir, circle):
    print("%s: analysis start" % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
    format = checkInfileFormat(inFile)
    pdifDB = scriptsDir + '/data/redundant.seed.fa'
    if format:
        if blastnPath:
            resistanceGenePosList = findResistanceGene(inFile, outdir, blastnPath, scriptsDir)
            print("%s resistance genes found" % len(resistanceGenePosList))
            startTime1 = time.time()
            if resistanceGenePosList != []:
                make_sort(outdir + '/AMRgene1.txt', "amr")
                flag = findPdif(inFile, outdir, blastnPath, pdifDB, resistanceGenePosList)
                if flag:
                    if os.path.exists(outdir + '/pdif_pair.txt'):
                        finalFilter(outdir)
                if os.path.exists(outdir + '/pdif_site1.txt'):
                    if os.path.getsize(outdir + '/pdif_site1.txt'):
                        # ISfinder(inFile, outdir, blastnPath, scriptsDir)
                        # make_sort(outdir + '/ISfinder.filter.xls', "is")
                        if circle == "true":
                            angularPlasmid(outdir)
                            # with open(outdir + '/plasmid.txt', 'w') as w:
                            #     w.write(plasmid)
                        # if circle == "false":
                        #     html_file = os.path.join(outdir, "structure.html")
                        #     initEnd, text = echartsSequence(outdir, html_file, scriptsDir)
                        #     with open(outdir + '/e1.txt', 'w') as w:
                        #         w.write(str(initEnd))
                        #     with open(outdir + '/e2.txt', 'w') as w:
                        #         w.write(text.replace("\\", ""))
                        changepdifname(inFile, outdir, scriptsDir)
                        Getpdifmoduleseq(inFile, outdir)
                        if os.path.exists(outdir + '/pdifmodule_list.txt'):
                            pdifmoduleechart(outdir, scriptsDir)
                            print("%s: analysis end\n" % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
                        else:
                            print("No pdif module found!")
            else:
                print("No pdif site found because of lack of resistance gene!")
        else:
            print("Program terminated because lack of dependencies!")
    else:
        print("Please check your input file because of invalid character!")
    shutil.rmtree("%s/tmp" % outdir)
    # shutil.rmtree("%s/tmp" % outdir)
    # os.system("rm -r %s/tmp"%outdir)




def getFileFormat(infile):
    try:
        for rec in SeqIO.parse(infile, 'fasta'):
            seq = rec.seq
            if seq:
                return "fasta"
        for rec in SeqIO.parse(infile, 'genbank'):
            seq = rec.seq
            # print(seq)
            if seq:
                return "genbank"
    except:
        return "unknown"
    return "unknown"


def makeOutdir(outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        os.mkdir(outdir + '/tmp')
        os.mkdir(outdir + '/tmp/pairSearch')
        os.mkdir(outdir + '/tmp/seedSearch')
    else:
        shutil.rmtree(outdir)
        os.mkdir(outdir)
        os.mkdir(outdir + '/tmp')
        os.mkdir(outdir + '/tmp/pairSearch')
        os.mkdir(outdir + '/tmp/seedSearch')


def make_sort(infile, sign):
    if sign == "amr":
        outfile = infile + ".temp"
        f = open(infile)
        datas = f.readlines()
        sorted_datas = sorted(datas, key=lambda data: int(data.split("\t")[1]))
        f.close()
        with open(outfile, 'w') as w:
            w.write(''.join(sorted_datas))
        os.remove(infile)
        os.rename(outfile, infile)
    elif sign == "is":
        outfile = infile + ".temp"
        f = open(infile)
        datas = f.readlines()
        sorted_datas = sorted(datas[1:], key=lambda data: int(data.split("\t")[6]))
        f.close()
        with open(outfile, 'w') as w:
            w.write(datas[0])
            w.write(''.join(sorted_datas))
        os.remove(infile)
        # os.remove(outfile)
        os.rename(outfile, infile)



def compare(str1, scriptDir):
    pdifdatabase = scriptDir + '/data/pdifdatabase.fasta'
    with open(pdifdatabase, 'r') as file:
        pdifdblines = file.readlines()
        for line1 in pdifdblines:
            pdifdblist = line1.strip().split('\t')
            # pdiftype = re.split('\t', pdifdblist[0])[0]
            pdifname2 = pdifdblist[1]
            pdif = pdifdblist[2]
            if str1 == pdif:
                return pdifname2




def changepdifname(inFile, outdir, scriptDir):
    site1 = outdir + '/pdif_site1.txt'
    site = outdir + '/pdif_site.txt'
    pdifdatabase = scriptDir + '/data/pdifdatabase.fasta'
    AMRgene1 = outdir + '/AMRgene1.txt'
    AMRgene = outdir + '/AMRgene.txt'
    # file1 = dir1 + '\\3.txt'
    # file2 = dir1 + '\\4.txt'
    for seq_rec in SeqIO.parse(inFile, "fasta"):
        with open(AMRgene1, 'r') as a:
            AMRgeneline = a.readlines()
            for line2 in AMRgeneline:
                with open(AMRgene, 'a') as b:
                    b.write(seq_rec.id + '\t' + line2)
        with open(site1, 'r') as f:
            pdiflines = f.readlines()
            if len(pdiflines) == 0:
                pass
            else:
                for line in pdiflines:
                    pdifsitelist = line.strip().split('\t')
                    # pdifname1 = re.split('\t', pdifsitelist[0])[0]
                    pdifseq = pdifsitelist[3] + pdifsitelist[4] + pdifsitelist[5]
                    pdifseq2 = str()
                    if pdifsitelist[6] == 'C|D':
                        pdifseq2 = pdifseq
                    elif pdifsitelist[6] == 'D|C':
                        pdifseq1 = Seq(pdifseq)
                        pdifseq2 = str(pdifseq1.reverse_complement())
                    flag = compare(pdifseq2, scriptDir)
                    if flag:
                        newline = seq_rec.id + '\t' + line.replace('\n', '') + '\t' + flag + '\n'
                        with open(site, 'a') as w:
                            w.write(newline)
                    else:
                        with open(pdifdatabase, 'r') as file:
                            pdifdblist = file.readlines()
                            num = len(pdifdblist)
                        with open(pdifdatabase, 'a') as file:
                            file.write('[single]' + '\t' + 'pdif' + str(num) + '\t' + pdifseq2 + '\t' + seq_rec.id + '\n')
                            newline = seq_rec.id + '\t' + line.replace('\n', '') + '\t' + 'pdif' + str(len(pdifdblist)) + '\n'
                            with open(site, 'a') as w:
                                w.write(newline)
    os.remove(site1)
    os.remove(AMRgene1)

def Getpdifmoduleseq(inFile, outdir):
    AMRlist = outdir + '/AMRgene.txt'
    site = outdir + '/pdif_site.txt'
    pdifmoduleseq = outdir + '/pdifmoduleseq.fasta'
    pdifmodulelist = outdir + '/pdifmodule_list.txt'
    pdifmodulelist1 = outdir + '/pdifmodule_list1.txt'
    comparefile = outdir + '/genelist.txt'
    for seq_rec in SeqIO.parse(inFile, "fasta"):
        with open(AMRlist, 'r') as file:
            Genelines = file.readlines()
            with open(site, 'r') as f:
                pdiflines = f.readlines()
                if len(Genelines) == 0 or len(pdiflines) == 0:
                    pass
                else:
                    for i in range(len(pdiflines) - 1):
                        line1 = pdiflines[i]
                        upsteampdifsitelist = line1.strip().split('\t')
                        upsteampdifname = upsteampdifsitelist[-1]
                        upsteampdifstartPos = upsteampdifsitelist[2]
                        upsteampdifendPos = upsteampdifsitelist[3]
                        upsteampdifCD = upsteampdifsitelist[7]
                        upsteampdif = upsteampdifsitelist[3] + upsteampdifsitelist[4] + upsteampdifsitelist[5]
                        line2 = pdiflines[i + 1]
                        downsteampdifsitelist = line2.strip().split('\t')
                        downsteampdifname = downsteampdifsitelist[-1]
                        downsteampdifstartPos = downsteampdifsitelist[2]
                        downsteampdifendPos = downsteampdifsitelist[3]
                        downsteampdifCD = downsteampdifsitelist[7]
                        downsteampdif = downsteampdifsitelist[3] + downsteampdifsitelist[4] + downsteampdifsitelist[5]
                        ARGname = str()
                        ARGnamen = str()
                        ARGPos = str()
                        dis11 = []
                        dis21 = []
                        for line in Genelines:
                            AMRgenelist = line.strip().split('\t')
                            AMRname1 = re.split('\t', AMRgenelist[1])[0]
                            AMRstartPos = AMRgenelist[2]
                            AMRendPos = AMRgenelist[3]
                            if int(AMRstartPos) > int(upsteampdifendPos) and int(AMRendPos) < int(
                                    downsteampdifstartPos):
                                ARGname = ARGname + '-' + AMRname1
                                ARGnamen = ARGnamen + '+' + AMRname1
                                ARGPos1 = AMRstartPos + '-' + AMRendPos
                                ARGPos = ARGPos + '\t' + ARGPos1
                                dis11.append(int(AMRstartPos) - int(upsteampdifendPos) - 1)
                                dis21.append(int(downsteampdifstartPos) - int(AMRendPos) - 1)
                                with open(comparefile, 'a') as a:
                                    a.write(line)
                        if dis11 and dis21:
                            dis1 = min(dis11)
                            dis2 = min(dis21)
                            seqlength = int(downsteampdifendPos) - int(upsteampdifstartPos)
                            myseq = str(seq_rec.seq)[int(upsteampdifstartPos) - 1:int(downsteampdifendPos)]
                            with open(pdifmoduleseq, 'a') as w:
                                w.write(
                                    '>' + seq_rec.id + '\t' + upsteampdifname + ARGname + '+' + downsteampdifname + '\n' + myseq + '\n')
                            with open(pdifmodulelist, 'a') as w:
                                w.write(
                                    '>' + seq_rec.id + '\t' + upsteampdifname + ARGname + '-' + downsteampdifname + '\t' + upsteampdifCD + '\t' + upsteampdifstartPos + '\t' + str(
                                        dis1) + ARGPos + '\t' + str(dis2) + '\t' + downsteampdifendPos + '\t' + downsteampdifCD + '\n')
                            with open(pdifmodulelist1, 'a') as w:
                                w.write(
                                    '>' + seq_rec.id + '\t' + upsteampdifname + ARGnamen + '+' + downsteampdifname + '\t' + upsteampdifCD + '\t' + upsteampdifstartPos + '\t' + str(
                                        dis1) + ARGPos + '\t' + str(dis2) + '\t' + downsteampdifendPos + '\t' + downsteampdifCD + '\n')
                        else:
                            pass


def split(filename, outdir):
    os.mkdir(outdir + '/allfile')
    indir = outdir + '/allfile'
    formatin = getFileFormat(filename)
    outfile = []
    if formatin == "fasta":
        with open(filename, 'r') as f1:
            lines = f1.readlines()
            for line in lines:
                if line.startswith(">"):
                    if (outfile != []): outfile.close()
                    names = line.split(' ')[0][1:]
                    name = names.split('\n')[0]
                    filename1 = indir + '/' + name + ".fasta"
                    outfile = open(filename1, 'w')
                    outfile.write(line)
                else:
                    outfile.write(line)
        outfile.close()

def process(outdir, blastnPath, scriptsDir, circle):
    indir1 = outdir + '/allfile'
    if indir1 != None:
        for file in os.listdir(indir1):
            infile1 = os.path.join(indir1, file)
            name = ".".join(file.strip().split('.')[0:-1])
            outdir1 = outdir + '/' + name
            makeOutdir(outdir1)
            formatin = ''
            inFile = ''
            formatin = getFileFormat(infile1)
            if formatin == 'fasta':
                inFile = getSeqFromFastaFile(infile1, outdir1)
                os.remove(infile1)
            elif formatin == 'genbank':
                inFile = getSeqFromGenbankFile(infile1, outdir1)
                os.remove(infile1)
            singleThread(inFile, outdir1, blastnPath, scriptsDir, circle)
    shutil.rmtree(indir1)

def main():
    # print("%s: analysis start"%(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())))
    startTime = time.time()
    # script abspath
    scriptsDir = sys.path[0] + '/../PdifFinder'
    fastaFile, genbankFile, indir, outdir, circle = arg_parse()
    blastnPath = check_dependencies()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        shutil.rmtree(outdir)
        os.mkdir(outdir)
    if genbankFile:
        print("current file format: ", "genbank")
        inFile = getSeqFromGenbankFile(genbankFile, outdir)
        split(inFile, outdir)
        os.remove(inFile)
        process(outdir, blastnPath, scriptsDir, circle)
    elif fastaFile != None:
        print("current file format: ", "fasta")
        split(fastaFile, outdir)
        process(outdir, blastnPath, scriptsDir, circle)
    elif indir != None:
        for file in os.listdir(indir):
            infile1 = os.path.join(indir, file)
            name = file.split('.')[0]
            outdir1 = outdir + '/' + time.strftime('%Y%m%d%H%M%S', time.localtime()) + '_' + name
            os.mkdir(outdir1)
            format = ""
            format = getFileFormat(infile1)
            print("current file format: ", format)
            if format == "genbank":
                inFile = getSeqFromGenbankFile(infile1, outdir)
                split(inFile, outdir1)
                os.remove(inFile)
                process(outdir1, blastnPath, scriptsDir, circle)
            elif format == "fasta":
                split(infile1, outdir1)
                process(outdir1, blastnPath, scriptsDir, circle)
            else:
                shutil.rmtree(outdir1)
                print("please check this file format")

    endTime = time.time()
    print('All jobs finished, consumed %ss' % (endTime - startTime))
    # print("%s: analysis end"%(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())))


if __name__ == "__main__":
    main()
