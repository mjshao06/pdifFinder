from Bio import SeqIO
import sys,os,re
from collections import OrderedDict


def read_file(outdir):
    pdifSiteFile = outdir + '/pdif_site.txt'
    resGeneFile = outdir + '/AMRgene.txt'
    # isFile = outdir + '/ISfinder.filter.xls'
    posDict = OrderedDict()
    posList = []
    index = 0
    if os.path.exists(pdifSiteFile) and os.path.getsize(pdifSiteFile) != 0:
        f = open(pdifSiteFile)
        for lines in f:
            lines = lines.strip()
            if lines:
                data = lines.split('\t')
                start = data[1]
                end = data[2]
                strand = '0'
                name = data[0]
                index = index + 1
                posDict[index] = start+'|'+end+'|'+strand+'|'+name
                posList.append(int(start))
                posList.append(int(end))
        f.close()
    if os.path.exists(resGeneFile) and os.path.getsize(resGeneFile) != 0:
        f = open(resGeneFile)
        for lines in f:
            lines = lines.strip()
            if lines:
                data = lines.split('\t')
                start = data[1]
                end = data[2]
                strandSign = data[3]
                if strandSign == 'plus':
                    strand = '1'
                else:
                    strand = '-1'
                name = data[0]
                index = index + 1
                posDict[index] = start+'|'+end+'|'+strand+'|'+name
                posList.append(int(start))
                posList.append(int(end))
        f.close()
    # if os.path.exists(isFile) and os.path.getsize(isFile) != 0:
    #     n = 0
    #     f = open(isFile)
    #     for lines in f:
    #         n = n + 1
    #         if n == 1:
    #             continue
    #         lines = lines.strip()
    #         if lines:
    #             data = lines.split('\t')
    #             start = data[6]
    #             end = data[7]
    #             s_start = int(data[8])
    #             s_end = int(data[9])
    #             if s_start > s_end:
    #                 strand = '-2'
    #             else:
    #                 strand = '2'
    #             name = data[1].split('|')[0]
    #             index = index + 1
    #             posDict[index] = start+'|'+end+'|'+strand+'|'+name
    #             posList.append(int(start))
    #             posList.append(int(end))
    #     f.close()
    initStart = sorted(posList)[0]
    initEnd = sorted(posList)[-1]
    return posDict,initStart,initEnd
    
def echartsSequence(indir,outfile,scriptsDir):
    if os.path.exists(indir + '/pdif_site.txt'):
        posDict,initStart,initEnd = read_file(indir)
        text = ''
        for k,v in posDict.items():
            name1 = v.split('|')[3].replace("(","\(").replace(")","\)").replace("-","\-")
            name = name1.replace("'","\\'")
            start0 = v.split('|')[0]
            end0 = v.split('|')[1]
            strand = v.split('|')[2]
            position = "top"
            if strand == '0':
                color = 'green'
                start = start0
                end = end0
                position = "bottom"
            elif strand == '1':
                color = 'red'
                start = start0
                end = end0
            elif strand == "-1":
                color = 'red'
                start = end0
                end = start0
            elif strand == '2':
                strand = '1'
                color = 'blue'
                start = start0
                end = end0
            elif strand == '-2':
                strand = '-1'
                color = 'blue'
                start = end0
                end = start0
            text = text + '\{name:\\"%s\\"'%name+',value:\[0,%s,%s,%s\],itemStyle:\{color:\\"%s\\"\},label:\{color:\\"black\\",rotate:45,show:true,distance:25,position:\\"%s\\",textStyle:\{fontSize:12\},formatter:function\(params\)\{return params\.name\;\}\}\},'%(str(start),str(end),strand,color,position)
        
        os.system("cp %s/demo.html %s"%(scriptsDir,outfile))
        os.system("sed -i \"s#var data = \[\]#var data = \[%s\]#\" %s"%(text,outfile))
        os.system("sed -i 's/var endpos = 100000/var endpos = %s/' %s"%(initEnd,outfile))
        return initEnd,text
        #tmpfile = outfile.replace('structure.html','tmp.txt')
        #with open(tmpfile,'a') as w:
            #w.write('0\n%s\n%s'%(initEnd,text.replace("\\","")))

        

    