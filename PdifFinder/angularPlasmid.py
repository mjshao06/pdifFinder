import os
import sys
import json
from Bio import SeqIO
import math
from collections import OrderedDict
import copy

def writeHtml(outfile,content):
    with open(outfile,'a') as w:
        w.write(content+'\n')
        
def read_file(outdir):
    pdifSiteFile = outdir + '/pdif_site1.txt'
    resGeneFile = outdir + '/AMRgene1.txt'
    # isFile = outdir + '/ISfinder.filter.xls'
    infoDict = {}
    infoDict4 = {}
    numList = []
    j = 0
    # f = open(isFile)
    # for lines in f:
    #     lines = lines.strip()
    #     j = j + 1
    #     if j == 1:
    #         continue
    #     if lines:
    #         data = lines.split('\t')
    #         name = data[1].split('|')[0]
    #         if not name:
    #             name = data[1].split('|')[1]
    #         start = data[6]
    #         end = data[7]
    #         sstart = int(data[8])
    #         send = int(data[9])
    #         if sstart < send:
    #             strand = "1"
    #         else:
    #             strand = "-1"
    #         info = name + '|' + start+'|'+end+'|'+strand + '|IS'
    #         infoDict4[int(start)] = info
    #         numList.append(int(start))
    #         numList.append(int(end))
    # f.close()
    f = open(pdifSiteFile)
    for lines in f:
        lines = lines.strip()
        if lines:
            data = lines.split('\t')
            start = data[1]
            end = data[2]
            seq = data[3]
            name = data[-1].replace('|','~')
            strand = '0'
            info = name + '|' + start+'|'+end+'|'+strand + '|pdif'
            infoDict[int(start)] = info
            numList.append(int(start))
            numList.append(int(end))
    f.close()
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
            info = name + '|' + start+'|'+end+'|'+strand + '|resGene'
            infoDict[int(start)] = info
            numList.append(int(start))
            numList.append(int(end))
    f.close()
    max = sorted(numList)[-1]
    # os.remove(isFile)
    return infoDict,max
 
def angularPlasmid(outdir):
    plasmid = ""
    alpha = 0.05
    fastaFile = outdir+'/inputFile.fasta'
    for rec in SeqIO.parse(fastaFile,'fasta'):
        seq = str(rec.seq)
        length = len(seq)
        break
    infoDict,max = read_file(outdir)
    plasmidheight = "1200"
    plasmidwidth = "1200"
    SeqID = rec.id
    radius = 400
    if length >= 10000:
        interval1 = str(math.floor((length-5000)/5)/10)
        interval2 = str(math.floor((length-5000)/5))
    else:
        interval1 = str(math.floor((length-1000)/5)/10)
        interval2 = str(math.floor((length-1000)/5))
    #legend and angular init
    #header = '<div class="col-md-2" style="position:relative;left:70%">\n<svg>\n<rect x="60" y="100" width="50" height="15" style="fill:green;stroke:green;stroke-width:5;" />\n<text x="150" y="105" font-size="14" >Insertion element</text>\n<rect x="60" y="60" width="50" height="15" style="fill:#666;stroke:#666;stroke-width:5;" />\n<text x="150" y="75" font-size="14">pdif site</text>\n<rect x="60" y="20" width="50" height="15" style="fill:blue;stroke:blue;stroke-width:5;" />\n<text x="150" y="35" font-size="14">AMR gene</text>\n</svg>\n</div>\n<div  class="col-md-10" style="position:relative;right:10%;bottom:10%">\n'
    header = '<div class="col-md-2" style="position:relative;left:70%">\n<svg>\n<rect x="60" y="60" width="50" height="15" style="fill:#666;stroke:#666;stroke-width:5;" />\n<text x="150" y="75" font-size="14">pdif site</text>\n<rect x="60" y="20" width="50" height="15" style="fill:blue;stroke:blue;stroke-width:5;" />\n<text x="150" y="35" font-size="14">AMR gene</text>\n</svg>\n</div>\n<div  class="col-md-10" style="position:relative;right:10%;bottom:10%">\n'
    header1 =  '<div style="position:relative;left:50%;top:5%">\n<svg>\n<rect x="60" y="100" width="50" height="15" style="fill:green;stroke:green;stroke-width:5;" />\n<text x="150" y="105" font-size="14" >Insertion element</text>\n<rect x="60" y="60" width="50" height="15" style="fill:#666;stroke:#666;stroke-width:5;" />\n<text x="150" y="75" font-size="14">pdif site</text>\n<rect x="60" y="20" width="50" height="15" style="fill:blue;stroke:blue;stroke-width:5;" />\n<text x="150" y="35" font-size="14">AMR gene</text>\n</svg>\n</div>\n<div   style="position:relative;margin-top:-15%;right:10%">\n'
    outfile = outdir + '/plasmid.html'
    with open(outfile,'w') as w:
        content = """
            <html>

            <head>
                <style>
                .labelline {
                    stroke: black;
                    stroke-dasharray: 2, 2;
                    stroke-width: 2px;
                }
                </style>
                <script src="https://cdn.staticfile.org/jquery/2.1.1/jquery.min.js"></script>
                <script src="http://bacant.net/static/js/angularplasmid.complete.min.js"></script>
                <link href="https://cdn.staticfile.org/twitter-bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">
                <script src="https://cdn.staticfile.org/twitter-bootstrap/3.3.7/js/bootstrap.min.js"></script>
            </head>

            <body>
        """  
        w.write(content)
        w.write(header)
        plasmid = plasmid + header1
        content = '<plasmid sequencelength="{sequencelength}" plasmidheight="{plasmidheight}" plasmidwidth="{plasmidwidth}">\n<plasmidtrack width="1" trackstyle="fill:none;stroke:#999;" radius="{radius}">\n<tracklabel labelstyle="font-size:20px;font-weight:700;fill:#666;" text="{seqID}"></tracklabel>\n<tracklabel labelstyle="font-size:12px;font-weight:400;fill:#999;" text="{sequencelength}bp" vadjust="18"></tracklabel>\n<trackscale interval="{interval1}" style="stroke:#999;" vadjust="8"></trackscale>\n<trackscale interval="{interval2}" showlabels="1" style="stroke:#333;stroke-width:2px" ticksize="5" vadjust="8" labelstyle="font-size:9px;fill:#999;font-weight:300" labelvadjust="15">\n</trackscale>\n'.format(sequencelength=length,plasmidheight=plasmidheight,plasmidwidth=plasmidwidth,radius=radius,interval1=interval1,interval2=interval2,seqID=SeqID)
        w.write(content)
        plasmid = plasmid + content
    hadjust = '0'
    vadjust = '120'
    valign = "outer"
    linevadjust = "-20"
    #resistance gene and pdif site
    for info in infoDict.values():
        data = info.split('|')
        name = data[0].replace('~','|')
        start = data[1]
        end = data[2]
        strand = data[3]
        flag = data[4]
        length1 = abs(int(end)-int(start)) + 1   
        if flag == "resGene":
            color = "blue"
            archiveLen = (length1/length)*2*math.pi*radius*alpha
            if archiveLen <= 10:
                arrowendlength = str(archiveLen)
                arrowendwidth = "5"
                arrowstartlength = str(archiveLen)
                arrowstartwidth = "5"
            else:
                arrowendlength = "10"
                arrowendwidth = "5"
                arrowstartlength = "10"
                arrowstartwidth = "5"
            if strand == "1":
                content = '<!-- {name} -->\n<trackmarker start="{start}" end="{end}" markerstyle="stroke:#000;fill:{color};" arrowendlength="{arrowendlength}" arrowendwidth="{arrowendwidth}" wadjust="10" vadjust="-5">\n<markerlabel labelstyle="font-size:10px" text="{name}" vadjust="{vadjust}" hadjust="{hadjust}" valign="{valign}" showline="1" linevadjust="{linevadjust}"  lineclass="labelline"></markerlabel>\n</trackmarker>\n'.format(name=name,start=start,end=end,arrowendlength=arrowendlength,arrowendwidth=arrowendwidth,vadjust=vadjust,hadjust=hadjust,color=color,valign=valign,linevadjust=linevadjust)
            elif strand == "-1":
                content = '<!-- {name} -->\n<trackmarker start="{start}" end="{end}" markerstyle="stroke:#000;fill:{color};" arrowstartlength="{arrowstartlength}" arrowstartwidth="{arrowstartwidth}" wadjust="10" vadjust="-5">\n<markerlabel labelstyle="font-size:10px" text="{name}" vadjust="{vadjust}" hadjust="{hadjust}" valign="{valign}" showline="1" linevadjust="{linevadjust}"  lineclass="labelline"></markerlabel>\n</trackmarker>\n'.format(name=name,start=start,end=end,arrowstartlength=arrowstartlength,arrowstartwidth=arrowstartwidth,vadjust=vadjust,hadjust=hadjust,color=color,valign=valign,linevadjust=linevadjust)
        else:
            content = '<!-- {name} -->\n<trackmarker start="{start}" end="{end}" markerstyle="stroke:#666;fill:#666;" wadjust="20" vadjust="-10">\n<markerlabel labelstyle="font-size:11px"  text="{name}" vadjust="{vadjust}" hadjust="{hadjust}" valign="{valign}" linevadjust="{linevadjust}" showline="1"  lineclass="labelline"></markerlabel>\n</trackmarker>\n'.format(name=name,start=start,end=end,vadjust="-60",hadjust="hadjust",valign="inner",linevadjust="20")
        writeHtml(outfile,content)
        plasmid = plasmid + content
    #IS
    # vadjust = '60'
    # for info in ISDict.values():
    #     data = info.split('|')
    #     name = data[0].replace('~','|')
    #     start = data[1]
    #     end = data[2]
    #     strand = data[3]
    #     flag = data[4]
    #     length1 = abs(int(end)-int(start)) + 1
    #     if flag == "IS":
    #         color = "green"
    #         archiveLen = (length1/length)*2*math.pi*radius*alpha
    #         if archiveLen <= 10:
    #             arrowendlength = str(archiveLen)
    #             arrowendwidth = "5"
    #             arrowstartlength = str(archiveLen)
    #             arrowstartwidth = "5"
    #         else:
    #             arrowendlength = "10"
    #             arrowendwidth = "5"
    #             arrowstartlength = "10"
    #             arrowstartwidth = "5"
    #         if strand == "1":
    #             content = '<!-- {name} -->\n<trackmarker start="{start}" end="{end}" markerstyle="stroke:#000;fill:{color};" arrowendlength="{arrowendlength}" arrowendwidth="{arrowendwidth}" wadjust="10" vadjust="-5">\n<markerlabel labelstyle="font-size:10px" text="{name}" vadjust="{vadjust}" hadjust="{hadjust}" valign="{valign}" showline="1" linevadjust="{linevadjust}"  lineclass="labelline"></markerlabel>\n</trackmarker>\n'.format(name=name,start=start,end=end,arrowendlength=arrowendlength,arrowendwidth=arrowendwidth,vadjust=vadjust,hadjust=hadjust,valign=valign,linevadjust=linevadjust,color=color)
    #         elif strand == "-1":
    #             content = '<!-- {name} -->\n<trackmarker start="{start}" end="{end}" markerstyle="stroke:#000;fill:{color};" arrowstartlength="{arrowstartlength}" arrowstartwidth="{arrowstartwidth}" wadjust="10" vadjust="-5">\n<markerlabel labelstyle="font-size:10px" text="{name}" vadjust="{vadjust}" hadjust="{hadjust}" valign="{valign}" showline="1" linevadjust="{linevadjust}"  lineclass="labelline"></markerlabel>\n</trackmarker>\n'.format(name=name,start=start,end=end,arrowstartlength=arrowstartlength,arrowstartwidth=arrowstartwidth,vadjust=vadjust,hadjust=hadjust,valign=valign,linevadjust=linevadjust,color=color)
    #     writeHtml(outfile,content)
    #     plasmid = plasmid + content
    contentEnd = '</plasmidtrack>\n</plasmid>\n</div>\n</body>\n</html>\n'
    writeHtml(outfile,contentEnd)
    plasmid = plasmid + '</plasmidtrack>\n</plasmid>\n</div>\n'
    # return plasmid
if __name__ == "__main__":
    angularPlasmid(sys.argv[1])