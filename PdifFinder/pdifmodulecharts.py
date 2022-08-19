import os
import random
import shutil
import string
def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

# 得到基因方向及对应颜色
def direction(indir):
    geneli = indir + '/genelist.txt'
    orientation = {}
    with open(geneli, 'r') as g:
        geneli = g.readlines()
        for line1 in geneli:
            genelist = line1.strip().split('\t')
            seqID = genelist[0]
            genename = genelist[1]
            start = genelist[2]
            ori = genelist[-3]
            orientation[seqID+genename+start] = ori
    return orientation




def pdifmoduleechart(indir, scriptsDir):
    pdifmoduleli = indir + '/pdifmodule_list1.txt'
    outfile = indir + '/pdifmodule.svg'
    genecolor = scriptsDir + '/data/genecolor.txt'
    geneli = indir + '/genelist.txt'
    colorset = {}
    with open(genecolor, 'r') as g:
        genecolorlist = g.readlines()
        for linec in genecolorlist:
            genename1 = linec.strip().split('\t')[0]
            gc = linec.strip().split('\t')[1]
            colorset[genename1] = gc
    content1 = """<defs>
                <style type='text/css'><![CDATA[
                    .svglite line, .svglite polyline, .svglite polygon, .svglite path, .svglite rect, .svglite circle {
                      fill: none;
                      stroke: #000000;
                      stroke-linecap: round;
                      stroke-linejoin: round;
                      stroke-miterlimit: 10.00;
                    }
                    .svglite text {
                      white-space: pre;
                    }
                ]]></style>\n</defs>\n"""
    with open(pdifmoduleli) as f:
        pdifmodulelist = f.readlines()
        Len = len(pdifmodulelist)
        width = 60 + 60 * (Len - 1)
        content0 = f"""<?xml version='1.0' encoding='UTF-8' ?>\n<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' class='svglite' width='504' height='{width}' viewBox='0 0 504 {width}'>\n"""
        content2 = f"""<rect width='100%' height='100%' style='stroke: none; fill: #FFFFFF;'/>\n<defs>\n<clipPath id='cpMC4wMHw1MDQuMDB8MC4wMHw1MDQuMDA='>\n<rect x='0.00' y='0.00' width='{width}' height='{width}' />\n</clipPath>\n</defs>\n<g clip-path='url(#cpMC4wMHw1MDQuMDB8MC4wMHw1MDQuMDA=)'>\n<rect x='0.00' y='-0.000000000000057' width='{width}' height='{width}' style='stroke-width: 1.07; stroke: #FFFFFF; fill: #FFFFFF;' />\n</g>\n"""
        content = content0 + content1 + content2
        orientation = direction(indir)
        contenth = ""
        for n in range(0, len(pdifmodulelist)):
            line = pdifmodulelist[n]
            list = line.strip().split('\t')
            ID = list[0][1:]
            pdifmodule = list[1].split('+')
            name = "-".join(pdifmodule)
            sites = list[2]
            sitee = list[-1]
            modulest = int(list[3])
            moduleend = int(list[-2])
            modulelen = moduleend - modulest + 1
            str1 = "".join(random.sample(string.ascii_letters + string.digits, 35))
            height = 27.31 + 60 * n
            y = height + 13
            recy = 5.48 + 62 * n
            content3 = f"""<defs>\n<clipPath id='{str1}'>\n<rect x='50.23' y='{recy}' width='400.00' height='43.66' />\n</clipPath>\n</defs>\n<g clip-path='url(#{str1})'>\n<polyline points='66.23,{height} 429.23,{height} ' style='stroke-width: 2.13; stroke: #BEBEBE; stroke-linecap: butt;' />\n<text x='66.23' y='{height}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{sites}</text>\n<text x='429' y='{height}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{sitee}</text>\n<text x='66.23' y='{height + 10}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{modulest}</text>\n<text x='429' y='{height + 10}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{moduleend}</text>\n<text x='247.5' y='{y}' text-anchor='middle' style='font-size: 6.00px; font-weight:bold; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{name}</text>\n"""
            contentm = ""
            for i in range(1, len(pdifmodule) - 1):
                gene = pdifmodule[i]
                start = int(list[i + 4].strip().split('-')[0])
                end = int(list[i + 4].strip().split('-')[1])
                genelen = end - start + 1
                color = colorset[gene]
                ori = orientation[ID + gene + str(start)]
                str2 = "".join(random.sample(string.ascii_letters + string.digits, 35))
                if ori == "minus":
                    H = height - 4.25
                    h = height + 4.25
                    ws = ((start - modulest) / modulelen) * 363 + 61.23
                    we = (end - modulest) / modulelen * 363 + 61.23
                    lpx = (we - ws) / 2 + ws
                    lpy = height - 7
                    hl = h + 1.42
                    Hh = H - 1.42
                    al = ws + ((we - ws) / 3)
                    content4 = f"""<polygon points='{we},{H} {we},{h} {al},{h} {al},{hl} {ws},{height} {al},{Hh} {al},{H} ' style='stroke-width: 0.64; fill: #{color};' />\n<text x='{lpx}' y='{lpy}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{gene}</text>\n"""
                elif ori == "plus":
                    H = height - 4.25
                    h = height + 4.25
                    ws = ((start - modulest) / modulelen) * 363 + 61.23
                    we = ((end - modulest) / modulelen) * 363 + 61.23
                    lpx = (we - ws) / 2 + ws
                    lpy = height - 7
                    hl = h + 1.42
                    Hh = H - 1.42
                    al = we - ((we - ws) / 3)
                    content4 = f"""<polygon points='{ws},{H} {ws},{h} {al},{h} {al},{hl} {we},{height} {al},{Hh} {al},{H} ' style='stroke-width: 0.64; fill: #{color};' />\n<text x='{lpx}' y='{lpy}' text-anchor='middle' style='font-size: 6.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{gene}</text>\n"""
                contentm = contentm + content4
            ending = f"""</g>\n<text x='35' y='{height}' text-anchor='middle' style='font-size: 6.00px; font-weight:bold; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>{ID}</text>\n<g clip-path = 'url(#cpMC4wMHw1MDQuMDB8MC4wMHw1MDQuMDA=)' >\n</g>"""
            contenth = contenth + content3 + contentm + ending
        # a = f"""\n<text transform='translate(13.50,{width/2}) rotate(-90)' text-anchor='middle' style='font-size: 11.00px; font-family: "DejaVu Sans";' lengthAdjust='spacingAndGlyphs'>pdifmodule</text>\n</svg>"""
        a = f"""\n</svg>"""
        content = content + contenth + a
        with open(outfile, 'w') as o:
            o.write(content)
    os.remove(pdifmoduleli)
    os.remove(geneli)
