import os


def ISfinder(inFile,outdir,blastnPath,scriptsDir,cov=60,ident=90):
    database=scriptsDir+"/ISDatabase/IS"
    os.system(blastnPath+" -num_threads 8 -task blastn -query "+inFile+" -db "+database+" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq slen\" -out "+outdir+"/tmp/ISfinder.xls")
    file1 = outdir+"/tmp/ISfinder.xls"
    f = open(file1)
    lines1 = f.readlines()
    f.close()
    with open(file1,'w') as w:
        w.write('query\tsubject|family|group|organism type\tidentity\tlength\tmismatch\tgapopen\tquery start\tquery end\tsubject start\tsubject end\tevalue\tbitscore\tquery seq\tsubject seq\tsubject length\n')
    for line1 in lines1:
        with open(file1,'a') as w:
            w.write(line1)
    blast_result = outdir+'/tmp/ISfinder.xls'
    filter_result_temp = outdir+'/tmp/ISfinder.filter.temp.xls'
    filter_result = outdir+'/ISfinder.filter.xls'
    filter_blast_result(blast_result,filter_result_temp,cov,ident)
    find_best_match(filter_result_temp,outdir)
    os.remove(filter_result_temp)
    os.remove(blast_result)
    # print("ISfinder done")
    
    
def filter_blast_result(blast_result,filter_result_temp,coverage,identity):
    if os.path.exists(filter_result_temp):
        os.remove(filter_result_temp)
    f = open(blast_result)
    i = 0
    for lines in f:
        i = i + 1
        if i == 1:
            with open(filter_result_temp,"a") as w:
                title = '\t'.join(lines.split("\t")[0:14])
                w.write(title+'\tslen\n')
            continue
        lines = lines.strip()
        if lines:
            line = lines.split("\t")
            ident = float(line[2])
            slen = float(line[14])
            match_length = abs(int(line[9])-int(line[8])) + 1
            cov = round(match_length*100/slen,2)
            if ident >= identity and  cov >= coverage and match_length>= 50:
                newline = '\t'.join(line[0:2]) + '\t' + str(round(ident,2)) + '\t' + '\t'.join(line[3:15])
                with open(filter_result_temp,"a") as w:
                    w.write(newline+'\n')
    f.close()
    
def ret_start_end(start1,end1,start2,end2):
    if int(start1) <= int(start2):
        start = start1
    else:
        start = start2
    if int(end1) >= int(end2):
        end = end1
    else:
        end = end2
    return start,end
    
def check_interaction(ini_pos,prepos):
    rstart = int(ini_pos.split(":")[0])
    rend = int(ini_pos.split(":")[1])
    for po2 in prepos:
        start = int(po2.split(":")[0])
        end = int(po2.split(":")[1])
        if start <= rend and end >= rstart:
            flag = 1
            min_start,max_end = ret_start_end(start,end,rstart,rend)
            position = str(min_start) + ':' + str(max_end)
            return flag,position
            break
        else:
            if po2 == prepos[-1]:
                flag = 0
                return flag,ini_pos
            else:
                continue
    
def find_best_match(filter_result_temp,resultdir):
    outfile = resultdir+"/ISfinder.filter.xls"
    f = open(filter_result_temp)
    match_dict = {}
    i = 0
    for lines in f:
        i = i + 1
        if i == 1:
            with open(outfile,"a") as w:
                w.write(lines)
            continue
        lines = lines.strip()
        line = lines.split("\t")
        query = line[0]
        start = line[6]
        end = line[7]

        pos = start + ":" + end
        if query not in match_dict:
            match_dict[query] = {}
            if pos not in match_dict[query]:
                pre_pos = []
                for k1 in match_dict[query].keys():
                    pre_pos.append(k1)
                if pre_pos != []:
                    ini_pos = pos
                    flag,ret_pos = check_interaction(ini_pos,pre_pos)
                    if flag == 1:
                        if ret_pos in match_dict[query].keys():
                            match_dict[query][ret_pos].append(lines)
                    else:
                        match_dict[query][ret_pos] = []
                        match_dict[query][ret_pos].append(lines)
                else:
                    match_dict[query][pos] = []
                    match_dict[query][pos].append(lines)
            else:
                match_dict[query][pos].append(lines)
        else:
            if pos not in match_dict[query]:
                pre_pos = []
                for k1 in match_dict[query].keys():
                    pre_pos.append(k1)
                if pre_pos != []:
                    ini_pos = pos
                    flag,ret_pos  = check_interaction(ini_pos,pre_pos)
                    if flag == 1:
                        if ret_pos in match_dict[query].keys():
                            match_dict[query][ret_pos].append(lines)
                    else:
                        match_dict[query][ret_pos] = []
                        match_dict[query][ret_pos].append(lines)
            else:
                match_dict[query][pos].append(lines)
    f.close()
    for k,v in match_dict.items():
        for k1,v1 in v.items():
            max_ident = 0
            max_cov = 0
            j = -1
            for vs in v1:
                j = j + 1
                v_list = vs.split("\t")
                query = v_list[0]
                pos = v_list[6] + ":" + v_list[7]
                ident = float(v_list[2])
                cov = len(v_list[12])*100/float(v_list[-1])
                if ident >= max_ident and cov >= max_cov:
                    max_ident = ident
                    max_cov = cov
                    max_po = k1
                    max_query = query
                    max_pos = j
            with open(outfile,"a") as w:
                w.write(match_dict[max_query][max_po][max_pos]+'\n')
    
