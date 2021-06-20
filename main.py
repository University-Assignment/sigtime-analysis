from flask import Flask, request
from flask_pymongo import PyMongo
from bson.objectid import ObjectId
import smtplib
import ssl
import os
import shutil
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerExtractor import sigpro as sig
from Bio import SeqIO
import pandas as pd
import numpy as np
import scipy.stats as stat
from pdf2image import convert_from_path

def sig_extract(id):
    sig.sigProfilerExtractor("matrix", "C:/Users/ghks0/Desktop/Project/sigtime-server/uploads/%s/results" % (id), "C:/Users/ghks0/Desktop/Project/sigtime-server/uploads/%s/output/SBS/test.SBS96.all" % (id), reference_genome="GRCh38",
                             minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=-1)


if __name__ == "__main__":
    port = 465
    context = ssl.create_default_context()

    app = Flask(__name__)
    app.config["MONGO_URI"] = "mongodb://localhost:27017/sigtime"
    mongo = PyMongo(app)

    db_results = mongo.db.results
    db_files = mongo.db.files

    # genInstall.install('GRCh38', rsync=False, bash=False)

    @app.route('/v1/analysis', methods=['POST'])
    def analysis():
        req_data = request.get_json(force=True)
        id = req_data['id']

        result_db = db_results.find_one({"_id": ObjectId(id)})
        files = result_db["files"]
        vcf_list = []
        tsv_list = []
        for file_id in files:
            try:
                file = db_files.find_one({"_id": ObjectId(file_id)})
                ext = file["extension"]
                if ext == ".vcf":
                    vcf_list.append(file)
                elif ext == ".tsv":
                    tsv_list.append(file)
            except:
                pass

        prefix = "C:/Users/ghks0/Desktop/Project/sigtime-server"
        uploads_folder = "%s/uploads/%s" % (prefix, id)
        input_folder = "%s/input" % uploads_folder
        deseq_folder = "%s/DESeq" % uploads_folder

        if not os.path.isdir(uploads_folder):
            os.makedirs(uploads_folder)
        if not os.path.isdir(input_folder):
            os.makedirs(input_folder)
        if not os.path.isdir(deseq_folder):
            os.makedirs(deseq_folder)

        for vcf in vcf_list:
            shutil.copy2("%s/%s" % (prefix, vcf["url"]), "%s/%s%s" % (input_folder, vcf["originalname"], vcf["extension"]))

        for tsv in tsv_list:
            shutil.copy2("%s/%s" % (prefix, tsv["url"]), "%s/%s%s" % (deseq_folder, tsv["originalname"], tsv["extension"]))

        matGen.SigProfilerMatrixGeneratorFunc("test", "GRCh38", "%s" % uploads_folder, plot=False,
                                                         exome=False, bed_file=None, chrom_based=False, tsb_stat=False,
                                                         seqInfo=False, cushion=100)
        sig_extract(id)

        best_signature = 3

        db_result = db_results.update_one({"_id": ObjectId(id)}, {"$set": {"signature": best_signature}})
        if db_result is None:
            print("insert best_signature error")

        poppler_path = "C:/Users/ghks0/Desktop/Project/sigtime-analysis/poppler-21.03.0/Library/bin"

        pages = convert_from_path(
            '%s/results/SBS96/SBS96_selection_plot.pdf' % uploads_folder,
            poppler_path=poppler_path) # poppler_path는 윈도우에서만 필요

        for page in pages:
            page.save(
                '%s/results/SBS96/SBS96_selection_plot.jpg' % uploads_folder,
                'JPEG')

        pages = convert_from_path(
            '%s/results/SBS96/All_Solutions/SBS96_%d_Signatures/Activities/SBS96_S%d_NMF_Activity_Plots.pdf' % (uploads_folder, best_signature, best_signature),
            poppler_path=poppler_path) # poppler_path는 윈도우에서만 필요

        for page in pages:
            page.save(
                '%s/results/SBS96/All_Solutions/SBS96_%d_Signatures/Activities/SBS96_S%d_NMF_Activity_Plots.jpg' % (uploads_folder, best_signature, best_signature),
                'JPEG')

        pages = convert_from_path(
            '%s/results/SBS96/All_Solutions/SBS96_%d_Signatures/Signatures/Signature_plotSBS_96_plots_S%d.pdf' % (uploads_folder, best_signature, best_signature),
            poppler_path=poppler_path) # poppler_path는 윈도우에서만 필요

        page_i = 0
        for page in pages:
            page.save(
                '%s/results/SBS96/All_Solutions/SBS96_%d_Signatures/Signatures/Signature_plotSBS_96_plots_S%d_%d.jpg' % (uploads_folder, best_signature, best_signature, page_i),
                'JPEG')
            page_i += 1

        result_folder = "%s/results/SBS96" % uploads_folder
        output_folder = uploads_folder  # <= 최종 결과를 여기다 저장

        def changeValue(before, after):
            if (before == 'G' and after == 'C') or (before == 'A' and after == 'T'):
                return after, before
            if before == 'G' and after == 'A':
                return 'C', 'T'
            if before == 'G' and after == 'T':
                return 'C', 'A'
            if before == 'A' and after == 'G':
                return 'T', 'C'
            if before == 'A' and after == 'C':
                return 'T', 'G'
            return before, after

        genes = 'gene.txt'
        seqs = list(SeqIO.parse('hg38.fa', 'fasta'))
        chrs = [s.id for s in seqs]

        flist = []

        for fname in os.listdir(input_folder):
            if fname[-4:] == '.vcf':
                flist.append(fname)

        g = pd.read_csv(genes, sep='\t', header=None)
        t1 = pd.read_csv('snvtype01.txt', header=None)
        t2 = pd.read_csv('snvtype02.txt', header=None)
        types = t1.loc[:, 0] + t2.loc[:, 0].str.split('>').str.get(1)
        columns = g.loc[:, 0]

        yMatList = []
        geneList = []

        with open(genes, 'r') as f1:
            for line1 in f1.readlines():
                tmp1 = line1.split()
                gene = {
                    "startPos": int(tmp1[1]),
                    "endPos": int(tmp1[2]),
                    "chr": tmp1[3]
                }
                geneList.append(gene)

        for i in range(len(flist)): # vcf 파일 반복
            df = pd.DataFrame(0, index=types, columns=columns)
            fname = flist[i]
            with open(input_folder + '/' + fname, 'r') as f2:
                for line2 in f2.readlines():
                    if line2[0] == '#': continue
                    tmp2 = line2.split()
                    if tmp2[0] not in chrs:
                        continue

                    Pos = int(tmp2[1])
                    loc = int(tmp2[1])
                    subseq = str(seqs[chrs.index(tmp2[0])][loc - 2:loc + 1].seq)  # loc이 index인지 순서인지
                    before, after = changeValue(subseq[1], tmp2[4])
                    subseq2 = subseq[0] + before + subseq[2]
                    if subseq2 + tmp2[4] not in df.index: continue

                    j = 0
                    for gene in geneList:
                        colname = columns[j]
                        j += 1
                        startPos = gene["startPos"]
                        endPos = gene["endPos"]
                        Chr = gene["chr"]
                        if tmp2[0] != Chr or len(tmp2[3]) > 1 or len(tmp2[4]) > 1 or Pos < startPos or Pos > endPos:
                            continue
                        if subseq[1] == tmp2[3]:
                            df.loc[subseq2 + after, colname] += 1
                            break
                        else:
                            print('Not match! Fname: {}, subseq: {}, vcf: {} to {}'.format(fname, subseq, tmp2[3],
                                                                                           tmp2[4]))
            print(fname)
            yMatList.append(df)
            df.to_csv("C:/Users/ghks0/Desktop/" + fname + ".csv")

        # dirname = "yMat"
        # for fname in os.listdir(dirname):
        #     if fname[-4:] == '.csv':
        #         yMatList.append(pd.read_csv(dirname + "/" + fname, index_col=0))

        P = pd.read_csv(result_folder + "/All_Solutions/SBS96_%d_Signatures/Signatures/SBS96_S%d_Signatures.txt" % (best_signature, best_signature), sep='\t',
                        index_col=0)  # Extractor 결과 중 process 사용
        E = pd.read_csv(result_folder + "/All_Solutions/SBS96_%d_Signatures/Activities/SBS96_S%d_NMF_Activities.txt" % (best_signature, best_signature), sep='\t',
                        index_col=0)  # Extractor 결과 중 exposure 사용

        samples = E.index
        sig = P.columns
        types = P.index

        # total -> 분자 합

        contri_list = []
        contri_name_list = []

        for s in range(len(samples)):
            df = pd.DataFrame()
            for j in range(len(types)):
                numerator = []
                sum = 0.0
                for i in range(len(sig)):
                    p = P.iloc[j, i]
                    e = E.iloc[s, i]
                    result = p * e
                    numerator.append(result)
                    sum += result
                for ss in range(len(sig)):
                    if sum == 0.0:
                        df.loc[j, ss] = 0
                    else:
                        df.loc[j, ss] = numerator[ss] / sum
                # print(samples[s])
            # df.to_csv('output/' + samples[s] + '_C.csv') # contribution 결과 파일
            contri_list.append(df)
            contri_name_list.append(samples[s][4:11])

        cList = []
        DESeq = []
        column = ['gene', 'sig']

        for i in range(len(contri_name_list)):
            cList.append(contri_list[i])
            column.append(contri_name_list[i])

        column.extend(column[2:7])
        column.extend(['c', 'p'])

        for fName in os.listdir(deseq_folder):
            if fName[-4:] == '.tsv':
                DESeq.append(pd.read_csv(deseq_folder + "/" + fName, sep='\t')) # 입력받는 파일

        df = pd.DataFrame(columns=column)

        index = 0
        genes = yMatList[0].columns
        sig = cList[0].columns
        types = yMatList[0].index

        count = 0

        for gene in genes:  # 정수인지, 합 확인
            for s in sig:
                df.loc[index, 0:2] = [gene, str(s)]
                v1 = []
                for i in range(len(yMatList)):
                    sum = 0
                    type_i = 0
                    for type in types:
                        sum += yMatList[i].loc[type, gene] * cList[i].iloc[type_i, int(s)]
                        type_i += 1
                    v1.append(sum)
                df.loc[index, 2:7] = tuple(map(str, v1))
                log_ = []
                for i in range(len(DESeq)):
                    log_.append(DESeq[i].loc[gene, 'log2FoldChange'])
                df.loc[index, 7:12] = tuple(map(str, log_))
                r = np.nan
                p_value = np.nan
                if np.isfinite(log_[0]):
                    r, p_value = stat.pearsonr(v1, log_)
                df.loc[index, 12:14] = [r, p_value]
                index += 1
            count += 1
            if count % 1000 == 0:
                print(count, "/ 58440", str(count / 58440 * 100) + "%") # 진행률 나타냄

        df.to_csv(output_folder + '/NewMat.csv')  # 이게 제공하는 전체 파일

        # df = pd.read_csv(output_folder + '/NewMat.csv', index_col=0)

        AllResult = df

        PartRe_list = []

        for s in range(best_signature):  # 시그니쳐 개수에 따라 range 변경
            df = pd.DataFrame(columns=AllResult.columns)  # columns=columns
            idx = 0
            for i in range(len(AllResult)):  # 총 행만큼 반복하기
                if AllResult.iloc[i, 12] == '' or AllResult.iloc[i, 13] == '': continue
                # Num = int(AllResult.iloc[i,0])
                Sig = int(AllResult.iloc[i, 1])
                Cor = float(AllResult.iloc[i, 12])
                P_val = float(AllResult.iloc[i, 13])
                if Sig == s and Cor <= -0.9 and P_val <= 0.05:
                    # print(Sig, Cor, P_val) # 확인을 위해서 출력
                    df.loc[idx, :] = AllResult.iloc[i, :]
                    idx += 1
            df.drop(df.columns[1], axis=1)
            df.to_csv(output_folder + '/NewMat_Sig' + str(s) + '.csv', sep=',', index=True, header=True)  # 이거 표로 출력
            PartRe_list.append(df)
        #     print('Done Sig' + str(s))

        for i in range(len(PartRe_list)):
            columns = []

            df = PartRe_list[i].iloc[:, 2:12]
            df = df.T
            df = df.astype(float).mean(axis=1)

            data1 = df[:5]
            data2 = df[5:]
            columns.append(PartRe_list[i].iloc[0, 2:7])

            csv = pd.DataFrame()

            data2.index = PartRe_list[i].columns[2:7]

            # csv = csv.append(columns)
            csv = csv.append(data1, ignore_index=True)
            csv = csv.append(data2, ignore_index=True)
            # csv.iloc[-1] = data2

            csv = csv.T

            csv.columns = ["mutation", "DESeq"]

            csv.to_csv(output_folder + '/FinalMat_Sig' + str(i) + '_data.csv')

        receiver_email = result_db["email"]
        db_result = db_results.update_one({"_id": ObjectId(id)}, {"$set": {"complete": True}})
        if db_result is not None:
            with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
                to = [receiver_email]
                subject = 'test url'
                email_text = """\
From: %s
To: %s
Subject: %s
    
http://localhost/result/%s""" % ("보낸 사람", ", ".join(to), subject, id)

                server.login("아이디", "비번")
                server.sendmail("보낸 사람", receiver_email, email_text)
        return 'complete'

    app.run(port=5000)