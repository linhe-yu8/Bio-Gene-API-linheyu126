from fastapi import FastAPI
import pymysql
import requests
import csv
from pydantic import BaseModel

app = FastAPI(title="基因序列分析")
class GeneQuery(BaseModel):
    gene_id:str

class CompareQuery(BaseModel):
    normal_id:str
    patient_id:str

codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
amino_acid_names = {
    'A': '丙氨酸 (Ala)', 'R': '精氨酸 (Arg)', 'N': '天冬酰胺 (Asn)', 'D': '天冬氨酸 (Asp)',
    'C': '半胱氨酸 (Cys)', 'Q': '谷氨酰胺 (Gln)', 'E': '谷氨酸 (Glu)', 'G': '甘氨酸 (Gly)',
    'H': '组氨酸 (His)', 'I': '异亮氨酸 (Ile)', 'L': '亮氨酸 (Leu)', 'K': '赖氨酸 (Lys)',
    'M': '甲硫氨酸 (Met)', 'F': '苯丙氨酸 (Phe)', 'P': '脯氨酸 (Pro)', 'S': '丝氨酸 (Ser)',
    'T': '苏氨酸 (Thr)', 'W': '色氨酸 (Trp)', 'Y': '酪氨酸 (Tyr)', 'V': '缬氨酸 (Val)'
}
amino_acid_weights = {
    'G': 75.07, 'A': 89.09, 'V': 117.15, 'L': 131.17, 'I': 131.17,
    'P': 115.13, 'M': 149.21, 'F': 165.19, 'W': 204.23, 'S': 105.09,
    'T': 119.12, 'C': 121.16, 'N': 132.12, 'Q': 146.15, 'Y': 181.19,
    'D': 133.10, 'E': 147.13, 'K': 146.19, 'R': 174.20, 'H': 155.16
}
hydrophobicity_scores = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
corresponding = {'A':'T','G':'C','T':'A','C':'G'}

#翻译逻辑
def translate(n):
    protein = ""
    n = n.upper()
    start_index = n.find("ATG")
    if start_index == -1:
        return "找不到起始密码子，无法翻译"
    coding = n[start_index:]
    for i in range(0,len(coding)-2,3):
        codon = coding[i:i+3]
        median = codon_table.get(codon,"?")
        if median == "*":
            break
        protein += median
    return protein
#对比正常与患病的异常点位
def get_seq(tid):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": "nucleotide", "id": tid, "rettype": "fasta", "retmode": "text"}
    res = requests.get(url, params=params)
    if res.status_code == 200:
        return "".join(res.text.strip().split('\n')[1:]).upper()
    return ""
#主要功能实现
def Execution_Tool(Gen_id,L,GC,MOW,HY_num,HY_IN):
    try:
        db=pymysql.connect(host='localhost',user='root',password='密码',database='fastapi_db',charset='utf8mb4')
        cursor=db.cursor()
        sql = "insert into history(gene_id,length,gc,end_MOW,num_Hydrophilicity,intutive_Hydrophilicity) values (%s,%s,%s,%s,%s,%s)"
        cursor.execute(sql,(Gen_id,L,GC,MOW,HY_num,HY_IN))
        db.commit()
        db.close()
        print(f"{Gen_id}的相关数据已经存入MySQL!")
    except Exception as f:
        print("数据库连接失败，原因：",f)
#各部分具体功能
def main_functions(gene_id:str):
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db":"nucleotide","id":gene_id,"rettype":"fasta","retmode":"text"}
    response = requests.get(url,params=params)
    if response.status_code == 200:
        while_text = response.text.strip().split('\n')
        hea = while_text[0]
        content = "".join(while_text[1:])
    #算G C含量
    L_content = len(content)
    g_content = content.count("G")
    c_content = content.count("C")
    if L_content > 0:
        co_GC = (g_content+c_content)/L_content
        CO_GC = f"GC含量：{co_GC:.2%}"
    print(CO_GC)

    #翻译氨基酸
    translate_reslut = translate(content)
    L_total = len(translate_reslut)
    print(translate_reslut)

    #反向互补链
    print(content)
    collect = []
    for k in content:
        Complementary = corresponding.get(k,"N")
        collect.append(Complementary)
    Reverse = "".join(collect)[::-1]
    print(Reverse)
    Initial = "".join([corresponding.get(lll,"N") for lll in Reverse[::-1]])
    if Initial == content:
        print("反向互补链绝对正确！")
    else:
        print("互补配对出错了！")
    #氨基酸组成成分分析

    super_dict = {}
    for Amino_acid in translate_reslut:
        if Amino_acid not in super_dict:
            Cn = amino_acid_names.get(Amino_acid,"未知")
            super_dict[Amino_acid] = {"中文名":Cn,"数量":0,"占比":"0%"}
        super_dict[Amino_acid]["数量"] += 1
    for Amino_acid,info in super_dict.items():
        percentage = (info["数量"] / L_total) * 100
        info["占比"] = f"{percentage:.2f}%"

    print(super_dict)
    #算这串氨基酸的分子量以及亲水性
    #分子量：

    Mow_dict = {}
    total_MOW = 0.0
    for eve_amino in translate_reslut:
        js = amino_acid_weights.get(eve_amino,0)
        total_MOW += js
    if len(translate_reslut) > 1:
        loss_water = (len(translate_reslut)-1) * 18.015
        Molecular_weight = total_MOW - loss_water
    print(f"真实分子量：{Molecular_weight:.2f}")
    #亲水性：

    ccnt = 0.0
    for anmino in translate_reslut:
        ccnt +=  hydrophobicity_scores.get(anmino,0)
    avg = ccnt / len(translate_reslut)
    hydrophobicity_final = f"亲水性：{avg:.3f}"
    print(hydrophobicity_final)
    result_Hydrophilicity = ""
    if avg>0:
        print("疏水")
        result_Hydrophilicity = "疏水"
    else:

        print("亲水")
        result_Hydrophilicity = "亲水"
    Execution_Tool(gene_id,L_content,CO_GC,f"{Molecular_weight:.2f}",f"{avg:.3f}",result_Hydrophilicity)
    return {
        "基因ID": gene_id,
        "序列长度": L_content,
        "GC含量": CO_GC,
        "蛋白质总长度": L_total,
        "蛋白质序列(前50)": translate_reslut[:50] + "...",
        "反向互补链(前50)": Reverse[:50] + "...",
        "分子量(Da)": f"{Molecular_weight:.2f}",
        "亲水性得分(GRAVY)": f"{avg:.3f}",
        "理化性质": result_Hydrophilicity,
        "氨基酸成分分析": super_dict
    }
#接受上述return结果
@app.post("/analyze")
def analyze_gene(query:GeneQuery):
    return main_functions(query.gene_id)
#呈现数据库历史记录
@app.get("/history")
def get_history():
    try:
        db=pymysql.connect(host='localhost',user='root',password='密码',database='fastapi_db',charset='utf8mb4')
        cursor=db.cursor(pymysql.cursors.DictCursor)
        cursor.execute("select * from history order by id desc")
        results = cursor.fetchall()
        db.close()
        return {"status": "success", "data": results}
    except Exception as f:
        return {"status": "error", "message": str(f)}

@app.post("/compare")
def compare_mutations(query:CompareQuery):
    normal_seq = get_seq(query.normal_id)
    patient_seq = get_seq(query.patient_id)
    if not normal_seq or not patient_seq:
        return {"status": "error", "message": "ID或者网络有问题"}
    mutations = []
    min_L = min(len(normal_seq),len(patient_seq))

    for i in range(min_L):
        if normal_seq[i] != patient_seq[i]:
            mutations.append(f"第{i+1}位：{normal_seq[i]} -> {patient_seq[i]}")
    mutation_total = len(mutations)
    try:
        db = pymysql.connect(host='localhost', user='root', password='密码', database='fastapi_db',
                             charset='utf8mb4')
        cursor = db.cursor()
        sql = "insert into compare_history (normal_id, patient_id, mutations_count) values (%s, %s, %s)"
        cursor.execute(sql, (query.normal_id, query.patient_id, mutation_total))
        db.commit()
        db.close()
        print(f"已成功存入compare_history表")

    except Exception as f:
        print("入库失败:", f)

    return{
        "正常ID:":query.normal_id,
        "异常ID:":query.patient_id,
        "突变数量总计:":len(mutations),
        "突变详情(前10个):":mutations[:10]
    }

