"""
肛周脓肿 与 血栓栓塞性疾病的孟德尔随机化分析
数据来源: IEU OpenGWAS 数据库 (https://gwas.mrcieu.ac.uk/)
使用 TwoSampleMR 包进行分析
"""

import os
import json
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy import stats

# 尝试设置中文字体
try:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
except:
    pass

from TwoSampleMR import MRData,harmonise,mr,mr_causal_direction,synthesis,forest,leave_one_out
from ieugwaspy import GWAS

# ============================================================
# 第一步：查找暴露和结局的GWAS数据
# ============================================================
print("=" * 70)
print("第一步：搜索IEU OpenGWAS数据库中的GWAS数据")
print("=" * 70)

# 使用ieugwaspy搜索GWAS数据
def search_gwas(keyword, poolsize=1000):
    """搜索GWAS数据"""
    try:
        res = GWAS().gdainfo(keyword=keyword, poolsize=poolsize)
        return res
    except Exception as e:
        print(f"  搜索 '{keyword}' 失败: {e}")
        return None

# 搜索暴露：肛周脓肿 (Perianal abscess)
print("\n[暴露] 搜索：Perianal abscess / Anal abscess / Anorectal abscess")
exp_search = search_gwas("perianal abscess")
if exp_search is not None and len(exp_search) > 0:
    print(f"  找到 {len(exp_search)} 个相关GWAS研究")
    print(exp_search[["id","trait","year","nsnp","sample_size","author"]].head(10).to_string())
else:
    # 备选搜索
    exp_search2 = search_gwas("anal abscess")
    if exp_search2 is not None and len(exp_search2) > 0:
        print(f"  备选搜索找到 {len(exp_search2)} 个相关GWAS研究")
        print(exp_search2[["id","trait","year","nsnp","sample_size","author"]].head(10).to_string())
        exp_search = exp_search2
    else:
        exp_search3 = search_gwas("abscess")
        if exp_search3 is not None and len(exp_search3) > 0:
            print(f"  扩展搜索 'abscess' 找到 {len(exp_search3)} 个相关GWAS研究")
            print(exp_search3[["id","trait","year","nsnp","sample_size","author"]].head(10).to_string())
            exp_search = exp_search3

# 搜索结局
outcomes = {
    "脑梗": "cerebral infarction",
    "心梗": "myocardial infarction",
    "深静脉血栓": "deep vein thrombosis",
    "动脉栓塞": "arterial embolism"
}

outcome_results = {}
for name_cn, name_en in outcomes.items():
    print(f"\n[结局] 搜索：{name_cn} ({name_en})")
    res = search_gwas(name_en)
    if res is not None and len(res) > 0:
        print(f"  找到 {len(res)} 个相关GWAS研究")
        print(res[["id","trait","year","nsnp","sample_size","author"]].head(5).to_string())
        outcome_results[name_cn] = res
    else:
        print(f"  未找到直接匹配的GWAS研究，尝试备选词...")
        alt_keywords = {
            "脑梗": ["ischemic stroke", "stroke", "brain infarction"],
            "心梗": ["heart attack", "coronary heart disease", "CHD"],
            "深静脉血栓": ["venous thrombosis", "DVT", "venous thromboembolism"],
            "动脉栓塞": ["peripheral arterial disease", "atherosclerosis", "arterial thrombosis"]
        }
        found_alt = False
        for alt_kw in alt_keywords.get(name_cn, []):
            res2 = search_gwas(alt_kw)
            if res2 is not None and len(res2) > 0:
                print(f"  备选 '{alt_kw}' 找到 {len(res2)} 个GWAS研究")
                print(res2[["id","trait","year","nsnp","sample_size","author"]].head(5).to_string())
                outcome_results[name_cn] = res2
                found_alt = True
                break
        if not found_alt:
            outcome_results[name_cn] = None
