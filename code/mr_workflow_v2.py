"""
肛周脓肿 与 血栓栓塞性疾病的孟德尔随机化分析
数据来源: IEU OpenGWAS 数据库 (https://gwas.mrcieu.ac.uk/)
通过直接HTTP API获取真实GWAS数据
"""

import os
import sys
import json
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import requests

# 设置matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Microsoft YaHei, SimHei, Arial Unicode MS"
plt.rcParams["axes.unicode_minus"] = False

from scipy import stats

OUT_DIR = r"D:\肛周脓肿与血栓性疾病"
os.makedirs(OUT_DIR, exist_ok=True)

IEU_API = "https://gwas.mrcieu.ac.uk"

def api_get(endpoint, params=None):
    """发送GET请求到IEU OpenGWAS API"""
    url = f"{IEU_API}{endpoint}"
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r
    except Exception as e:
        print(f"  API请求失败 [{endpoint}]: {e}")
        return None

def search_gwas_api(keyword, max_results=20):
    """通过API搜索GWAS研究"""
    print(f"  搜索关键词: '{keyword}'")
    r = api_get("/api/v2/gwasinfo", {"keyword": keyword, "max_results": max_results})
    if r is None:
        return []
    try:
        data = r.json()
        if isinstance(data, list):
            return data
        elif isinstance(data, dict) and "result" in data:
            return data["result"]
        return []
    except Exception as e:
        print(f"  解析搜索结果失败: {e}")
        return []

def get_gwas_variants(gwas_id):
    """获取指定GWAS的所有变异位点"""
    print(f"  获取GWAS变异: {gwas_id}")
    r = api_get(f"/api/v2/variant", {"id": gwas_id})
    if r is None:
        return pd.DataFrame()
    try:
        data = r.json()
        if isinstance(data, list):
            df = pd.DataFrame(data)
        elif isinstance(data, dict):
            if "result" in data:
                df = pd.DataFrame(data["result"])
            else:
                df = pd.DataFrame([data])
        else:
            df = pd.DataFrame()
        print(f"    获得 {len(df)} 个变异")
        return df
    except Exception as e:
        print(f"  解析变异数据失败: {e}")
        return pd.DataFrame()

def get_associations(gwas_id, pval=5e-8):
    """获取GWAS中达到显著性的SNP位点"""
    print(f"  获取 {gwas_id} 的显著SNP (p<{pval})...")
    r = api_get(f"/api/v2/associations", {"id": gwas_id, "pval": pval, "size": 10000})
    if r is None:
        return pd.DataFrame()
    try:
        data = r.json()
        if isinstance(data, list):
            df = pd.DataFrame(data)
        elif isinstance(data, dict) and "result" in data:
            df = pd.DataFrame(data["result"])
        else:
            df = pd.DataFrame()
        print(f"    获得 {len(df)} 个显著SNP")
        return df
    except Exception as e:
        print(f"  解析失败: {e}")
        return pd.DataFrame()

# ============================================================
# 第一步：搜索暴露和结局的GWAS数据
# ============================================================
print("\n" + "="*70)
print("第一步：搜索 IEU OpenGWAS 数据库")
print("="*70)

# 暴露搜索策略
exp_keywords = [
    "perianal abscess",
    "anal abscess",
    "anorectal abscess",
    "perianal disease",
    "anorectal disease",
    "abscess perianal",
]

exposure_candidates = []
for kw in exp_keywords:
    results = search_gwas_api(kw)
    if results:
        exposure_candidates.extend(results)
        if len(results) >= 3:
            break

# 去重
if exposure_candidates:
    ids_seen = set()
    exposure_candidates = [x for x in exposure_candidates if not (x.get("id") in ids_seen or ids_seen.add(x.get("id")))]

print(f"\n暴露候选GWAS: {len(exposure_candidates)} 个")
if exposure_candidates:
    cols = ["id", "trait", "year", "nsnp", "sample_size", "author"]
    for r in exposure_candidates[:10]:
        row = [f"{r.get(c, 'NA')}" for c in cols]
        print(f"  {row}")

# 选择最佳暴露GWAS
exposure_id = None
if exposure_candidates:
    # 选择SNP数量最多的
    def get_nsnp(x):
        try: return int(x.get("nsnp", 0))
        except: return 0
    exposure_candidates.sort(key=get_nsnp, reverse=True)
    for cand in exposure_candidates:
        eid = cand.get("id")
        if eid:
            exposure_id = eid
            print(f"\n选择的暴露GWAS: {eid} | {cand.get('trait')} | nsnp={cand.get('nsnp')} | n={cand.get('sample_size')}")
            break

# 结局搜索
outcome_configs = {
    "脑梗": ["cerebral infarction", "ischemic stroke", "lacunar stroke"],
    "心梗": ["myocardial infarction", "coronary heart disease", "heart attack"],
    "深静脉血栓": ["deep vein thrombosis", "venous thrombosis", "venous thromboembolism"],
    "动脉栓塞": ["peripheral arterial disease", "arterial embolism", "arterial thrombosis"]
}

outcome_ids = {}
outcome_info = {}
for name_cn, keywords in outcome_configs.items():
    print(f"\n--- {name_cn} ---")
    for kw in keywords:
        results = search_gwas_api(kw)
        if results:
            # 选择nsnp最大的
            def get_nsnp(x):
                try: return int(x.get("nsnp", 0))
                except: return 0
            results.sort(key=get_nsnp, reverse=True)
            cand = results[0]
            oid = cand.get("id")
            if oid:
                outcome_ids[name_cn] = oid
                outcome_info[name_cn] = cand
                print(f"  选择: {oid} | {cand.get('trait')} | nsnp={cand.get('nsnp')} | n={cand.get('sample_size')}")
                break
    if name_cn not in outcome_ids:
        print(f"  未找到合适的GWAS")

print("\n" + "="*70)
print("搜索结果汇总")
print("="*70)
print(f"暴露GWAS ID: {exposure_id}")
for k, v in outcome_ids.items():
    print(f"结局 {k}: {v} | {outcome_info[k].get('trait') if k in outcome_info else ''}")

# ============================================================
# 第二步：获取SNP工具变量
# ============================================================
print("\n" + "="*70)
print("第二步：提取SNP工具变量")
print("="*70)

exp_snps_df = pd.DataFrame()
if exposure_id:
    exp_snps_df = get_associations(exposure_id, pval=5e-8)
    if exp_snps_df.empty:
        # 尝试更宽松的p值
        for pv in [1e-5, 5e-6]:
            print(f"  p<5e-8无结果，尝试 p<{pv}...")
            exp_snps_df = get_associations(exposure_id, pval=pv)
            if not exp_snps_df.empty:
                break

out_snps = {}
for name_cn, oid in outcome_ids.items():
    df = get_associations(oid, pval=5e-8)
    if df.empty:
        for pv in [1e-5, 5e-6]:
            df = get_associations(oid, pval=pv)
            if not df.empty:
                print(f"  {name_cn} 使用 p<{pv}: {len(df)} SNP")
                break
    out_snps[name_cn] = df
    print(f"  {name_cn} ({oid}): {len(df)} 显著SNP")

print(f"\n暴露 {exposure_id}: {len(exp_snps_df)} 个显著SNP")
print(f"暴露SNP列名: {list(exp_snps_df.columns)[:15]}")

# ============================================================
# 第三步：孟德尔随机化分析
# ============================================================
print("\n" + "="*70)
print("第三步：执行MR分析")
print("="*70)

results_all = []

def ivw_meta_analysis(b_exp, b_out, se_out):
    """逆方差加权Meta分析 (IVW)"""
    weights = 1.0 / (se_out ** 2)
    b_ivw = np.sum(weights * b_out) / np.sum(weights)
    se_ivw = np.sqrt(1.0 / np.sum(weights))
    z = b_ivw / se_ivw
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return b_ivw, se_ivw, z, p

def simple_mr(b_exp, b_out, se_out):
    """简单MR (比率估计)"""
    ratios = b_out / b_exp
    se_ratios = se_out / np.abs(b_exp)
    # 使用中位数作为稳健估计
    b_med = np.median(ratios)
    # 计算标准误
    mad = np.median(np.abs(ratios - b_med))
    se_med = 1.483 * mad / np.sqrt(len(ratios))
    z = b_med / se_med if se_med > 0 else 0
    p_med = 2 * (1 - stats.norm.cdf(abs(z)))
    return b_med, se_med, z, p_med

def mr_egger(b_exp, b_out, se_out):
    """MR-Egger方法"""
    n = len(b_exp)
    # 权重
    w = 1.0 / (se_out ** 2)
    sum_w = np.sum(w)
    sum_wx = np.sum(w * b_exp)
    sum_wy = np.sum(w * b_out)
    sum_wxx = np.sum(w * b_exp ** 2)
    sum_wxy = np.sum(w * b_exp * b_out)
    
    denom = sum_w * sum_wxx - sum_wx ** 2
    if denom == 0:
        return None, None, None, None
    
    b_0 = (sum_wxy * sum_w - sum_wx * sum_wy) / denom
    b_1 = (sum_wxx * sum_wy - sum_wx * sum_wxy) / denom
    
    # Egger截距的标准误
    q = np.sum(w * (b_out - b_0 - b_1 * b_exp) ** 2)
    se_0 = np.sqrt(sum_wxx / denom)
    se_1 = np.sqrt(sum_w / denom)
    
    z = b_1 / se_1 if se_1 > 0 else 0
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return b_1, se_1, z, p

for name_cn, out_df in out_snps.items():
    print(f"\n>>> 分析: Perianal Abscess -> {name_cn}")
    
    if exp_snps_df.empty or out_df.empty:
        print(f"  跳过: 数据为空")
        results_all.append({
            "exposure": "Perianal abscess",
            "outcome": name_cn,
            "exposure_gwas_id": exposure_id,
            "outcome_gwas_id": outcome_ids.get(name_cn),
            "status": "no_data",
            "methods": {}
        })
        continue
    
    # 找SNP重叠
    # 尝试不同的SNP标识列
    snp_col = None
    for col in ["snp", "rsid", "variant", "id"]:
        if col in exp_snps_df.columns and col in out_df.columns:
            snp_col = col
            break
    
    if snp_col is None:
        print(f"  跳过: 找不到SNP标识列")
        print(f"  暴露列: {list(exp_snps_df.columns)}")
        print(f"  结局列: {list(out_df.columns)}")
        results_all.append({
            "exposure": "Perianal abscess",
            "outcome": name_cn,
            "status": "no_snp_column",
            "methods": {}
        })
        continue
    
    # 合并数据
    merged = pd.merge(exp_snps_df, out_df, on=snp_col, suffixes=("_exp", "_out"))
    print(f"  重叠SNP数量: {len(merged)}")
    
    if len(merged) < 3:
        print(f"  重叠SNP<3, 使用全部SNP进行汇总比较...")
        # 使用尽可能多的数据
        merged = pd.merge(exp_snps_df, out_df, left_on=exp_snps_df.columns[0],
                         right_on=out_df.columns[0], how="inner", suffixes=("_exp", "_out"))
        print(f"  重叠SNP(备选): {len(merged)}")
        if len(merged) < 3:
            results_all.append({
                "exposure": "Perianal abscess",
                "outcome": name_cn,
                "status": "insufficient_snps",
                "methods": {}
            })
            continue
    
    # 识别beta和se列
    beta_exp_col, se_exp_col = None, None
    beta_out_col, se_out_col = None, None
    for col in merged.columns:
        cl = col.lower()
        if "_exp" in col:
            if "beta" in cl or "effect" in cl or "ea_odds" in cl:
                if beta_exp_col is None: beta_exp_col = col
            if cl.endswith("_se") or "standard_error" in cl:
                if se_exp_col is None: se_exp_col = col
        if "_out" in col:
            if "beta" in cl or "effect" in cl or "ea_odds" in cl:
                if beta_out_col is None: beta_out_col = col
            if cl.endswith("_se") or "standard_error" in cl:
                if se_out_col is None: se_out_col = col
    
    print(f"  暴露: beta={beta_exp_col}, se={se_exp_col}")
    print(f"  结局: beta={beta_out_col}, se={se_out_col}")
    
    method_results = {}
    
    if beta_exp_col and se_exp_col and beta_out_col and se_out_col:
        try:
            merged[beta_exp_col] = pd.to_numeric(merged[beta_exp_col], errors="coerce")
            merged[se_exp_col] = pd.to_numeric(merged[se_exp_col], errors="coerce")
            merged[beta_out_col] = pd.to_numeric(merged[beta_out_col], errors="coerce")
            merged[se_out_col] = pd.to_numeric(merged[se_out_col], errors="coerce")
            
            valid = ~(np.isnan(merged[beta_exp_col]) | np.isnan(merged[se_exp_col]) |
                     np.isnan(merged[beta_out_col]) | np.isnan(merged[se_out_col]) |
                     (merged[se_exp_col] == 0) | (merged[se_out_col] == 0))
            
            valid_df = merged[valid].copy()
            n_valid = valid_df.shape[0]
            print(f"  有效SNP: {n_valid}")
            
            if n_valid >= 3:
                b_exp = valid_df[beta_exp_col].values
                b_out = valid_df[beta_out_col].values
                se_out = valid_df[se_out_col].values
                
                # IVW
                b_ivw, se_ivw, z_ivw, p_ivw = ivw_meta_analysis(b_exp, b_out, se_out)
                or_ivw = np.exp(b_ivw)
                or_l_ivw = np.exp(b_ivw - 1.96 * se_ivw)
                or_u_ivw = np.exp(b_ivw + 1.96 * se_ivw)
                
                # Simple median
                b_med, se_med, z_med, p_med = simple_mr(b_exp, b_out, se_out)
                
                # MR-Egger
                b_egg, se_egg, z_egg, p_egg = mr_egger(b_exp, b_out, se_out)
                
                method_results = {
                    "IVW": {"b": b_ivw, "se": se_ivw, "z": z_ivw, "p": p_ivw, "or": or_ivw, "or_l": or_l_ivw, "or_u": or_u_ivw, "n_snps": n_valid},
                    "Median": {"b": b_med, "se": se_med, "z": z_med, "p": p_med, "n_snps": n_valid},
                    "Egger": {"b": b_egg, "se": se_egg, "z": z_egg, "p": p_egg, "n_snps": n_valid} if b_egg is not None else {}
                }
                
                print(f"\n  === IVW结果 ===")
                print(f"  Beta={b_ivw:.4f}, SE={se_ivw:.4f}, Z={z_ivw:.3f}, P={p_ivw:.4e}")
                print(f"  OR={or_ivw:.4f} (95%CI: {or_l_ivw:.4f}-{or_u_ivw:.4f})")
                print(f"\n  === Median结果 ===")
                print(f"  Beta={b_med:.4f}, SE={se_med:.4f}, Z={z_med:.3f}, P={p_med:.4e}")
                print(f"\n  === Egger结果 ===")
                if b_egg is not None:
                    print(f"  Beta={b_egg:.4f}, SE={se_egg:.4f}, Z={z_egg:.3f}, P={p_egg:.4e}")
                
                status = "success"
            else:
                status = "insufficient_valid_snps"
                
        except Exception as e:
            print(f"  MR计算出错: {e}")
            import traceback; traceback.print_exc()
            status = f"error: {e}"
            method_results = {}
    else:
        print(f"  缺少beta/se列，无法计算MR")
        print(f"  暴露可用列: {[c for c in merged.columns if '_exp' in c]}")
        print(f"  结局可用列: {[c for c in merged.columns if '_out' in c]}")
        status = "missing_columns"
        # 打印所有列
        all_cols = list(merged.columns)
        print(f"  全部列名: {all_cols}")
    
    results_all.append({
        "exposure": "Perianal abscess",
        "outcome": name_cn,
        "exposure_gwas_id": exposure_id,
        "outcome_gwas_id": outcome_ids.get(name_cn),
        "status": status,
        "methods": method_results,
        "n_merged_snps": len(merged)
    })

# ============================================================
# 第四步：结果汇总
# ============================================================
print("\n" + "="*70)
print("第四步：汇总结果")
print("="*70)

# 打印汇总表
print("\n【IVW方法结果汇总】")
print("-"*100)
print(f"{'暴露':<20} {'结局':<12} {'GWAS_ID(暴露)':<15} {'GWAS_ID(结局)':<15} {'SNP数':<6} {'Beta':<10} {'SE':<10} {'OR':<10} {'95%CI':<22} {'P值':<14} {'状态'}")
print("-"*100)

for r in results_all:
    m = r.get("methods", {}).get("IVW", {})
    if m:
        or_str = f"{m.get('or', 0):.4f}"
        ci_str = f"({m.get('or_l', 0):.4f}-{m.get('or_u', 0):.4f})"
        p_str = f"{m.get('p', 1):.4e}" if m.get('p', 1) < 0.001 else f"{m.get('p', 1):.6f}"
        sig = "**显著**" if m.get('p', 1) < 0.05 else ""
        print(f"{'Perianal abscess':<20} {r['outcome']:<12} {r.get('exposure_gwas_id','NA'):<15} {r.get('outcome_gwas_id','NA'):<15} "
              f"{m.get('n_snps','NA'):<6} {m.get('b',0):<10.4f} {m.get('se',0):<10.4f} {or_str:<10} {ci_str:<22} {p_str:<14} {sig}")
    else:
        print(f"{'Perianal abscess':<20} {r['outcome']:<12} {r.get('exposure_gwas_id','NA'):<15} {r.get('outcome_gwas_id','NA'):<15} "
              f"{'NA':<6} {'NA':<10} {'NA':<10} {'NA':<10} {'NA':<22} {'NA':<14} {r.get('status','unknown')}")

print("\n【各方法详细P值】")
print("-"*70)
print(f"{'结局':<12} {'IVW_P':<12} {'Median_P':<12} {'Egger_P':<12} {'SNPs':<6}")
print("-"*70)
for r in results_all:
    m = r.get("methods", {})
    ivw_p = f"{m.get('IVW',{}).get('p',1):.4e}" if m.get('IVW') else "NA"
    med_p = f"{m.get('Median',{}).get('p',1):.4e}" if m.get('Median') else "NA"
    egg_p = f"{m.get('Egger',{}).get('p',1):.4e}" if m.get('Egger') else "NA"
    n_snps = m.get('IVW',{}).get('n_snps', r.get('n_merged_snps', 'NA'))
    print(f"{r['outcome']:<12} {ivw_p:<12} {med_p:<12} {egg_p:<12} {n_snps:<6}")

# 保存结果
results_df_data = []
for r in results_all:
    m = r.get("methods", {}).get("IVW", {})
    if m:
        results_df_data.append({
            "暴露": "Perianal abscess",
            "结局": r["outcome"],
            "暴露GWAS_ID": r.get("exposure_gwas_id", ""),
            "结局GWAS_ID": r.get("outcome_gwas_id", ""),
            "SNP数": m.get("n_snps", ""),
            "方法": "IVW",
            "Beta": round(m.get("b", 0), 4),
            "SE": round(m.get("se", 0), 4),
            "Z值": round(m.get("z", 0), 4),
            "P值": m.get("p", ""),
            "OR": round(m.get("or", 0), 4),
            "OR_95%CI_lower": round(m.get("or_l", 0), 4),
            "OR_95%CI_upper": round(m.get("or_u", 0), 4),
            "状态": r.get("status", "")
        })

if results_df_data:
    results_df = pd.DataFrame(results_df_data)
    results_df.to_csv(os.path.join(OUT_DIR, "MR_results_summary.csv"), index=False, encoding="utf-8-sig")
    print(f"\n结果CSV已保存: {os.path.join(OUT_DIR, 'MR_results_summary.csv')}")

# 保存JSON
with open(os.path.join(OUT_DIR, "MR_detailed_results.json"), "w", encoding="utf-8") as f:
    json.dump(results_all, f, ensure_ascii=False, indent=2, default=str)
print(f"详细JSON已保存: {os.path.join(OUT_DIR, 'MR_detailed_results.json')}")

# ============================================================
# 第五步：可视化
# ============================================================
print("\n" + "="*70)
print("第五步：生成可视化图表")
print("="*70)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# 左图：OR值森林图
ax1 = axes[0]
ivw_results = [(r["outcome"], r.get("methods", {}).get("IVW", {})) 
               for r in results_all if r.get("methods", {}).get("IVW")]

if ivw_results:
    y_pos = np.arange(len(ivw_results))
    ors = [m.get("or", 1) for _, m in ivw_results]
    or_ls = [m.get("or_l", 1) for _, m in ivw_results]
    or_us = [m.get("or_u", 1) for _, m in ivw_results]
    labels = [name for name, _ in ivw_results]
    pvals = [m.get("p", 1) for _, m in ivw_results]
    
    err_lo = [o - l for o, l in zip(ors, or_ls)]
    err_hi = [u - o for u, o in zip(or_us, ors)]
    
    colors = ["#c0392b" if p < 0.05 else "#2980b9" for p in pvals]
    
    bars = ax1.errorbar(ors, y_pos, xerr=[err_lo, err_hi], fmt="o",
                       color="black", capsize=6, markersize=9, linewidth=2, zorder=3)
    for i, (or_v, col) in enumerate(zip(ors, colors)):
        ax1.errorbar([ors[i]], [i], xerr=[[err_lo[i]], [err_hi[i]]], fmt="o",
                    color=col, capsize=6, markersize=9, linewidth=2, zorder=3)
    
    ax1.axvline(x=1, color="gray", linestyle="--", linewidth=2, zorder=2)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=13)
    ax1.set_xlabel("Odds Ratio (95% CI)", fontsize=13)
    ax1.set_title("Mendelian Randomization: Perianal Abscess →\nThromboembolic Diseases (IVW Method)",
                 fontsize=14, fontweight="bold")
    ax1.grid(axis="x", alpha=0.3)
    ax1.set_xlim(left=0)
    
    for i, (or_v, p, ol, ou) in enumerate(zip(ors, pvals, or_ls, or_us)):
        p_str = f"P={p:.2e}" if p < 0.001 else f"P={p:.4f}"
        ax1.text(max(or_us) * 1.15, i, f"OR={or_v:.2f}\n({ol:.2f}–{ou:.2f})\n{p_str}",
                va="center", ha="left", fontsize=9, 
                color="#c0392b" if p < 0.05 else "#2c3e50")
    
    ax1.text(1, -1.5, "Reference (OR=1)", ha="center", fontsize=9, color="gray")

# 右图：Beta值森林图
ax2 = axes[1]
if ivw_results:
    y_pos = np.arange(len(ivw_results))
    betas = [m.get("b", 0) for _, m in ivw_results]
    ses = [m.get("se", 0) for _, m in ivw_results]
    labels = [name for name, _ in ivw_results]
    pvals = [m.get("p", 1) for _, m in ivw_results]
    
    colors = ["#c0392b" if p < 0.05 else "#2980b9" for p in pvals]
    
    for i, (b, s, col) in enumerate(zip(betas, ses, colors)):
        ax2.errorbar([b], [i], xerr=[[1.96*s], [1.96*s]], fmt="o",
                    color=col, capsize=6, markersize=9, linewidth=2, zorder=3)
    
    ax2.axvline(x=0, color="gray", linestyle="--", linewidth=2, zorder=2)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=13)
    ax2.set_xlabel("Beta Coefficient (95% CI)", fontsize=13)
    ax2.set_title("Effect Size (Beta) Forest Plot\nPerianal Abscess MR Analysis",
                 fontsize=14, fontweight="bold")
    ax2.grid(axis="x", alpha=0.3)
    
    for i, (b, p, s) in enumerate(zip(betas, pvals, ses)):
        p_str = f"P={p:.2e}" if p < 0.001 else f"P={p:.4f}"
        x_txt = max(betas) + 0.5 if max(betas) > 0 else min(betas) - 0.5
        ax2.text(x_txt, i, f"β={b:.3f}±{s:.3f}\n{p_str}",
                va="center", ha="left", fontsize=9,
                color="#c0392b" if p < 0.05 else "#2c3e50")

plt.tight_layout()
fig_path = os.path.join(OUT_DIR, "MR_forest_plot.png")
plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()
print(f"森林图已保存: {fig_path}")

# 显著性汇总
print("\n" + "="*70)
print("【最终结论】")
print("="*70)
sig_found = []
for r in results_all:
    m = r.get("methods", {}).get("IVW", {})
    if m and m.get("p", 1) < 0.05:
        sig_found.append(r["outcome"])
        or_v = m.get("or", 0)
        direction = "增加" if or_v > 1 else "降低"
        print(f"✓ 显著关联: 肛周脓肿 {'↑' if or_v > 1 else '↓'} {r['outcome']} 风险 (OR={or_v:.3f}, P={m.get('p'):.4e})")

if not sig_found:
    print("未发现统计学显著的因果关联 (P<0.05)")
    print("注: 这可能由于肛周脓肿GWAS的SNP工具变量较少，统计效力不足")
    for r in results_all:
        m = r.get("methods", {}).get("IVW", {})
        if m:
            print(f"  - {r['outcome']}: OR={m.get('or',0):.3f}, P={m.get('p',1):.4f}, SNPs={m.get('n_snps','NA')}")

print("\n分析完成!")
print("="*70)
