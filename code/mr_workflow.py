"""
肛周脓肿 与 血栓栓塞性疾病的孟德尔随机化分析
数据来源: IEU OpenGWAS 数据库 (https://gwas.mrcieu.ac.uk/)
分析工具: Python (TwoSampleMR, ieugwaspy)
"""

import os
import sys
import json
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np

# 设置matplotlib中文字体
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Microsoft YaHei, SimHei, Arial Unicode MS, DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

# 尝试导入MR分析包
try:
    from ieugwaspy import GWAS
    HAS_IEU = True
    print("ieugwaspy 导入成功")
except ImportError:
    HAS_IEU = False
    print("警告: ieugwaspy 未安装，部分功能可能受限")

try:
    import TwoSampleMR
    from TwoSampleMR import ao, ld, radim
    HAS_TWOSAMPLEMR = True
    print("TwoSampleMR 导入成功")
except ImportError:
    HAS_TWOSAMPLEMR = False
    print("警告: TwoSampleMR 未安装")

# ============================================================
# 配置
# ============================================================
OUT_DIR = r"D:\肛周脓肿与血栓性疾病"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# 第一步：通过IEU API搜索GWAS数据ID
# ============================================================
print("\n" + "="*70)
print("第一步：搜索 IEU OpenGWAS 数据库中的GWAS数据")
print("="*70)

def search_gwas_keyword(keyword, poolsize=1000):
    """使用ieugwaspy搜索GWAS摘要信息"""
    try:
        gdainfo = GWAS().gdainfo(keyword=keyword, poolsize=poolsize)
        df = pd.DataFrame(gdainfo)
        return df
    except Exception as e:
        print(f"  搜索关键词 '{keyword}' 出错: {e}")
        return pd.DataFrame()

def pick_best_gwas(df, priority_id=None):
    """从搜索结果中选择最优GWAS数据集"""
    if df.empty:
        return None
    cols = ["id", "trait", "year", "nsnp", "sample_size", "author", "unit"]
    avail = [c for c in cols if c in df.columns]
    # 优先选择有SNP数量且样本量最大的
    if "nsnp" in avail and "sample_size" in avail:
        df2 = df.copy()
        df2["nsnp"] = pd.to_numeric(df2["nsnp"], errors="coerce")
        df2["sample_size"] = pd.to_numeric(df2["sample_size"], errors="coerce")
        df2 = df2.dropna(subset=["nsnp", "sample_size"])
        df2 = df2.sort_values(["nsnp", "sample_size"], ascending=[False, False])
        print(df2[avail].head(10).to_string())
        return df2.iloc[0]["id"] if len(df2) > 0 else None
    print(df[avail].head(10).to_string())
    return df.iloc[0]["id"] if len(df) > 0 else None

# 搜索暴露
print("\n--- 搜索暴露：Perianal abscess ---")
exp_df = search_gwas_keyword("perianal abscess")
if exp_df.empty:
    print("  未找到结果，尝试 'anal abscess'")
    exp_df = search_gwas_keyword("anal abscess")
if exp_df.empty:
    print("  尝试 'anorectal disease'")
    exp_df = search_gwas_keyword("anorectal disease")

exposure_id = pick_best_gwas(exp_df) if not exp_df.empty else None
print(f"\n选择的暴露GWAS ID: {exposure_id}")

# 搜索结局
outcome_configs = {
    "脑梗": ["cerebral infarction", "ischemic stroke", "stroke"],
    "心梗": ["myocardial infarction", "coronary heart disease", "heart attack"],
    "深静脉血栓": ["deep vein thrombosis", "venous thrombosis", "venous thromboembolism"],
    "动脉栓塞": ["arterial embolism", "peripheral arterial disease", "arterial thrombosis"]
}

outcome_ids = {}
for name_cn, keywords in outcome_configs.items():
    print(f"\n--- 搜索结局：{name_cn} ---")
    chosen_id = None
    for kw in keywords:
        df2 = search_gwas_keyword(kw)
        if not df2.empty:
            chosen_id = pick_best_gwas(df2)
            if chosen_id:
                print(f"  使用关键词 '{kw}' 选择ID: {chosen_id}")
                break
    outcome_ids[name_cn] = chosen_id
    if not chosen_id:
        print(f"  警告：未找到 {name_cn} 的合适GWAS数据")

print("\n搜索完成。")
print(f"暴露ID: {exposure_id}")
print("结局ID:", json.dumps(outcome_ids, ensure_ascii=False, indent=2))

# ============================================================
# 第二步：提取SNP信息
# ============================================================
print("\n" + "="*70)
print("第二步：从选中GWAS中提取SNP工具变量")
print("="*70)

def extract_snps_from_gwas(gwas_id, pval_threshold=5e-8):
    """从GWAS中提取全基因组显著SNP"""
    try:
        # 使用ieugwaspy获取变异位点信息
        variants = GWAS().gwasinfo(gwas_id=gwas_id)
        if variants is None or (isinstance(variants, dict) and not variants):
            print(f"  无法获取GWAS {gwas_id} 的数据")
            return pd.DataFrame()
        df = pd.DataFrame(variants)
        print(f"  GWAS {gwas_id} 共有 {len(df)} 个变异位点")
        # 尝试找p值列
        pcol = None
        for col in ["p", "pval", "pvalue", "p_value"]:
            if col in df.columns:
                pcol = col
                break
        if pcol:
            df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
            sig = df[df[pcol] < pval_threshold]
            print(f"  P<5e-8 的显著SNP: {len(sig)} 个")
            return sig
        print(f"  可用列: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"  提取SNP失败: {e}")
        return pd.DataFrame()

if exposure_id and HAS_IEU:
    print(f"\n提取暴露 {exposure_id} 的显著SNP...")
    exp_snps = extract_snps_from_gwas(exposure_id)
    print(f"暴露SNP数量: {len(exp_snps)}")
else:
    exp_snps = pd.DataFrame()
    print("跳过暴露SNP提取")

for name_cn, oid in outcome_ids.items():
    if oid and HAS_IEU:
        print(f"\n提取结局 {name_cn} ({oid}) 的显著SNP...")
        out_snps = extract_snps_from_gwas(oid)
        print(f"结局SNP数量: {len(out_snps)}")
    else:
        print(f"\n{name_cn} 无可用GWAS ID")

# ============================================================
# 第三步：孟德尔随机化分析 (基于已有GWAS数据)
# ============================================================
print("\n" + "="*70)
print("第三步：执行孟德尔随机化分析")
print("="*70)

# 保存汇总结果
results_all = []

def run_mr_analysis(exp_id, out_id, exp_name, out_name):
    """对一对暴露-结局执行MR分析"""
    print(f"\n>>> MR分析: {exp_name} -> {out_name}")
    print(f"    暴露: {exp_id}")
    print(f"    结局: {out_id}")
    
    try:
        # 使用TwoSampleMR的标准分析流程
        # 获取暴露和结局的汇总统计数据
        exp_data = ao.get_metadata(exp_id)
        if exp_data is None or (isinstance(exp_data, (list, dict)) and len(exp_data) == 0):
            print(f"    无法获取暴露 {exp_id} 元数据")
            return None
        
        print(f"    暴露元数据获取成功")
        
        # 尝试从数据库获取SNP信息
        # 先用exp_id搜索相关SNP
        try:
            # 获取暴露的SNP列表
            exp_variants = GWAS().gwasinfo(gwas_id=exp_id)
            if exp_variants is not None and len(exp_variants) > 0:
                exp_df = pd.DataFrame(exp_variants)
                print(f"    暴露 {exp_id}: {len(exp_df)} 个变异")
            else:
                print(f"    暴露 {exp_id} 无变异数据")
                exp_df = pd.DataFrame()
        except Exception as e:
            print(f"    获取暴露变异出错: {e}")
            exp_df = pd.DataFrame()
        
        # 获取结局的SNP
        try:
            out_variants = GWAS().gwasinfo(gwas_id=out_id)
            if out_variants is not None and len(out_variants) > 0:
                out_df = pd.DataFrame(out_variants)
                print(f"    结局 {out_id}: {len(out_df)} 个变异")
            else:
                print(f"    结局 {out_id} 无变异数据")
                out_df = pd.DataFrame()
        except Exception as e:
            print(f"    获取结局变异出错: {e}")
            out_df = pd.DataFrame()
        
        # 检查是否有足够的重叠SNP进行分析
        if exp_df.empty or out_df.empty:
            print(f"    数据不完整，跳过此配对分析")
            return None
        
        # 尝试找SNP重叠
        snp_col = None
        for col in ["snp", "variant", "rsid", "id"]:
            if col in exp_df.columns and col in out_df.columns:
                snp_col = col
                break
        
        if snp_col is None:
            print(f"    找不到SNP标识列")
            return None
        
        # 找到重叠SNP
        common_snps = set(exp_df[snp_col]) & set(out_df[snp_col])
        print(f"    重叠SNP数量: {len(common_snps)}")
        
        if len(common_snps) < 3:
            print(f"    重叠SNP数量不足(<3)，使用更多SNP...")
            # 使用更多的SNP（放宽p值阈值）
            pcol_exp = None
            for col in ["p", "pval", "pvalue"]:
                if col in exp_df.columns:
                    pcol_exp = col
                    break
            if pcol_exp:
                try:
                    exp_df[pcol_exp] = pd.to_numeric(exp_df[pcol_exp], errors="coerce")
                    # 逐步放宽阈值
                    for thresh in [1e-5, 1e-3, 0.05]:
                        exp_top = exp_df[exp_df[pcol_exp] < thresh]
                        if len(exp_top) >= 3:
                            print(f"    使用 p<{thresh} 的 {len(exp_top)} 个SNP")
                            break
                except:
                    pass
        
        # 构建分析用的数据框
        # 合并暴露和结局数据
        merged = pd.merge(exp_df, out_df, on=snp_col, suffixes=("_exp", "_out"))
        print(f"    合并后SNP数量: {len(merged)}")
        
        if len(merged) < 3:
            print(f"    合并后SNP数量不足，跳过")
            return None
        
        # 计算MR
        # 需要beta和se列
        beta_col_exp = None
        se_col_exp = None
        beta_col_out = None
        se_col_out = None
        
        for col in merged.columns:
            cn = col.lower()
            if "_exp" in col:
                if "beta" in cn or "effect" in cn:
                    beta_col_exp = col
                if "se" in cn or "stderr" in cn:
                    se_col_exp = col
            if "_out" in col:
                if "beta" in cn or "effect" in cn:
                    beta_col_out = col
                if "se" in cn or "stderr" in cn:
                    se_col_out = col
        
        print(f"    暴露beta列: {beta_col_exp}, se列: {se_col_exp}")
        print(f"    结局beta列: {beta_col_out}, se列: {se_col_out}")
        
        # 如果没有找到合适的beta/se，尝试其他方式
        if not beta_col_exp or not se_col_exp:
            # 尝试从列名推断
            for col in merged.columns:
                if "exp" in col.lower():
                    print(f"    暴露列: {col} -> {merged[col].dtype}")
        if not beta_col_out or not se_col_out:
            for col in merged.columns:
                if "out" in col.lower():
                    print(f"    结局列: {col} -> {merged[col].dtype}")
        
        # MR分析 - 使用基本方法
        result_entry = {
            "exposure": exp_name,
            "outcome": out_name,
            "exposure_gwas_id": exp_id,
            "outcome_gwas_id": out_id,
            "n_snps": len(merged),
            "method": "IVW",
            "b": None,
            "se": None,
            "pval": None,
            "or": None,
            "or_lower": None,
            "or_upper": None,
            "status": "insufficient_data"
        }
        
        # 尝试提取beta和se（如果有的话）
        if beta_col_exp and se_col_exp and beta_col_out and se_col_out:
            try:
                merged[beta_col_exp] = pd.to_numeric(merged[beta_col_exp], errors="coerce")
                merged[se_col_exp] = pd.to_numeric(merged[se_col_exp], errors="coerce")
                merged[beta_col_out] = pd.to_numeric(merged[beta_col_out], errors="coerce")
                merged[se_col_out] = pd.to_numeric(merged[se_col_out], errors="coerce")
                
                # 简单的IVW方法
                b_exp = merged[beta_col_exp].values
                se_exp = merged[se_col_exp].values
                b_out = merged[beta_col_out].values
                se_out = merged[se_col_out].values
                
                valid = ~(np.isnan(b_exp) | np.isnan(se_exp) | np.isnan(b_out) | np.isnan(se_out) | (se_exp == 0) | (se_out == 0))
                if valid.sum() >= 3:
                    b_exp = b_exp[valid]
                    se_exp = se_exp[valid]
                    b_out = b_out[valid]
                    se_out = se_out[valid]
                    
                    # IVW meta-analysis
                    # b_mr = sum(b_out/se_out^2) / sum(1/se_out^2)
                    # var_b_mr = 1 / sum(1/se_out^2)
                    weights = 1.0 / (se_out ** 2)
                    b_ivw = np.sum(weights * b_out) / np.sum(weights)
                    var_ivw = 1.0 / np.sum(weights)
                    se_ivw = np.sqrt(var_ivw)
                    z = b_ivw / se_ivw
                    p_ivw = 2 * (1 - stats.norm.cdf(abs(z)))
                    
                    # 计算OR值
                    or_ivw = np.exp(b_ivw)
                    or_lower = np.exp(b_ivw - 1.96 * se_ivw)
                    or_upper = np.exp(b_ivw + 1.96 * se_ivw)
                    
                    result_entry.update({
                        "method": "IVW (fixed-effects)",
                        "b": round(b_ivw, 4),
                        "se": round(se_ivw, 4),
                        "pval": p_ivw,
                        "or": round(or_ivw, 4),
                        "or_lower": round(or_lower, 4),
                        "or_upper": round(or_upper, 4),
                        "status": "success"
                    })
                    print(f"    IVW结果: b={b_ivw:.4f}, SE={se_ivw:.4f}, p={p_ivw:.4e}")
                    print(f"    OR={or_ivw:.4f} (95%CI: {or_lower:.4f}-{or_upper:.4f})")
            except Exception as e:
                print(f"    MR计算出错: {e}")
        
        return result_entry
        
    except Exception as e:
        print(f"    MR分析过程出错: {e}")
        import traceback
        traceback.print_exc()
        return {
            "exposure": exp_name,
            "outcome": out_name,
            "exposure_gwas_id": exp_id,
            "outcome_gwas_id": out_id,
            "status": f"error: {str(e)}"
        }

# 对每个结局执行MR分析
for name_cn, out_id in outcome_ids.items():
    if exposure_id and out_id:
        res = run_mr_analysis(exposure_id, out_id, "Perianal abscess", name_cn)
        if res:
            results_all.append(res)
    else:
        print(f"\n跳过 {name_cn}: 缺少GWAS ID (暴露={exposure_id}, 结局={out_id})")

# ============================================================
# 第四步：结果汇总与可视化
# ============================================================
print("\n" + "="*70)
print("第四步：汇总分析结果")
print("="*70)

# 打印结果表格
print("\n【孟德尔随机化分析结果汇总】")
print("-" * 90)
header = f"{'暴露':<20} {'结局':<12} {'方法':<20} {'SNP数':<6} {'Beta':<10} {'SE':<10} {'OR':<10} {'95%CI':<20} {'P值':<12}"
print(header)
print("-" * 90)

results_df = pd.DataFrame(results_all)
for _, row in results_all:
    status = row.get("status", "unknown")
    if status == "success":
        or_ci = f"{row.get('or', 'NA')} ({row.get('or_lower', '')}-{row.get('or_upper', '')})"
        print(f"{'Perianal abscess':<20} {row['outcome']:<12} {row['method']:<20} {row.get('n_snps', 'NA'):<6} "
              f"{row.get('b', 'NA'):<10} {row.get('se', 'NA'):<10} {row.get('or', 'NA'):<10} {or_ci:<20} {row.get('pval', 'NA'):<12.4e}")
    else:
        print(f"{'Perianal abscess':<20} {row['outcome']:<12} {'--':<20} {'--':<6} {'--':<10} {'--':<10} {'--':<10} {'--':<20} {status:<12}")

# 保存结果
results_df.to_csv(os.path.join(OUT_DIR, "MR_results_summary.csv"), index=False, encoding="utf-8-sig")
print(f"\n结果已保存: {os.path.join(OUT_DIR, 'MR_results_summary.csv')}")

# 绘制森林图
if results_all:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 左图：OR值森林图
    ax1 = axes[0]
    success_results = [r for r in results_all if r.get("status") == "success"]
    
    if success_results:
        y_pos = np.arange(len(success_results))
        ors = [r.get("or", 1) for r in success_results]
        lowers = [r.get("or_lower", 1) for r in success_results]
        uppers = [r.get("or_upper", 1) for r in success_results]
        labels = [r["outcome"] for r in success_results]
        
        errors = [[o - l for o, l in zip(ors, lowers)],
                  [u - o for u, o in zip(uppers, ors)]]
        
        colors = ["#e74c3c" if r.get("pval", 1) < 0.05 else "#3498db" for r in success_results]
        
        ax1.errorbar(ors, y_pos, xerr=errors, fmt="o", color="black",
                    capsize=5, markersize=8, linewidth=2)
        for i, (col, or_val) in enumerate(zip(colors, ors)):
            ax1.errorbar([or_val], [y_pos[i]], xerr=[[errors[0][i]], [errors[1][i]]],
                        fmt="o", color=col, capsize=5, markersize=8, linewidth=2)
        
        ax1.axvline(x=1, color="gray", linestyle="--", linewidth=1.5)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(labels, fontsize=12)
        ax1.set_xlabel("Odds Ratio (95% CI)", fontsize=12)
        ax1.set_title("Perianal Abscess -> Thromboembolic Outcomes\nMendelian Randomization", fontsize=13)
        ax1.grid(axis="x", alpha=0.3)
        
        # 添加OR值标注
        for i, r in enumerate(success_results):
            or_txt = f"OR={r['or']:.2f} ({r['or_lower']:.2f}-{r['or_upper']:.2f})"
            p_txt = f"P={r['pval']:.3e}" if r.get("pval") < 0.001 else f"P={r['pval']:.3f}"
            ax1.text(max(uppers) * 1.1, i, f"{or_txt}\n{p_txt}", va="center", fontsize=9)
    
    # 右图：Beta值森林图
    ax2 = axes[1]
    if success_results:
        y_pos = np.arange(len(success_results))
        betas = [r.get("b", 0) for r in success_results]
        ses = [r.get("se", 0) for r in success_results]
        labels = [r["outcome"] for r in success_results]
        
        colors = ["#e74c3c" if r.get("pval", 1) < 0.05 else "#3498db" for r in success_results]
        
        for i, (b, s, col) in enumerate(zip(betas, ses, colors)):
            ax2.errorbar([b], [i], xerr=[[1.96*s], [1.96*s]], fmt="o",
                        color=col, capsize=5, markersize=8, linewidth=2)
        
        ax2.axvline(x=0, color="gray", linestyle="--", linewidth=1.5)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(labels, fontsize=12)
        ax2.set_xlabel("Beta (95% CI)", fontsize=12)
        ax2.set_title("Effect Size (Beta) Forest Plot\nPerianal Abscess MR Analysis", fontsize=13)
        ax2.grid(axis="x", alpha=0.3)
    
    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "MR_forest_plot.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"森林图已保存: {fig_path}")

print("\n分析完成!")
print("="*70)

# 最终状态摘要
print("\n【分析状态摘要】")
if results_all:
    success_count = sum(1 for r in results_all if r.get("status") == "success")
    print(f"  总分析配对数: {len(results_all)}")
    print(f"  成功完成: {success_count}")
    print(f"  失败/数据不足: {len(results_all) - success_count}")
    for r in results_all:
        status_mark = "✓" if r.get("status") == "success" else "✗"
        print(f"  {status_mark} {r['exposure']} -> {r['outcome']}: {r.get('status', 'unknown')}")
else:
    print("  警告: 未能完成任何MR分析配对")

# 保存详细结果JSON
with open(os.path.join(OUT_DIR, "MR_detailed_results.json"), "w", encoding="utf-8") as f:
    json.dump(results_all, f, ensure_ascii=False, indent=2, default=str)

print(f"\n所有结果文件:")
print(f"  1. {os.path.join(OUT_DIR, 'MR_results_summary.csv')}")
print(f"  2. {os.path.join(OUT_DIR, 'MR_forest_plot.png')}")
print(f"  3. {os.path.join(OUT_DIR, 'MR_detailed_results.json')}")
