"""诊断 IEUGWASPY + TwoSampleMR API 连通性"""
import sys, os

# 确保用户 site-packages 在路径最前面
user_pkgs = r'C:\Users\zhushengm\AppData\Roaming\Python\Python311\site-packages'
if user_pkgs not in sys.path:
    sys.path.insert(0, user_pkgs)

os.chdir(r'D:\肛周脓肿与血栓性疾病')

print('=== 1. API 状态检查 ===')
from ieugwaspy import config, query
print('Base URL:', config.env['base_url'])
try:
    st = query.status()
    print(f"API 版本: {st['API__VERSION']} | 数据集数: {st['STATISTICS__N_DATASETS']}")
except Exception as e:
    print('API 状态失败:', e)

print('\n=== 2. JWT Token 检查 ===')
try:
    from TwoSampleMR.verifytoken import get_opengwas_jwt
    jwt = get_opengwas_jwt()
    print('JWT token:', (jwt[:20] + '...') if jwt else '未设置')
except Exception as e:
    print('JWT 错误:', e)

print('\n=== 3. 测试 ieugwaspy.gwasinfo (已知 ID) ===')
try:
    res = query.gwasinfo(['ieu-a-2'])
    if isinstance(res, dict) and 'message' in res:
        print('需要 JWT:', res['message'][:100])
    else:
        print('结果类型:', type(res))
        if isinstance(res, dict):
            for k, v in list(res.items())[:2]:
                print(f'  {k}: {str(v)[:120]}')
except Exception as e:
    print('gwasinfo 失败:', e)

print('\n=== 4. available_outcomes (过滤血栓/脓肿相关) ===')
try:
    df = query.available_outcomes()
    print(f'总数据集: {df.shape[0]}')
    keywords = ['abscess','stroke','myocardial','venous','thrombosis','embolism','arterial','infarction']
    cn_names = {'abscess':'肛周脓肿','stroke':'脑梗/卒中','myocardial':'心梗','venous':'深静脉血栓',
                'thrombosis':'血栓','embolism':'栓塞','arterial':'动脉','infarction':'梗死'}
    for kw in keywords:
        sub = df[df['trait'].str.lower().str.contains(kw, na=False)]
        cn = cn_names.get(kw, kw)
        print(f'\n  [{cn}] ({kw}): {len(sub)} 个')
        if not sub.empty:
            print(sub[['id','trait','nsnp','sample_size','author','year']].head(3).to_string())
except Exception as e:
    print('available_outcomes 失败:', e)
    import traceback; traceback.print_exc()
