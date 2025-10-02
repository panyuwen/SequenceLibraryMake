#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: python pool_assist.py <input_file.xlsx> <output_file.txt>
输入文件（Excel）列（从第2行起，跳过表头）：
  <sample name> <qubit conc> <average size, /bp> <ratio, will be normalized> <lane ID> <aimed total volume for each lane, /ul>

主要逻辑：
  1) 按每条lane的 ratio 计算各样本的目标分子数 M_i = (aim_V * 10 nM) * (ratio_i / sum_ratio_lane)
  2) 直接用原液体积 v_raw = M_i / C_i（C_i为“Calculated nM”）
  3) 若 v_raw < MIN_PIPETTE_UL 且 C_i >= 10 nM，则做 d=ceil(MIN_PIPETTE_UL / v_raw) 倍预稀释，
     工作液浓度降为 C_i/d，实际取用体积 v_use = v_raw * d（>= MIN_PIPETTE_UL）
  4) 组内统一补水到 aim_V：water_for_group = max(0, aim_V - sum(v_use))
  5) 若 sum(v_use) > aim_V：
        - SCALE_WHEN_OVERFLOW=False（默认）：不缩放，提示溢出与建议（改总量或后浓缩）
        - SCALE_WHEN_OVERFLOW=True：等比例缩放到 aim_V，并提示可能出现 <MIN_PIPETTE_UL 的体积
"""

import pandas as pd
import sys
import os
import math

# ===== 可调参数 =====
MIN_PIPETTE_UL = 1.0          # 最小可移液体积（µL）
SCALE_WHEN_OVERFLOW = False   # 溢出时是否进行等比例缩放到 aim_V
ROUND_NDIGITS = 3             # 体积等输出保留的小数位

def ceil_div_volume(min_ul, v_raw):
    """返回使得 v_raw * d >= min_ul 的最小整数 d"""
    if v_raw <= 0:
        return 1
    return max(1, int(math.ceil(min_ul / v_raw)))

def main():
    if len(sys.argv) < 3:
        print("Usage: python generate_nM.py <input_file.xlsx> <output_file.txt>")
        sys.exit(1)

    filename = sys.argv[1]
    output = sys.argv[2]

    # 读取输入
    if filename.endswith('.xlsx'):
        data = pd.read_excel(filename, skiprows=1, header=None, usecols=[0,1,2,3,4,5])
    elif filename.endswith('.csv'):
        data = pd.read_csv(filename, skiprows=1, header=None, usecols=[0,1,2,3,4,5])
    else:
        data = pd.read_csv(filename, skiprows=1, header=None, usecols=[0,1,2,3,4,5], sep='\t')

    data.columns = ['Sample Name (Before Pooling)', "Qbit Con'C", 'Average size (bp)', 'ratio', 'group', 'aim_V']

    # 可选：补充 index 信息（与原脚本一致）
    data['index_name'] = data['Sample Name (Before Pooling)'].apply(lambda x: 'NEB-index' + x.split('_ID')[-1])
    if os.path.exists('Index primer, adaptor and Universal primer.xlsx'):
        index = pd.read_excel('Index primer, adaptor and Universal primer.xlsx', sheet_name='Sheet3')
        index['index_ID'] = index['Product'].apply(lambda x: x.split(' ')[-1])
        index = dict(zip(index['index_ID'], index['Index barcode']))
        data['Adapter 1 (I7_Index_ID) Sequence'] = data['Sample Name (Before Pooling)'].apply(
            lambda x: index.get(x.split('_ID')[-1])
        )
        data['Adapter 2 (I5_Index_ID) Sequence'] = 'universal'
    else:
        print('Index primer, adaptor and Universal primer.xlsx not found')

    # 计算分子浓度（nM）
    data['Calculated nM'] = data.apply(
        lambda x: x["Qbit Con'C"] / (660 * x['Average size (bp)']) * 1e6, axis=1
    )
    data['Calculated nM'] = data['Calculated nM'].round(3)

    # 每个 group 的 ratio 总和
    unit_count = data.groupby('group')['ratio'].sum().to_dict()

    # 目标分子数（用 nM*µL 表示，等价于纳摩尔·微升的单位；计算方便）
    # x['aim_V'] * 10.0 是 pool 总分子数
    data['Mi_nM_ul'] = data.apply(
        lambda x: (x['aim_V'] * 10.0) * (x['ratio'] / unit_count[x['group']]),
        axis=1
    )

    # 原液所需体积（µL）
    # 若 Calculated nM 很低，会得到较大的体积；若为0或缺失，置为 NaN 并提示
    def safe_divide(m, c):
        if pd.isna(c) or c <= 0:
            return float('nan')
        return m / c

    data['v_raw'] = data.apply(lambda x: safe_divide(x['Mi_nM_ul'], x['Calculated nM']), axis=1)

    # 预稀释判定（仅当 v_raw < 1µL 且 C>=10nM）
    def decide_predilute(row):
        v_raw = row['v_raw']
        C = row['Calculated nM']
        if pd.isna(v_raw):
            return pd.Series({'pre_dilute_required': 'NA', 'dilute_factor': 1, 'v_use': float('nan'), 'working_nM': float('nan')})
        if (v_raw < MIN_PIPETTE_UL) and (C >= 10.0):
            d = ceil_div_volume(MIN_PIPETTE_UL, v_raw)
            v_use = v_raw * d
            working_nM = C / d
            return pd.Series({
                'pre_dilute_required': 'Yes',
                'dilute_factor': d,
                'v_use': v_use,
                'working_nM': working_nM
            })
        else:
            # 不预稀释，直接使用原液
            return pd.Series({
                'pre_dilute_required': 'No',
                'dilute_factor': 1,
                'v_use': v_raw,
                'working_nM': C
            })

    extra = data.apply(decide_predilute, axis=1)
    data = pd.concat([data, extra], axis=1)

    # 分组统计混合体积、补水体积
    def group_stats_fn(df):
        aimV = df['aim_V'].iloc[0]
        V_mix = df['v_use'].sum(skipna=True)
        overflow = (V_mix > aimV)
        if overflow and SCALE_WHEN_OVERFLOW:
            scale = aimV / V_mix
            df['v_use_scaled'] = df['v_use'] * scale
            V_final = df['v_use_scaled'].sum(skipna=True)
            water = max(0.0, aimV - V_final)
            mode = 'scaled_to_aimV'
        else:
            df['v_use_scaled'] = df['v_use']  # 不缩放
            V_final = df['v_use_scaled'].sum(skipna=True)
            water = max(0.0, aimV - V_final)
            mode = 'no_scale_overflow' if overflow else 'fit_or_under'

        return df, pd.Series({
            'V_mix_before_scale': round(V_mix, ROUND_NDIGITS),
            'V_final': round(V_final, ROUND_NDIGITS),
            'aim_V': aimV,
            'water_for_group': round(water, ROUND_NDIGITS),
            'mode': mode
        })

    # 为输出文件清空/创建
    with open(output, 'w') as f:
        pass

    # 按 group 输出
    all_group_stats = []
    for group in data['group'].unique():
        df_g = data[data['group'] == group].copy()
        df_g, stats = group_stats_fn(df_g)
        all_group_stats.append((group, stats))

        # 打印分组头
        with open(output, 'a') as f:
            f.write('=' * 10 + '\n')
            f.write(f'group {group} | aim_V={stats["aim_V"]} µL | V_mix(before scale)={stats["V_mix_before_scale"]} µL | V_final={stats["V_final"]} µL | mode={stats["mode"]}\n')
            f.write(f'Add water: {stats["water_for_group"]} µL\n')

        # 预稀释样本
        tmp_yes = df_g[df_g['pre_dilute_required'] == 'Yes'].copy()
        if not tmp_yes.empty:
            tmp_yes = tmp_yes[['Sample Name (Before Pooling)', 'Calculated nM', 'Mi_nM_ul', 'v_raw', 'dilute_factor', 'working_nM', 'v_use_scaled']]
            tmp_yes = tmp_yes.rename(columns={
                'v_use_scaled': 'final_volume_to_use(µL)'
            })
            tmp_yes = tmp_yes.round(ROUND_NDIGITS)
            with open(output, 'a') as f:
                f.write('--- pre-dilute: Yes (prepare working solution; then take final volume) ---\n')
                f.write('Columns: Sample\tC(nM)\tM(nM*µL)\tv_raw(µL)\tdilute_factor\tworking_nM\ttake(µL)\n')
            tmp_yes.to_csv(output, sep='\t', mode='a', header=False, index=False)

        # 不预稀释样本
        tmp_no = df_g[df_g['pre_dilute_required'] == 'No'].copy()
        if not tmp_no.empty:
            tmp_no = tmp_no[['Sample Name (Before Pooling)', 'Calculated nM', 'Mi_nM_ul', 'v_use_scaled']]
            tmp_no = tmp_no.rename(columns={
                'v_use_scaled': 'final_volume_to_use(µL)'
            })
            tmp_no = tmp_no.round(ROUND_NDIGITS)
            with open(output, 'a') as f:
                f.write('--- pre-dilute: No (use stock directly) ---\n')
                f.write('Columns: Sample\tC(nM)\tM(nM*µL)\ttake(µL)\n')
            tmp_no.to_csv(output, sep='\t', mode='a', header=False, index=False)

        # 无法计算（例如浓度为0或缺失）
        tmp_na = df_g[df_g['pre_dilute_required'] == 'NA'].copy()
        if not tmp_na.empty:
            tmp_na = tmp_na[['Sample Name (Before Pooling)', 'Calculated nM', 'Mi_nM_ul', 'v_raw']]
            tmp_na = tmp_na.round(ROUND_NDIGITS)
            with open(output, 'a') as f:
                f.write('--- WARNING: Missing/zero concentration; cannot compute volume ---\n')
                f.write('Columns: Sample\tC(nM)\tM(nM*µL)\tv_raw(µL)\n')
            tmp_na.to_csv(output, sep='\t', mode='a', header=False, index=False)

    # 追加一个总览表
    with open(output, 'a') as f:
        f.write('=' * 10 + '\n')
        f.write('=== Group Summary ===\n')
        f.write('group\tmode\taim_V(µL)\tV_mix_before_scale(µL)\tV_final(µL)\tAdd_water(µL)\n')
    for group, stats in all_group_stats:
        with open(output, 'a') as f:
            f.write(f'{group}\t{stats["mode"]}\t{stats["aim_V"]}\t{stats["V_mix_before_scale"]}\t{stats["V_final"]}\t{stats["water_for_group"]}\n')

    # 最后把详表也贴到文件末尾（便于排查/记录）
    printable = data.copy()
    for col in ['v_raw', 'v_use', 'v_use_scaled', 'working_nM', 'Mi_nM_ul', 'Calculated nM']:
        if col in printable.columns:
            printable[col] = printable[col].round(ROUND_NDIGITS)
    printable.sort_values(by=['group', 'pre_dilute_required', 'Sample Name (Before Pooling)'], inplace=True)
    with open(output, 'a') as f:
        f.write('=' * 10 + '\n')
    printable.to_csv(output, sep='\t', mode='a', index=False)

if __name__ == '__main__':
    main()