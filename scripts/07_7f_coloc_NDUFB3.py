import pandas as pd
from liftover import get_lifter

df = pd.read_csv("E:/retina_tools/aging_rgc_axon/AMD/results/07_7f_coloc_NDUFB3/ndufb3_gwas_region.csv")
print(f"读入：{len(df)} 个 SNP")
print(df.head())

converter = get_lifter('hg19', 'hg38', one_based=True)

hg38_pos = []
for _, row in df.iterrows():
    chrom = str(int(row['chr']))
    pos = int(row['pos'])
    result = converter[chrom][pos]
    if result:
        hg38_pos.append(result[0][1])
    else:
        hg38_pos.append(None)

df['pos_hg38'] = hg38_pos
df_hg38 = df.dropna(subset=['pos_hg38'])
print(f"转换成功：{len(df_hg38)} / {len(df)} 个 SNP")

colstats = df_hg38[['chr', 'pos_hg38', 'ES', 'SE']].copy()
colstats.columns = ['chrom', 'pos', 'beta', 'se']
colstats.to_csv("E:/retina_tools/aging_rgc_axon/AMD/results/07_7f_coloc_NDUFB3NDUFB3_earlyAMD_COLSTATS_hg38.tsv", sep='\t', index=False)
print("已保存：E:/retina_tools/aging_rgc_axon/AMD/results/07_7f_coloc_NDUFB3")