import sys
import pandas as pd

data = sys.argv[1]
sample_info = sys.argv[2]
out = sys.argv[3]

df_data = pd.read_csv(data, sep='\t', index_col=0)
df_data = df_data.T
# print(df_data.index)
df_sampleinfo = pd.read_csv(sample_info, sep='\t', index_col=0)
df_sampleinfo = pd.DataFrame(df_sampleinfo['Group'])
df_sampleinfo.index = [str(i) for i in df_sampleinfo.index]
# print(df_sampleinfo.index)

new_df = pd.concat([df_sampleinfo, df_data], axis=1, sort=False)
new_df.index.names = ['SampleID']
new_df.to_csv(out, sep='\t')
