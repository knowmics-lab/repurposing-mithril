import pandas as pd
import sys
filename=sys.argv[1]
df=pd.read_csv('MITHRIL/input/'+filename+'.mi2',sep=';',header=None)
print(df.iloc[0])
df.to_csv('MITHRIL/input/'+filename+'.mi',sep='\t', header=None, index=False )
