import os
from os.path import join, getsize
for root, dirs, files in os.walk('/pnfs/roma1.infn.it/data/cms/store/user/fpreiato/JetMET/MC/01Jun15/G_HT-100to200/'):
    print(root) 
    print(dirs) 
    print(files)
    for name in files:
        print(name) 
        dccp -H root/name crab_MC_G_HT-100to200/

 #       done 
#    os.path.join(dirpath, name)
#    print(root, "consumes", end=" ")
#    print(sum(getsize(join(root, name)) for name in files), end=" ")

#fede        print("bytes in", len(names), "non-directory files")
        
        



