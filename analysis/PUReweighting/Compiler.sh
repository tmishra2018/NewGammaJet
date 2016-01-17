echo Compilo PUWeight_NvtxBased.c
g++ -o PUWeight_NvtxBased.exe PUWeight_NvtxBased.c `root-config --cflags  --glibs`
echo Compilo generate_mc_pileup.c
g++ -o generate_mc_pileup.exe generate_mc_pileup.c `root-config --cflags  --glibs`
echo Compilo Ratio_Test.c
g++ -o Ratio_Test.exe Ratio_Test.c `root-config --cflags  --glibs`

rm *~