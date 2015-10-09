echo Compilo PUWeight_NvtxBased.c
g++ -o PUWeight_NvtxBased.exe PUWeight_NvtxBased.c `root-config --cflags  --glibs`

rm *~