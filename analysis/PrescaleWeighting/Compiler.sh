echo Compilo MC_ptPhot_scaled.c
g++ -o MC_ptPhot_scaled.exe MC_ptPhot_scaled.c `root-config --cflags  --glibs`
echo Compilo PrescaleWeight_After.c
g++ -o PrescaleWeight_After.exe PrescaleWeight_After.c `root-config --cflags  --glibs`
echo Compilo PrescaleWeight_Before.c
g++ -o PrescaleWeight_Before.exe PrescaleWeight_Before.c `root-config --cflags  --glibs`
