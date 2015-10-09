for i in `list_pnfn.txt`
do
    echo "$i" 
    dccp -H "$i" crab_MC_G_HT-100to200
done
