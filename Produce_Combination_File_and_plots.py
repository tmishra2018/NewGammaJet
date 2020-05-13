#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime

parser = argparse.ArgumentParser(description='Create combination Files and corresponding plots')
parser.add_argument('--Pu_profile',type=str,dest="Pu_profile",default=1, help="Pu_profile")
parser.add_argument('--Input_data',type=str,dest="Input_data",default=1, help="directory where to find the lists")
parser.add_argument('--Input_mc',type=str,dest="Input_mc",default=1, help="directory where to find the lists")
parser.add_argument('--output',type=str,dest="output",default=1, help="output suffix")
parser.add_argument('--IOV',type=str,dest="IOV",default=1, help="IOV suffix")
args = parser.parse_args()
Pu_profile = args.Pu_profile
Input_data = args.Input_data
Input_mc   = args.Input_mc
output     = args.output
IOV        = args.IOV
today = datetime.date.today()
today.strftime('%d-%m-%Y')

List_alpha = [0.3,0.25,0.20,0.15,0.1]
List_alpha_name = ['0_3','0_25','0_20','0_15','0_1']
i_alpha = 0
cmdFinal = "hadd -f Gjet_combinationfile_"+output+"_"+IOV+".root"

for alpha in List_alpha:
	
	cmd1 = "gammaJetFinalizer -i "+Input_data+" -d DATA_"+output+"_"+List_alpha_name[i_alpha]+"_"+str(today)+" --type pf --algo ak4  --runera "+Pu_profile+" --alpha "+str(alpha)
	dataname = "DATA_"+output+"_"+List_alpha_name[i_alpha]+"_"+str(today)+"_"+Pu_profile
	os.system(cmd1)
	
	
	cmd2 = "gammaJetFinalizer --input-list "+Input_mc+" -d MC_"+output+"_"+List_alpha_name[i_alpha]+"_"+str(today)+" --type pf --algo ak4  --runera "+Pu_profile+" --mc --alpha "+str(alpha)
	mcname = "MC_"+output+"_"+List_alpha_name[i_alpha]+"_"+str(today)+"_"+Pu_profile
	os.system(cmd2)
	
	cmd3 = "./analysis/draw/drawPhotonJet_2bkg "+dataname+" "+mcname+" "+mcname+"  pf ak4 LUMI"
	os.system(cmd3)
	
	cmd4 = "./analysis/draw/drawPhotonJetExtrap --type pf --algo ak4 "+dataname+" "+mcname+" "+mcname  
	os.system(cmd4)
	
	cmd5 = "./analysis/draw/draw_ratios_vs_pt "+dataname+" "+mcname+" "+mcname+"  pf ak4"
	os.system(cmd5)
	
	cmd6 = "./analysis/draw/drawFlavorFractions "+mcname+"  pf ak4"
	os.system(cmd6)
	
	cmdFinal +=" PhotonJetPlots_"+dataname+"_vs_"+mcname+"_PFlowAK4chs_LUMI_vs_pt/plots.root PhotonJetPlots_"+mcname+"_PFlowAK4chs_FlavorFractions/*.root"
	
	#cmd6 = "analysis/draw/draw_all_methods_vs_pt "+dataname+" "+mcname+" "+mcname+"  pf ak4"
	#os.system(cmd6)
	
	i_alpha += 1

print(cmdFinal)
os.system(cmdFinal)
	
