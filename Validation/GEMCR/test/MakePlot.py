#!/usr/bin/env python

"""
Copied from GEMCode/GEMValidation
"""

from ROOT import *

import os
import sys
#sys.argv.append( '-b' )
import ROOT
from array import array
ROOT.gROOT.SetBatch(1)

ablue = array("d", [0.51, 1.00, 0.12, 0.00, 0.00])
ared = array("d", [0.00, 0.00, 0.87, 1.00, 0.51])
agreen = array("d", [0.00, 0.81, 1.00, 0.20, 0.00])
astop = array("d", [0.00, 0.34, 0.61, 0.84, 1.00])
myPalette = []
fi = ROOT.TColor.CreateGradientColorTable(5, astop, ared, agreen, ablue, 100 )
for x in xrange(100): myPalette.append(fi+x)
ROOT.gStyle.SetPalette(100, array("i", myPalette))

import configureRun_cfi as runConfig

rootF = "DQM_V0001_R%09d__Global__CMSSW_X_Y_Z__RECO.root"%runConfig.RunNumber
run = runConfig.RunNumber
tDir = "DQMData/Run %d/MuonGEMRecHitsV/Run summary/GEMRecHitsTask"%(run)
oDir = "Run%06d_Plots/"%run

SLOTLIST=[]
VFATLIST=[] 
COLUMNLIST=[] 
ROWLIST=[]
LAYERLIST=[]
chamber=[]
columnStand=[]
rowStand=[]
layerSC=[]

runConfig.configureRun(SLOTLIST, VFATLIST, COLUMNLIST, ROWLIST, LAYERLIST, chamber, columnStand, rowStand, layerSC)

chamList = {}
for i, cn in enumerate(chamber):
  chamList["{}_{}_{}".format(COLUMNLIST[i*24], ROWLIST[i*24], LAYERLIST[i*24])] = cn

chIndex = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29]
chi = 0
chamGEO = {}
for c in [1,2,3]:
  for r in [1,2,3,4,5]:
    chamGEO["chamber_{}".format(chIndex[chi])] = [c,r]
    chi+=1
def findName(name):
  if not name.startswith("chamber") : return name
  tname = name.split("_")
  layer = int(tname[3])
  name = tname[0]+"_"+tname[1]
  colum = chamGEO[name][0]
  row = chamGEO[name][1]
  return chamList["{}_{}_{}".format(colum, row, layer)]  

def makeSummary():
  head = """
\documentclass{beamer}
\usefonttheme[onlylarge]{structurebold}
\setbeamerfont*{frametitle}{size=\\normalsize,series=\\bfseries}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{itemize items}[circle]
\setbeamerfont{page number in head/foot}{size=\small}
\setbeamertemplate{footline}[frame number]


\usepackage[english]{babel}
\usepackage[latin1]{inputenc}
\usepackage{times}
\usepackage{graphicx}

\usepackage[T1]{fontenc}
\usepackage{alltt}
%\usepackage{tikz}
%\usetikzlibrary{arrows}
%\\tikzstyle{block}=[draw opacity=0.7,line width=1.4cm]


\\newcommand{\\baseLoc}{./}

\\newcommand{\imageOne}[1]{
\scalebox{0.3}{
\includegraphics{\\baseLoc#1}
}
}

\\newcommand{\imageSix}[6]{
\scalebox{0.18}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
}
\\
\scalebox{0.18}{
\includegraphics{\\baseLoc#4}
\includegraphics{\\baseLoc#5}
\includegraphics{\\baseLoc#6}
}
}

\\begin{document}
"""
  tmp = """
\\begin{frame}[plain]{%s}
\imageSix{%s_gemDigi.png}{%s_recHit.png}{%s_trxy_eff.png}{%s_recHit_size.png}{%s_recHit_size_map.png}{%s_recHit_efficiency.png}
\end{frame}
"""
  tmp2 = """
\\begin{frame}[plain]{%s vs Det\_N\_LocalX for all N != %s}
\imageOne{%s_local_x.png}
\end{frame}

"""
  tmp3 = """
\\begin{frame}[plain]{Run Info.}
\\begin{itemize}
  \item RAWFileName : %s
  \item OutPutFile : %s
  \item MaxEvents : %d
  \item maxClusterSize : %d
  \item maxResidual : %1.2f cm
\end{itemize}
\end{frame}
"""
  import os
  os.chdir(oDir)
  outF = open(runConfig.OutputFileName.replace(".root", ".tex"), "w")
  outF.write(head)
  outF.write(tmp3%(runConfig.RAWFileName.split("/")[-1].replace("_","\_"), runConfig.OutputFileName.replace("_","\_"), runConfig.MaxEvents, runConfig.maxClusterSize, runConfig.maxResidual))
  for x in chamber:
    t = x.replace("GE1/1", "GE11")
    x = t+"/"+t
    outF.write(tmp%(t,x,x,x,x,x,x))
    if runConfig.makeTrack : outF.write(tmp2%(t,t,x))
  outF.write("\end{document}")
  outF.close() 
  os.system("latex --output-format=pdf "+runConfig.OutputFileName.replace(".root", ".tex"))
  import shutil
  shutil.copy2("./"+runConfig.OutputFileName.replace(".root", ".pdf"), "..")
def flipHist(hist):
  fhist = hist.Clone()
  fhist.Reset()
  ny = fhist.GetNbinsY()
  nx = fhist.GetNbinsX()
  for y in xrange(ny):
    for x in xrange(nx):
      tmpV = hist.GetBinContent(x+1,y+1)
      fhist.SetBinContent(x+1,ny-y,tmpV)
      fhist.GetYaxis().SetBinLabel(y+1,"{}".format(ny-y))
  return fhist
 
def setAxiNum(hist,axis,r, offSet=0):
  for x in xrange(r[0],r[1]+1):
    if axis == "x" or axis == "X":
      hist.GetXaxis().SetBinLabel(x,"{}".format(x+offSet))
    elif axis == "y" or axis == "Y":
      hist.GetYaxis().SetBinLabel(x,"{}".format(x+offSet))

def makeMapHist(hist):
  h2 = TH2D(hist.GetName()+"_map", "RecHit size mean value", 3,1,4,8,1,9) 
  h2.SetXTitle("vfat number")
  h2.SetYTitle("roll number")
  ny = hist.GetNbinsY()
  nx = hist.GetNbinsX()
  for y in xrange(ny):
    ent = 0
    val = 0
    for x in xrange(nx):
      tmpV =  hist.GetBinContent(x+1,y+1)
      ent += tmpV
      val += tmpV*(x+1)
    if ent == 0 : mean = 0
    else : mean = val/ent
    h2.SetBinContent(divmod(y,8)[0]+1, 8-divmod(y,8)[1], mean)
  return h2

def localXFitter(hist):
  #myFun = TF2("myfun", "[0]*x + [1] - y")
  myFun = TF1("myfun", "[0]*x + [1]")
  myFun.SetParameter(0,1)
  myFun.SetParameter(1,0)
  
  if hist.GetEntries() == 0: return 0,0
  hist.Fit("myfun")
  fitresult = TVirtualFitter.GetFitter()
  m = fitresult.GetParameter(0)
  b = fitresult.GetParameter(1)
  return m,b 
 
import optparse

def draw_occ(target_dir, h, ext =".png", opt = "colz"):
  gStyle.SetStatStyle(0)
  gStyle.SetOptStat(1110)
  name = findName(h.GetName())
  tname = h.GetName().split("_")
  etc = "_"
  for x in tname[4:]:
    etc += x+"_"
  title = name+" "+h.GetTitle()
  h.SetTitle(title)
  h.SetName(name.replace("GE1/1", "GE11")+etc[:-1])
  c = TCanvas(h.GetTitle(),h.GetName(),600,600)
  c_title = c.GetTitle()
  c.Clear()
  c.SetRightMargin(0.35)
  if not h:
    sys.exit('h does not exist')
  h.SetLineWidth(2)
  h.SetLineColor(kBlue)
  h.Draw(opt)
  c.SaveAs(target_dir + name.replace("GE1/1", "GE11")+"/"+ c_title + ext)



def draw_plot( file, tDir,oDir ) :
  c = TCanvas("c","c",600,600)
  dqm_file = TFile( file)
  d1 = dqm_file.Get(tDir)
  key_list =[]

  try :
    tlist = d1.GetListOfKeys()
  except :
    print oDir
    if ( oDir.find("Digi") != -1 ):
      tDir = "DQMData/Run 1/MuonGEMDigisV/Run summary/GEMDigiTask"
      d1 = dqm_file.Get(tDir)
      tlist = d1.GetListOfKeys()
    elif ( oDir.find("RecHit") != -1 ):
      tDir = "DQMData/Run 1/MuonGEMRecHitsV/Run summary/GEMRecHitTask"
      d1 = dqm_file.Get(tDir)
      tlist = d1.GetListOfKeys()
    else :
      print "error"
      exit(-1)
  for x in tlist :
    key_list.append(x.GetName())

  for hist in key_list :
    if ( hist.startswith("chamber") and hist.endswith("recHit_efficiency")):
      tmph = d1.Get(hist)
      thEff = d1.Get(hist.replace("recHit_efficiency", "th2D_eff"))
      tmph.Divide(thEff)
      tmph.SetXTitle("vfat number")
      tmph.SetYTitle("roll number")
      setAxiNum(tmph,"x",[1,3])
      setAxiNum(tmph,"y",[1,8])
      tmpf.GetZaxis().SetRange(0,1)
      tmpf = flipHist(tmph)
      draw_occ(oDir, tmpf, ".png", "colz text")

    if ( hist.startswith("chamber") and hist.endswith("trxy_eff")):
      tmph = d1.Get(hist)
      thEff = d1.Get(hist.replace("trxy_eff", "thxy_eff"))
      tmph.Divide(thEff)
      tmph.SetXTitle("x [cm]")
      tmph.SetYTitle("roll number")
      #setAxiNum(tmph,"x",[1,3])
      setAxiNum(tmph,"y",[1,8])
      tmpf = flipHist(tmph)
      tmpf.GetZaxis().SetRange(0,1)
      draw_occ(oDir, tmpf)
   
    elif ( hist.startswith("chamber") and hist.endswith("gemDigi")):
      tmph = d1.Get(hist)
      tmph.SetXTitle("Strip")
      tmph.SetYTitle("Roll Number (iEta)") 
      tmpf = flipHist(tmph)
      draw_occ(oDir, tmpf)
    elif ( hist.startswith("chamber") and hist.endswith("recHit")):
      tmph = d1.Get(hist)
      tmph.SetXTitle("x [cm]")
      tmph.SetYTitle("Roll Number (iEta)") 
      tmpf = flipHist(tmph)
      draw_occ(oDir, tmpf)
    elif ( hist.startswith("chamber") and hist.endswith("recHit_size")):
      tmph = d1.Get(hist)
      h2 = makeMapHist(tmph)
      tmph.SetXTitle("recHit size")
      tmph.SetYTitle("vfat number")
      setAxiNum(tmph,"y",[1,24],-1)
      draw_occ(oDir, tmph)
      setAxiNum(h2,"x",[1,3])
      tmpf = flipHist(h2)
      draw_occ(oDir, tmpf, ".png", "colz text")
    elif ( hist.startswith("chamber") and hist.endswith("residual")):
      tmph = d1.Get(hist)
      #tmph.SetAxisRange(-1,1,"X")
      #tmph.SetAxisRange(-1,1,"Y")
      tmph.SetXTitle("residual x [cm]")
      tmph.SetYTitle("residual y [cm]")
      draw_occ(oDir, tmph)
    elif ( hist.startswith("chamber") and hist.endswith("residual_r")):
      tmph = d1.Get(hist)
      tmph.SetXTitle("residual x [cm]")
      draw_occ(oDir, tmph)
    elif ( hist == "cluster_size"):
      tmph = d1.Get(hist)
      tmph.SetXTitle("cluster size")
      tmph.SetYTitle("count")
      draw_occ(oDir, tmph)
    elif (hist.startswith("chamber") and hist.endswith("local_x")):
      if not runConfig.makeTrack : continue
      c = TCanvas("local_X","local_x",600,600)
      tmph = d1.Get(hist)
      tmph.Draw()
      fitR = localXFitter(tmph) 
      tmph.SetXTitle("Det_1_LocalX [cm]")     
      tmph.SetYTitle("Det_N_LocalX for all N != 1 [cm]") 
      gStyle.SetStatStyle(0)
      gStyle.SetOptStat(1110)
      gStyle.SetStatStyle(0)
      gStyle.SetOptStat(1110)
      name = findName(hist)
      tname = hist.split("_")
      etc = "_"
      for x in tname[4:]:
        etc += x+"_"
      title = name+" "+tmph.GetTitle()
      tmph.SetTitle(title)
      tmph.SetName(name.replace("GE1/1", "GE11")+etc[:-1])
      #c.SetRightMargin(0.35)
      extraText = TLatex()
      extraText.SetNDC()
      extraText.SetTextFont(52)
      extraText.SetTextSize(0.03)
      extraText.DrawLatex(0.1,0.9,"fit result : y = mx + b (m = %1.2f, b = %1.2f)"%(fitR[0], fitR[1]))
      dName = name.replace("GE1/1", "GE11")+"/"
      c.SaveAs(oDir+dName+tmph.GetName()+".png")   

    else : continue
    #  draw_occ( oDir, d1.Get(hist) )

if __name__ == '__main__' :
  
 
  os.system("mkdir -p "+oDir )
  for c in chamber:
    os.system("mkdir -p "+oDir+"/"+c.replace("/",""))
  draw_plot(rootF,tDir,oDir)  
   
  makeSummary()
