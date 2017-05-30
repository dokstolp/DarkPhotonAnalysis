from ROOT import TFile, TTree

#Collects dataset trees and provides normalization information for each dataset.


#Collect root files for each dataset
pwd = "/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/RootFiles/"

darkFiles = {}
jetFiles = {}
treefin = {}
treepj = {}
treeTrig = {}

darkFactors = ["df1en0","df1en1","df1en2","df1en3"]
photonJets = ["GJets_HT100To200","GJets_HT200To400","GJets_HT400To600","GJets_HT600ToInf"]
diJet = ["QCD_Pt120to170","QCD_Pt170to300","QCD_Pt300to470","QCD_Pt470to600","QCD_Pt600to800","QCD_Pt800to1000","QCD_Pt1000to1400","QCD_Pt1400to1800","QCD_Pt1800to2400"]#,"QCD_Pt2400to3200","QCD_Pt3200toInf"]
diPho = ['DiPhoton']

for fi in darkFactors:
        darkFiles[fi] = TFile(pwd+fi+".root","OPEN")
        treepj[fi] = darkFiles[fi].Get("darkPho")
        treefin[fi] = darkFiles[fi].Get("darkTypes")
        treeTrig[fi] = darkFiles[fi].Get("Trigger")
for fij in photonJets:
        jetFiles[fij] = TFile(pwd+fij+".root","OPEN")
        treepj[fij] = jetFiles[fij].Get("darkPho")
        treeTrig[fij] = jetFiles[fij].Get("Trigger")
for fij in diJet:
        jetFiles[fij] = TFile(pwd+fij+".root","OPEN")
        treepj[fij] = jetFiles[fij].Get("darkPho")
        treeTrig[fij] = jetFiles[fij].Get("Trigger")
jetFiles['DiPhoton'] = TFile(pwd+"DiPhoton.root","OPEN")
treepj['DiPhoton'] = jetFiles['DiPhoton'].Get("darkPho")
treeTrig['DiPhoton'] = jetFiles['DiPhoton'].Get("Trigger")
dataFile = TFile(pwd+"Data.root","OPEN")
treepj['Data'] = dataFile.Get("darkPho")
treeTrig['Data'] = dataFile.Get("Trigger")

lumi = 23.22

PtBins = [145.,160.,190.,250.,400.,700.,1000]

ListBacks = photonJets+diJet+diPho
List = darkFactors+photonJets+diJet+diPho

#############################################################################################

#Signal Samples
colors = {}
colors['df1en0'] =4
colors['df1en1'] =6
colors['df1en2'] =16
colors['df1en3'] =2

darks = {}
darks['df1en0'] =1E0
darks['df1en1'] =1E-1
darks['df1en2'] =1E-2
darks['df1en3'] =1E-3

Nevents = {}
Nevents['df1en0'] = 9885.
Nevents['df1en1'] = 48650.
Nevents['df1en2'] = 188100.
Nevents['df1en3'] = 253497.

crossx = {}
crossx['df1en0']  =2.692E2*1E0
crossx['df1en1']  =2.692E2*1E-1
crossx['df1en2']  =2.692E2*1E-2
crossx['df1en3']  =2.692E2*1E-3

kFact = {}
kFact['df1en0'] = 1
kFact['df1en1'] = 1
kFact['df1en2'] = 1
kFact['df1en3'] = 1

#Photon + Jet background
Nevents['GJets_HT40To100']  = 4467939.
Nevents['GJets_HT100To200'] = 5131808.
Nevents['GJets_HT200To400'] = 10036339.
Nevents['GJets_HT400To600'] = 2529663.
Nevents['GJets_HT600ToInf'] = 2463751.

crossx['GJets_HT40to100']   = 23080E3
crossx['GJets_HT100To200']  = 9110E3
crossx['GJets_HT200To400']  = 2281E3
crossx['GJets_HT400To600']  = 273E3
crossx['GJets_HT600ToInf']  = 94.5E3

kFact['GJets_HT40To100']  = 1.0
kFact['GJets_HT100To200']  = 1.0
kFact['GJets_HT200To400']  = 1.0
kFact['GJets_HT400To600']  = 1.0
kFact['GJets_HT600ToInf'] = 1.0

#QCD background
Nevents['QCD_Pt120to170'] = 6708462.
Nevents['QCD_Pt170to300'] = 4150323.
Nevents['QCD_Pt300to470'] = 4150323.
Nevents['QCD_Pt470to600'] = 3959379.
Nevents['QCD_Pt600to800'] = 3895383.
Nevents['QCD_Pt800to1000'] = 3990472.
Nevents['QCD_Pt1000to1400'] = 2997099.
Nevents['QCD_Pt1400to1800'] = 396010.
Nevents['QCD_Pt1800to2400'] = 397083.
Nevents['QCD_Pt2400to3200'] = 399226.
Nevents['QCD_Pt3200toInf'] = 391904.

crossx['QCD_Pt120to170']  = 471100E3
crossx['QCD_Pt170to300']  = 117276E3
crossx['QCD_Pt300to470']  = 7823E3
crossx['QCD_Pt470to600']  = 648.2E3
crossx['QCD_Pt600to800']  = 186.9E3
crossx['QCD_Pt800to1000'] = 32.293E3
crossx['QCD_Pt1000to1400'] = 9.418E3
crossx['QCD_Pt1400to1800'] = 0.84E3
crossx['QCD_Pt1800to2400'] = 0.11E3
crossx['QCD_Pt2400to3200'] = 0.0068E3
crossx['QCD_Pt3200toInf'] = 0.00016E3

kFact['QCD_Pt120to170']  = 1.0
kFact['QCD_Pt170to300']  = 1.0
kFact['QCD_Pt300to470']  = 1.0
kFact['QCD_Pt470to600']  = 1.0
kFact['QCD_Pt600to800']  = 1.0
kFact['QCD_Pt800to1000']  = 1.0
kFact['QCD_Pt1000to1400'] = 1.0
kFact['QCD_Pt1400to1800'] = 1.0
kFact['QCD_Pt1800to2400'] = 1.0
kFact['QCD_Pt2400to3200'] = 1.0
kFact['QCD_Pt3200toInf'] = 1.0

#DiPhoton background
Nevents['DiPhoton'] = 9885.
crossx['DiPhoton'] = 139.3
kFact['DiPhoton'] = 1.0
