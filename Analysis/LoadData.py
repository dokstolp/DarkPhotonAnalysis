from ROOT import TFile, TTree

darkFiles = {}
jetFiles = {}
treefin = {}
treepj = {}
treeTrig = {}

darkFactors = ["df1en0","df1en1","df1en2","df1en3"]
#photonJets = ["pj_40to100","pj_100to200","pj_200to400","pj_400to600","pj_600toInf"]#,"photon_jet_1400to1800","photon_jet_1800"]
photonJets = ["pj_100to200","pj_200to400","pj_400to600","pj_600toInf"]#,"photon_jet_1400to1800","photon_jet_1800"]
diJet = ["qcd_120to170","qcd_170to300","qcd_300to470","qcd_470to600","qcd_600to800","qcd_800to1000","qcd_1000to1400","qcd_1400to1800","qcd_1800to2400","qcd_2400to3200","qcd_3200toInf"]
diPho = ['dipho']


for fi in darkFactors:
        darkFiles[fi] = TFile("../anaRoots/"+fi+".root","OPEN")
#        darkFiles[fi] = TFile("../withBDT/"+fi+".root","OPEN")
        treepj[fi] = darkFiles[fi].Get("darkPho")
        treefin[fi] = darkFiles[fi].Get("darkTypes")
for fij in photonJets:
        jetFiles[fij] = TFile("../anaRoots/"+fij+".root","OPEN")
#        jetFiles[fij] = TFile("../withBDT/"+fij+".root","OPEN")
        treepj[fij] = jetFiles[fij].Get("darkPho")
        treeTrig[fij] = jetFiles[fij].Get("darkTrig")
for fij in diJet:
        jetFiles[fij] = TFile("../anaRoots/"+fij+".root","OPEN")
#        jetFiles[fij] = TFile("../withBDT/"+fij+".root","OPEN")
        treepj[fij] = jetFiles[fij].Get("darkPho")
        treeTrig[fij] = jetFiles[fij].Get("darkTrig")
jetFiles['dipho'] = TFile("../anaRoots/dipho.root","OPEN")
#jetFiles['dipho'] = TFile("../withBDT/dipho.root","OPEN")
treepj['dipho'] = jetFiles['dipho'].Get("darkPho")
treeTrig['dipho'] = jetFiles['dipho'].Get("darkTrig")

dataFile = TFile("../anaRoots/data.root","OPEN")
#dataFile = TFile("../withBDT/data.root","OPEN")
treepj['data'] = dataFile.Get("darkPho")
treeTrig['data'] = dataFile.Get("darkTrig")

lumi = 12.9
#lumi = 26.3

PtBins = [145.,160.,190.,250.,400.,700.,1000]

ListBacks = photonJets+diJet+diPho
List = darkFactors+photonJets+diJet+diPho
#List = photonJets+diJet
updown = ['down']
#updown = ['up','down']

#############################################################################################


colors = {}
colors['df1en0'] =4
colors['df1en1'] =6
colors['df1en2'] =16
colors['df1en3'] =2
#colors['df9en1'] =28
#colors['df7en1'] =3
#colors['df5en1'] =46
#colors['df3en1'] =5
#colors['df8en2'] =7
#colors['df6en2'] =8
#colors['df4en2'] =9
#colors['df2en2'] =38
#colors['df5en3'] =30
#colors['df8en3'] =41

darks = {}
darks['df1en0'] =1.0
darks['df1en1'] =1E-1
darks['df1en2'] =1E-2
darks['df1en3'] =1E-3
#darks['df9en1'] =9E-1
#darks['df7en1'] =7E-1
#darks['df5en1'] =5E-1
#darks['df3en1'] =3E-1
#darks['df8en2'] =8E-2
#darks['df6en2'] =6E-2
#darks['df4en2'] =4E-2
#darks['df2en2'] =2E-2
#darks['df5en3'] =5E-3
#darks['df8en3'] =8E-3

Nevents = {}
#Nevents['df1en0'] =4026.0
#Nevents['df1en1'] =22023.0
#Nevents['df1en2'] =44597.0
#Nevents['df1en3'] =89088.0

Nevents['df1en0'] =9880.0
Nevents['df1en1'] =49899.0
Nevents['df1en2'] =199595.0
Nevents['df1en3'] =299394.0
#Nevents['df1en0'] =9980.0
#Nevents['df1en1'] =49699.0
#Nevents['df1en2'] =199195.0
#Nevents['df1en3'] =299394.0
#Nevents['df2en2'] =1000.0
#Nevents['df3en1'] =1000.0
#Nevents['df4en2'] =1000.0
#Nevents['df5en1'] =1000.0
#Nevents['df5en3'] =1000.0
#Nevents['df6en2'] =1000.0
#Nevents['df7en1'] =1000.0
#Nevents['df8en2'] =1000.0
#Nevents['df8en3'] =1000.0
#Nevents['df9en1'] =1000.0

crossx = {}
crossx['df1en0']  =6.918E2*1E0
crossx['df1en1']  =6.918E2*1E-1
crossx['df1en2']  =6.918E2*1E-2
crossx['df1en3']  =6.918E2*1E-3
#crossx['df2en2'] =2.52E2*2E-2
#crossx['df3en1'] =2.52E2*3E-1
#crossx['df4en2'] =2.52E2*4E-2
#crossx['df5en1'] =2.52E2*5E-1
#crossx['df5en3'] =2.52E2*5E-3
#crossx['df6en2'] =2.52E2*6E-2
#crossx['df7en1'] =2.52E2*7E-1
#crossx['df8en2'] =2.52E2*8E-2
#crossx['df8en3'] =2.52E2*8E-3
#crossx['df9en1'] =2.52E2*9E-1

kFact = {}
kFact['df1en0'] = 1
kFact['df1en1'] = 1
kFact['df1en2'] = 1
kFact['df1en3'] = 1
#kFact['df2en2'] = 1
#kFact['df3en1'] = 1
#kFact['df4en2'] = 1
#kFact['df5en1'] = 1
#kFact['df5en3'] = 1
#kFact['df6en2'] = 1
#kFact['df7en1'] = 1
#kFact['df8en2'] = 1
#kFact['df8en3'] = 1
#kFact['df9en1'] = 1


#Nevents['pj_40to100']  = 3506278
#Nevents['pj_100to200'] = 5116711
#Nevents['pj_200to400'] = 10462651
#Nevents['pj_400to600'] = 2507554
#Nevents['pj_600toInf'] = 2456253

#Nevents['qcd_120to170']  = 4797996
#Nevents['qcd_170to300']  = 6913886
#Nevents['qcd_300to470']  = 5968960
#Nevents['qcd_470to600']  = 3977770
#Nevents['qcd_600to800']  = 3979884
#Nevents['qcd_800to1000'] = 3973224


#Nevents['pj_40to100']  = 4458827
Nevents['pj_40to100']  = 4468703
#Nevents['pj_100to200'] = 5132757
Nevents['pj_100to200'] = 5142769
#Nevents['pj_200to400'] = 10312670
Nevents['pj_200to400'] = 10253085
#Nevents['pj_400to600'] = 2518528
Nevents['pj_400to600'] = 2528405
#Nevents['pj_600toInf'] = 2449258
Nevents['pj_600toInf'] = 2459246


crossx['pj_40to100']   = 23080E3
crossx['pj_100to200']  = 9110E3
crossx['pj_200to400']  = 2281E3
crossx['pj_400to600']  = 273E3
crossx['pj_600toInf']  = 94.5E3

kFact['pj_40to100']  = 1.0
kFact['pj_100to200']  = 1.0
kFact['pj_200to400']  = 1.0
kFact['pj_400to600']  = 1.0
kFact['pj_600toInf'] = 1.0


#Nevents['qcd_120to170']  = 6853860
#Nevents['qcd_170to300']  = 6904210
#Nevents['qcd_300to470']  = 575316
#Nevents['qcd_470to600']  = 3977770
#Nevents['qcd_600to800']  = 3979884
#Nevents['qcd_800to1000'] = 3973224
#Nevents['qcd_1000to1400'] = 9.418E3
#Nevents['qcd_1400to1800'] = 0.84E3
#Nevents['qcd_1800to2400'] = 0.11E3
#Nevents['qcd_2400to3200'] = 0.0068E3
#Nevents['qcd_3200toInf'] = 0.00016E3

Nevents['qcd_120to170']  = 6863805
Nevents['qcd_170to300']  = 6914063
Nevents['qcd_300to470']  = 5762821
Nevents['qcd_470to600']  = 3770179
Nevents['qcd_600to800']  = 3959746
Nevents['qcd_800to1000'] = 3976116 
Nevents['qcd_1000to1400'] = 2999055
Nevents['qcd_1400to1800'] = 396407
Nevents['qcd_1800to2400'] = 396093
Nevents['qcd_2400to3200'] = 399226
Nevents['qcd_3200toInf'] = 391904


crossx['qcd_120to170']  = 471100E3
crossx['qcd_170to300']  = 117276E3
crossx['qcd_300to470']  = 7823E3
crossx['qcd_470to600']  = 648.2E3
crossx['qcd_600to800']  = 186.9E3
crossx['qcd_800to1000'] = 32.293E3
crossx['qcd_1000to1400'] = 9.418E3
crossx['qcd_1400to1800'] = 0.84E3
crossx['qcd_1800to2400'] = 0.11E3
crossx['qcd_2400to3200'] = 0.0068E3
crossx['qcd_3200toInf'] = 0.00016E3

kFact['qcd_120to170']  = 1.0
kFact['qcd_170to300']  = 1.0
kFact['qcd_300to470']  = 1.0
kFact['qcd_470to600']  = 1.0
kFact['qcd_600to800']  = 1.0
kFact['qcd_800to1000']  = 1.0
kFact['qcd_1000to1400'] = 1.0
kFact['qcd_1400to1800'] = 1.0
kFact['qcd_1800to2400'] = 1.0
kFact['qcd_2400to3200'] = 1.0
kFact['qcd_3200toInf'] = 1.0

#Nevents['dipho'] = 1E4
#crossx['dipho'] = 123.5 
#kFact['dipho'] = 1.0
Nevents['dipho'] = 2036060
#Nevents['dipho'] = 3687041
crossx['dipho'] = 135100
kFact['dipho'] = 1.0
#kFact['dipho'] = 0.2
