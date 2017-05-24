from ROOT import TCanvas, TFile, TH1F, TF1, TPad, TLegend, TGraph, TAxis, TMultiGraph, THStack, TColor, TH1D, TLatex, gROOT, gBenchmark, gRandom, gSystem, Double, gPad, gStyle, TGraphErrors, TH2F, TH2D, TPie
from LoadData import *
#from plot_fun import *
from plot_variables import *
from array import array
import scipy

vari = [Jet_NJets,Pho_phi,Pho_eta,Jet_phi,Jet_eta,Pho_pt,Jet_pt]
vari = [Jet_NchargedHad1,Jet_NchargedHad5,Jet_NchargedHad7,Jet_NchargedHad10,Jet_NneutralHad1,Jet_NneutralHad5,Jet_NneutralHad7,Jet_NneutralHad10,Jet_Ngamma1,Jet_Ngamma5,Jet_Ngamma7,Jet_Ngamma10]
vari = [Jet_CEF,Jet_CHF,Jet_NConst,Jet_n60,Jet_n90,Jet_ptDist,Jet_etaWidthInECal,Jet_phiWidthInECal,Jet_etaWidthInHCal,Jet_phiWidthInHCal]
cuts = {'noCuts-DiScale':"(neg_weight*event_weight*(isOverlap==0))"}

def trigger(nbins,low,high):
	dataFile = TFile("../anaRoots/data.root","OPEN")
	treeTrig = dataFile.Get("Trigger")
	cutter = cuts[region]
	can = TCanvas("can","can",1000,900)
	can.SetBottomMargin(0.3)
	can.SetRightMargin(0.06)
	gStyle.SetOptStat(0)
	histNumerator = TH1F('histNumerator',";Photon p_T [GeV];Efficiency",nbins,low,high)
	treeTrig.Draw('TrigPhoton_pt>>histNumerator','(isHLT165==1) && (isPrescaled==1)')
	histDenominator = TH1F('histDenominator',";Photon p_T [GeV];Efficiency",nbins,low,high)
	treeTrig.Draw('TrigPhoton_pt>>histDenominator','(isPrescaled==1)')
	histNumerator.Divide(histDenominator)
	histNumerator.Draw('e1')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	can.SaveAs("plots/trigger.pdf")
	can.Close()

def piePlot():
	darknames = {'df1en0':'f_D = 1E0','df1en1':'f_D = 1E-1','df1en2':'f_D = 1E-2','df1en3':'f_D = 1E-3'}
	colors = [2,4,6]
	labels = ['Photon','Jet','MET']
	can = TCanvas('can','can')
	can.Divide(1,4)
	iter = 1
	pie = {}
	for tre in darkFactors:
		can.cd(iter)
		GenDarkPhoisPho = 0
		GenDarkPhoisJet = 0
		GenDarkPhoisMet = 0
		for event in treefin[tre]:
			GenDarkPhoisPho+=event.GenDarkPho_isPho
			GenDarkPhoisJet+=event.GenDarkPho_isJet
			GenDarkPhoisMet+=event.GenDarkPho_isMet
		vals = [GenDarkPhoisPho,GenDarkPhoisJet,GenDarkPhoisMet]
		pie[tre] = TPie('pie',darknames[tre],3,scipy.array(vals),scipy.array(colors))
		pie[tre].SetLabels(labels)
		pie[tre].Draw("3d")
		iter+=1
	can.SaveAs('plots/pieChart.pdf')	

def interactionLocation():
        can3 = {}
        for tre in darkFactors:
                can = TCanvas(tre,"can",1000,900)
                histo = TH2F('histo',"f_D = "+str(darks[tre])+";#gamma_D Interaction Z [cm];#gamma_D Interaction r [cm]",30,0,600,30,0,300)
                treepj[tre].Draw("darkRadius:fabs(darkZ)>>histo",weight)
                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
#                scale = 100.0/Variables[tre].Integral(0,xnbin+1,0,ynbin+1)
                Variables[tre].Scale(scale)
                Variables[tre].Draw("COLZ")
                can.SaveAs("plots/Signal_"+tre+"-2D.pdf")


def plot_2D():
        can = TCanvas("can","can", 1000,900)
        diphoHist = TH2F('diphoHist',";#Delta R(#gamma_{Gen},Jet);Number of Jet Constituents",20,0,4,20,0,60)
        treepj['dipho'].Draw("Jet_NConst:JetGenPho_dR>>diphoHist",'(neg_weight*event_weight)')
        diphoHist.Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
	diphoHist.Draw()
        can.SaveAs("plots/DiPhoton2D-DR_JetGenPho-Jet_NConst.pdf")

def dipho_NConst(nbins,low,high,isCorr):
	diphoCorr = ['',' && (JetPho_dR<0.1)']
	cutter = cuts[region]+diphoCorr[isCorr]
	can = TCanvas("can","can",1000,900)
	can.SetBottomMargin(0.3)
	can.SetRightMargin(0.06)
	gPad.SetLogy(1)
	gStyle.SetOptStat(0)
	histDiPhoton = TH1F('histDiPhoton',";Number of Jet Constituents;Events",nbins,low,high)
	treepj['dipho'].Draw('Jet_NConst>>histDiPhoton',cutter)
	scalable = lumi*crossx['dipho']*kFact['dipho']/Nevents['dipho']
	histDiPhoton.Scale(scalable)
	histDiPhoton.SetLineColor(3)
	histSignal = TH1F('histSignal',";Number of Jet Constituents;Events",nbins,low,high)
	treepj['df1en0'].Draw('Jet_NConst>>histSignal',cutter)
	scalable = lumi*crossx['df1en0']*kFact['df1en0']/Nevents['df1en0']
	histSignal.Scale(scalable)
	histSignal.SetLineColor(4)
	histDiPhoton.Draw('e1')
	histSignal.Draw('same')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	led = TLegend(0.7, 0.7, 0.94, 0.9)
	led.AddEntry(Variables['dipho'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"f_D = 1E0")
	led.SetFillColor(0)
	led.Draw("same")
	can.SaveAs("plots/dipho_NConst-"+isCorr+".pdf")
	can.Close()

def pileup(isCorr):
	Corr = ['','(event_weight)']
	can = TCanvas("can","can", 1000,900)
	gStyle.SetOptStat(0)
	phojet = TH1D('b','b',60,0,60)
	phojet.SetMarkerColor(7)
	qcd = TH1D('c','c',60,0,60)
	qcd.SetMarkerColor(5)
	Variables = {}
	for tre in List:
		Variables[tre] = TH1F(tre,";Number of Vertices;Events",60,0,60)
#		if tre == "dipho":
#			weight = cutter+"&& (JetPho_dR<0.1)"
		treeTrig[tre].Draw("NVerts>>"+tre,Corr[isCorr])
		Variables[tre].Sumw2()
		scalable = lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Variables[tre].Scale(scalable)
		if tre.startswith("df"):
			Variables[tre].SetMarkerColor(colors[tre])
			Variables[tre].SetLineWidth(3)
			scalable = Variables[tre].Integral()
			Variables[tre].Scale(1/scalable)
		if tre.startswith("pj"):
			phojet.Add(Variables[tre])
		if tre.startswith("qcd"):
			qcd.Add(Variables[tre])
		if tre.startswith("dipho"):
			Variables[tre].SetMarkerColor(3)
			scalable = Variables[tre].Integral()
			Variables[tre].Scale(1/scalable)
	scalable = phojet.Integral()
	phojet.Scale(1/scalable)
	scalable = qcd.Integral()
	qcd.Scale(1/scalable)
	data = TH1D('data','',60,0,60)
	dataFi = TFile("/uscms_data/d3/dokstolp/darkphoton/13TeV/Analysis/CMSSW_8_0_18_patch1/src/Runner/Weights/forPileup/data.root","OPEN")
	data = dataFi.Get("hist")
	scalar = data.Integral()
	data.Scale(1/scalar)
	data.SetLineColor(1)
	data.SetLineWidth(2)
	data.Draw("e1")
	Variables['df1en0'].Draw("e1same")		
	Variables['df1en1'].Draw("e1same")		
	Variables['df1en2'].Draw("e1same")		
	Variables['df1en3'].Draw("e1same")		
	Variables['dipho'].Draw('e1same')
	qcd.Draw('e1same')
	phojet.Draw('e1same')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	led = TLegend(0.7, 0.7, 0.94, 0.9)
	led.AddEntry(phojet,"#gamma + Jet",'f')
	led.AddEntry(qcd,"QCD DiJet",'f')
	led.AddEntry(Variables['dipho'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"f_D = 1E0}")
	led.AddEntry(Variables['df1en1'],"f_D = 1E-1}")
	led.AddEntry(Variables['df1en2'],"f_D = 1E-2}")
	led.AddEntry(Variables['df1en3'],"f_D = 1E-3}")
	led.AddEntry(data,"Data")
	led.SetFillColor(0)
	led.Draw("same")
	can.SaveAs("plots/pileup-"+isCorr+".pdf")
	can.Close()
	phojet.Delete()
	qcd.Delete()


def aXe():
        can = TCanvas("can","can", 1000,900)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        VariablesNum = {}
        VariablesDen = {}
        treeG = {}
        for tre in darkFactors:
                treeG[tre] = darkFiles[tre].Get("darkGens")
                VariablesNum[tre] = TH1F(tre+"-num",";Dark Photon Pt [GeV];Acceptance x Efficiency",26,150,800)
                VariablesDen[tre] = TH1F(tre+"-den",";Dark Photon Pt [GeV];Acceptance x Efficiency",26,150,800)
                VariablesNum[tre].SetLineColor(colors[tre])
                VariablesNum[tre].SetLineWidth(2)
                treepj[tre].Draw("GenDarkPho_pt>>"+tre+"-num")
                treeG[tre].Draw("GenDarkPho_pt>>"+tre+"-den")
                VariablesNum[tre].Sumw2()
                VariablesDen[tre].Sumw2()
                VariablesNum[tre].Divide(VariablesDen[tre])
		VariablesNum[tre].SetMaximum(5.0)
        VariablesNum['df1en0'].Draw('e1')
        VariablesNum['df1en1'].Draw('e1same')
        VariablesNum['df1en2'].Draw('e1same')
        VariablesNum['df1en3'].Draw('e1same')
        led = TLegend(0.6, 0.6, 0.9, 0.9)
        led.AddEntry(VariablesNum['df1en0'],"f_D = 1E0")
        led.AddEntry(VariablesNum['df1en1'],"f_D = 1E-1}")
        led.AddEntry(VariablesNum['df1en2'],"f_D = 1E-2}")
        led.AddEntry(VariablesNum['df1en3'],"f_D = 1E-3}")
        led.SetFillColor(0)
        led.Draw("same")
        can.SaveAs("plots/aXe-"+region+".pdf")


def noPull(var,title,xtitle,nbin,low,high,region):
	cutter = cuts[region]
	weight = cutter
	can = TCanvas("can","can", 1000,900)
	gPad.SetLogy(1)
	gStyle.SetOptStat(0)
	backs = TH1D('a','a',nbin,low,high)
	backs.SetFillStyle(3344)
	backs.SetFillColor(1)
	phojet = TH1D('b','b',nbin,low,high)
	phojet.SetFillColor(7)
	phojet.SetLineColor(7)
	qcd = TH1D('c','c',nbin,low,high)
	qcd.SetFillColor(5)
	qcd.SetLineColor(5)
	stack = THStack('d',title+";"+xtitle+";"+"Events")
	Variables = {}
	print "Variable: "+var
	for tre in List:
		histName = var+tre
		Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
		if tre == "dipho":
			weight = cutter+"&& (JetPho_dR<0.1)"
		treepj[tre].Draw(var+">>"+histName,weight)
		Variables[tre].Sumw2()
		scalable = lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Variables[tre].Scale(scalable)
		if tre.startswith("df"):
			Variables[tre].SetLineColor(colors[tre])
			Variables[tre].SetLineWidth(3)
		if tre.startswith("pj"):
			backs.Add(Variables[tre])
			phojet.Add(Variables[tre])
		if tre.startswith("qcd"):
			backs.Add(Variables[tre])
			qcd.Add(Variables[tre])
		if tre.startswith("dipho"):
			Variables[tre].SetFillColor(3)
			Variables[tre].SetLineColor(3)
			backs.Add(Variables[tre])
	Variables['data'] = TH1F('data',";"+xtitle+";"+"Events",nbin,low,high)
	treepj['data'].Draw(var+">>data",cutter,"goff")
	Variables['data'].SetLineColor(1)
	Variables['data'].SetLineWidth(2)
	stack.Add(Variables['dipho'])
	stack.Add(qcd)
	stack.Add(phojet)		
	stack.SetMaximum(1000000.0)
	stack.SetMinimum(0.001)
	stack.Draw("hist")
	stack.GetXaxis().SetLabelSize(0)
	backs.Draw("e2same")
	Variables['data'].Draw("e1same")	
	Variables['df1en0'].Draw("histsame")		
	Variables['df1en1'].Draw("histsame")		
	Variables['df1en2'].Draw("histsame")		
	Variables['df1en3'].Draw("histsame")		
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	led = TLegend(0.7, 0.7, 0.94, 0.9)
	led.AddEntry(phojet,"#gamma + Jet",'f')
	led.AddEntry(qcd,"QCD DiJet",'f')
	led.AddEntry(Variables['dipho'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"f_D = 1E0}")
	led.AddEntry(Variables['df1en1'],"f_D = 1E-1}")
	led.AddEntry(Variables['df1en2'],"f_D = 1E-2}")
	led.AddEntry(Variables['df1en3'],"f_D = 1E-3}")
	led.AddEntry(Variables['data'],"Data")
	led.SetFillColor(0)
	led.Draw("same")
	can.SaveAs("plots/noPull"+region+"-"+var+".pdf")
	can.Close()
	phojet.Delete()
	qcd.Delete()

def thresholds(var,title,xtitle,nbin,low,high,region):
	cutter = cuts[region]
	weight = cutter
	can = TCanvas("can","can", 1000,900)
	can.SetBottomMargin(0.3)
	can.SetRightMargin(0.06)
	can.cd()
	gPad.SetLogy(1)
	gStyle.SetOptStat(0)
	backs = TH1D('a','a',nbin,low,high)
	backs.SetFillStyle(3344)
	backs.SetFillColor(1)
	phojet = TH1D('b','b',nbin,low,high)
	phojet.SetFillColor(7)
	phojet.SetLineColor(7)
	qcd = TH1D('c','c',nbin,low,high)
	qcd.SetFillColor(5)
	qcd.SetLineColor(5)
	stack = THStack('d',title+";"+xtitle+";"+"Events")
	Variables = {}
	print "Variable: "+var
	for tre in List:
		histName = var+tre
		Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
		if tre == "dipho":
			weight = cutter+"&& (JetPho_dR<0.1)"
		treepj[tre].Draw(var+">>"+histName,weight)
		Variables[tre].Sumw2()
		scalable = lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Variables[tre].Scale(scalable)
		if tre.startswith("df"):
			Variables[tre].SetLineColor(colors[tre])
			Variables[tre].SetLineWidth(3)
		if tre.startswith("pj"):
			backs.Add(Variables[tre])
			phojet.Add(Variables[tre])
		if tre.startswith("qcd"):
			backs.Add(Variables[tre])
			qcd.Add(Variables[tre])
		if tre.startswith("dipho"):
			Variables[tre].SetFillColor(3)
			Variables[tre].SetLineColor(3)
			backs.Add(Variables[tre])
	Variables['data'] = TH1F('data',";"+xtitle+";"+"Events",nbin,low,high)
	treepj['data'].Draw(var+">>data",cutter,"goff")
	Variables['data'].SetLineColor(1)
	Variables['data'].SetLineWidth(2)
	stack.Add(Variables['dipho'])
	stack.Add(qcd)
	stack.Add(phojet)		
	stack.SetMaximum(1000000.0)
	stack.SetMinimum(0.001)
	stack.Draw("hist")
	stack.GetXaxis().SetLabelSize(0)
	backs.Draw("e2same")
	Variables['data'].Draw("e1same")	
	Variables['df1en0'].Draw("histsame")		
	Variables['df1en1'].Draw("histsame")		
	Variables['df1en2'].Draw("histsame")		
	Variables['df1en3'].Draw("histsame")		
	pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
	pad.SetTopMargin(0.7)
	pad.SetFillColor(0)
	pad.SetGridy(1)
	pad.SetFillStyle(0)
	pad.Draw()
	pad.cd(0)
	pad.SetRightMargin(0.06)
	Pull  = TH1D('Pull', 'Pull', nbin, low, high)
	Pull = Variables['data'].Clone()
	Pull.Divide(backs)
	Pull.SetMarkerStyle(20)
	Pull.SetMaximum(2.0)
	Pull.SetMinimum(0.0)
	Pull.SetFillColor(2)
	Pull.GetXaxis().SetTitle(xtitle)
	Pull.GetYaxis().SetTitleSize(0.04)
	Pull.GetYaxis().SetTitle('Data/MC')
	Pull.SetMarkerSize(0.7)
	Pull.GetYaxis().SetNdivisions(5)
	line = TF1("line","1",low,high)
	line.SetLineColor(2)
	line.SetLineWidth(2)
	Pull.Draw("e1")
	line.Draw('same')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	led = TLegend(0.7, 0.7, 0.94, 0.9)
	led.AddEntry(phojet,"#gamma + Jet",'f')
	led.AddEntry(qcd,"QCD DiJet",'f')
	led.AddEntry(Variables['dipho'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"f_D = 1E0")
	led.AddEntry(Variables['df1en1'],"f_D = 1E-1")
	led.AddEntry(Variables['df1en2'],"f_D = 1E-2")
	led.AddEntry(Variables['df1en3'],"f_D = 1E-3")
	led.AddEntry(Variables['data'],"Data")
	led.SetFillColor(0)
	led.Draw("same")
	can.SaveAs("plots/"+region+"-"+var+".pdf")
	can.Close()
	phojet.Delete()
	qcd.Delete()
