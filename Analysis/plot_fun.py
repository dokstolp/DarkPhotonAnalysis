import ROOT, sys, os, string, re
from ROOT import TCanvas, TFile, TH1F, TF1, TPad, TLegend, TGraph, TAxis, TMultiGraph, THStack, TColor, TH1D, TLatex, gROOT, gBenchmark, gRandom, gSystem, Double, gPad, gStyle, TGraphErrors, TH2F, TH2D
from array import array
import math
from math import sqrt
import scipy
from LoadData import *

def plot_2D(xvar,yvar,title,xtitle,ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh,region,cutter):
	if cutter == "":
		weight = "event_weight"
	else:
		weight = "("+cutter+")*event_weight"
        can = TCanvas("can","can", 1000,900)
        can.cd()
#	gPad.SetLogy(1)
#	gPad.SetLogz(1)
        gStyle.SetOptStat(0)
        backs = TH2D('a',"DiPhoton;"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
	diphob = TH1F('b',"DiPhoton;"+xtitle+";",xnbin,xlow,xhigh)
        Variables = {}
        for tre in ListBacks:
                histName = xvar+yvar+tre
		if (yvar == "Jet_NConst - Jet_NCH"):
			histName = xvar+"Jet_NNH"+tre
                Variables[tre] = TH2F(histName,";"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
                treepj[tre].Draw(yvar+":"+xvar+">>"+histName,weight)
                Variables[tre].Sumw2()
                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
                if tre.startswith("pj"):
                        backs.Add(Variables[tre])
                if tre.startswith("qcd"):
                        backs.Add(Variables[tre])
                if tre.startswith("dipho"):
                        backs.Add(Variables[tre])
			treepj[tre].Draw(xvar+">>b",weight)
			diphob.Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
#			diphob.Add(Variables[tre])
	backs.Draw("COLZ")
	if (yvar == "Jet_NConst - Jet_NCH"):
		can.SaveAs("plots/Background"+region+"-"+xvar+"-Jet_NNH.pdf")
	else:
		can.SaveAs("plots/Background"+region+"-"+xvar+"-"+yvar+".pdf")

#	can2 = TCanvas("can2","can",1000,900)
        Variables['data'] = TH2D('data',"Data/Background;"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
        treepj['data'].Draw(yvar+":"+xvar+">>data",cutter)
#	Variables['data'].Draw("COLZ")
#	if (yvar == "Jet_NConst - Jet_NCH"):
#		can2.SaveAs("plots/Data"+region+"-"+xvar+"-Jet_NNH.pdf")
#	else:
#		can2.SaveAs("plots/Data"+region+"-"+xvar+"-"+yvar+".pdf")
#	can3 = TCanvas("can3","can",1000,900)
#	diphob.Draw()
#	can3.SaveAs("plots/Hist1D-"+region+"-"+xvar+".pdf")

#	can3 = {}
#	for tre in darkFactors:
#		can3[tre] = TCanvas(tre,"can",1000,900)
#                histName = xvar+yvar+tre
#		if (yvar == "Jet_NConst - Jet_NCH"):
#			histName = xvar+"Jet_NNH"+tre
#                Variables[tre] = TH2F(histName,";"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
#                treepj[tre].Draw(yvar+":"+xvar+">>"+histName,weight)
#                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
#		Variables[tre].Draw("COLZ")
#		if (yvar == "Jet_NConst - Jet_NCH"):
#			can3[tre].SaveAs("plots/Signal"+tre+"-"+region+"-"+xvar+"-Jet_NNH.pdf")
#		else:
#			can3[tre].SaveAs("plots/Signal"+tre+"-"+region+"-"+xvar+"-"+yvar+".pdf")

        Pull  = TH2D('Pull',title+";"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
        Pull = Variables['data'].Clone()
#        Pull.Add(backs,-1)
        Pull.Divide(backs)
        Pull.SetMaximum(2.0)
        Pull.SetMinimum(0.0)
        for i in range(xnbin):
                i += 1
		for j in range(ynbin):
			j += 1
			if Variables['data'].GetBinContent(i,j) != 0:
#				Pull.SetBinContent(i,j,Pull.GetBinContent(i,j)/Pull.GetBinError(i,j))
                        	Pull.SetBinContent(i,j,Pull.GetBinContent(i,j))
				print "x: "+str(i)+" y: "+str(j)+" Data: "+str(Variables['data'].GetBinContent(i,j))+" Backgrounds: "+str(backs.GetBinContent(i,j))+" Data/Backgrounds: "+str(Pull.GetBinContent(i,j))
			else: Pull.SetBinContent(i,j,0)
	can4 = TCanvas("can4","can",1000,900)
        Pull.Draw("COLZ")
	can4.SaveAs("plots/Pull2D-"+region+"-"+xvar+"-"+yvar+".pdf")

def normalized(var,title,xtitle,nbin,low,high,region,cutter):
        if cutter == "":
                weight = "event_weight"
        else:
                weight = cutter+"*event_weight"
        can = TCanvas("can","can", 1000,900)
	gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        Variables = {}
        print "Variable: "+var
        for tre in darkFactors:
                histName = var+tre
                Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
                treepj[tre].Draw(var+">>"+histName,weight)
#                Variables[tre].Sumw2()
                scalable = Variables[tre].Integral(0,nbin+1)
                Variables[tre].Scale(1/scalable)
                if tre.startswith("df"):
                        Variables[tre].SetLineColor(colors[tre])
                        Variables[tre].SetLineWidth(3)
        Variables['df1en0'].SetLineColor(4)
        Variables['df1en1'].SetLineColor(6)
        Variables['df1en2'].SetLineColor(16)
        Variables['df1en3'].SetLineColor(2)
        Variables['df1en0'].Draw("e1")
        Variables['df1en1'].Draw("e1same")
        Variables['df1en2'].Draw("e1same")
        Variables['df1en3'].Draw("e1same")
        led = TLegend(0.6, 0.6, 0.9, 0.9)
        led.AddEntry(Variables['df1en0'],"#alpha_{Dark} = #alpha_{EM} x 10^{0}")
        led.AddEntry(Variables['df1en1'],"#alpha_{Dark} = #alpha_{EM} x 10^{-1}")
        led.AddEntry(Variables['df1en2'],"#alpha_{Dark} = #alpha_{EM} x 10^{-2}")
        led.AddEntry(Variables['df1en3'],"#alpha_{Dark} = #alpha_{EM} x 10^{-3}")
        led.SetFillColor(0)
        led.Draw("same")
        print "df = 10^0       "+str(Variables["df1en0"].Integral(0,nbin+1))
        print "df = 10^-1      "+str(Variables["df1en1"].Integral(0,nbin+1))
        print "df = 10^-2      "+str(Variables["df1en2"].Integral(0,nbin+1))
        print "df = 10^-3      "+str(Variables["df1en3"].Integral(0,nbin+1))
        can.SaveAs("plots/norm-"+region+"-"+var+".pdf")
	
def aXe(region,cutter):
	weight = cutter
#        if cutter == "":
#                weight = "event_weight"
#        else:
#                weight = cutter+"*event_weight"
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
                treepj[tre].Draw("GenDarkPho_pt>>"+tre+"-num",weight)
                treeG[tre].Draw("GenDarkPho_pt>>"+tre+"-den",weight)
                VariablesNum[tre].Sumw2()
                VariablesDen[tre].Sumw2()
#                scalable = Variables[tre].Integral(0,nbin+1)
                VariablesNum[tre].Divide(VariablesDen[tre])
#                if tre.startswith("df"):
#			VariablesNum[tre].SetMarkerColor(colors[tre])
#        VariablesNum['df1en0'].SetLineColor(4)
#        Variables['df1en1'].SetLineColor(6)
#        Variables['df1en2'].SetLineColor(16)
#        Variables['df1en3'].SetLineColor(2)
        VariablesNum['df1en0'].Draw('e1')
        VariablesNum['df1en1'].Draw('e1same')
        VariablesNum['df1en2'].Draw('e1same')
        VariablesNum['df1en3'].Draw('e1same')
        led = TLegend(0.7, 0.75, 0.9, 0.9)
        led.AddEntry(VariablesNum['df1en0'],"#alpha_{Dark} = #alpha_{EM} x 10^{0}")
        led.AddEntry(VariablesNum['df1en1'],"#alpha_{Dark} = #alpha_{EM} x 10^{-1}")
        led.AddEntry(VariablesNum['df1en2'],"#alpha_{Dark} = #alpha_{EM} x 10^{-2}")
        led.AddEntry(VariablesNum['df1en3'],"#alpha_{Dark} = #alpha_{EM} x 10^{-3}")
        led.SetFillColor(0)
        led.Draw("same")
#        print "df = 10^0       "+str(VariablesNum["df1en0"].Integral(0,27))
#        print "df = 10^-1      "+str(VariablesNum["df1en1"].Integral(0,27))
#        print "df = 10^-2      "+str(VariablesNum["df1en2"].Integral(0,27))
#        print "df = 10^-3      "+str(VariablesNum["df1en3"].Integral(0,27))
        can.SaveAs("plots/aXe-"+region+".pdf")


def weighter():
	var = "Jet_NConst"
	title = ""
	xtitle = "Number of PF Constituents"
	nbin = 150
	low = 0
	high = 300
	cutter = {"Neutrality":"(nJet == 1 && (Jet_NCH > 16 || Jet_CHF > 0.2))","ShowerShape":"(nJet==1 && (Jet_NCH < 16 && Jet_CHF < 0.2) &&"+Chi2+">5)"}
	errors = {"Neutrality":[],"ShowerShape":[]}
	weighted = {"Neutrality":[],"ShowerShape":[]}
	gROOT.ProcessLine('.L CFunk.C')
	for region in cutter:
		weight = "eve_weight*"+cutter[region]
		backs = TH1D('a','a',nbin,low,high)
		phojet = TH1D('b','b',nbin,low,high)
		qcd = TH1D('c','c',nbin,low,high)
		Variables = {}
		print "Variable: "+var
		if region == "ShowerShape":
			weight = weight+"*CFunk(Jet_NConst)"
		for tre in List:
			print tre
			histName = var+tre
			Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
			tree[tre].Draw(var+">>"+histName,weight,"goff")
			Variables[tre].Sumw2()
			Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
			if tre.startswith("df"):
				Variables[tre].Scale(5.5)
			if tre.startswith("photon_jet"):
				backs.Add(Variables[tre])
				phojet.Add(Variables[tre])
			if tre.startswith("qcd"):
				backs.Add(Variables[tre])
				qcd.Add(Variables[tre])
		Variables['data'] = TH1F('data',";"+xtitle+";"+"Events",nbin,low,high)
		tree['data'].Draw(var+">>data",cutter[region],"goff")
#	tree['data'].Draw(var+">>data",blindcuts,"goff")
		weights_ = Variables['data'].Clone("weights_")
		den = backs.Clone("den")

		deltaH = weights_.Integral(0,nbin+1)
		if (abs(1.0-deltaH) > 0.02):
			weights_.Scale(1/deltaH)
			print weights_.Integral(0,nbin+1)
		deltaMC = den.Integral(0,nbin+1)
		if (abs(1.0-deltaMC) > 0.02):
			den.Scale(1/deltaMC)
			print den.Integral(0,nbin+1)
		print "Data:"+str(deltaH)+"\tMC:"+str(deltaMC)
		weights_.Sumw2()
		den.Sumw2()
		weights_.Divide(den)
		for i in range(low,high,(high-low)/nbin):
			ibin = weights_.GetXaxis().FindBin(i)
			weighted[region].append(weights_.GetBinContent(ibin))
			errors[region].append(weights_.GetBinError(ibin))
	lister = []
	for ibin in range(nbin):
#		val = weighted["ShowerShape"][ibin]
#		val = weighted["Neutrality"][ibin]
		error = errors["ShowerShape"][ibin]
		if errors["ShowerShape"][ibin] == 0 and errors["Neutrality"][ibin] != 0:
			val = weighted["Neutrality"][ibin]
		elif errors["Neutrality"][ibin] == 0 and errors["ShowerShape"][ibin] != 0:
			val = weighted["ShowerShape"][ibin]
		elif errors["Neutrality"][ibin] == 0 and errors["ShowerShape"][ibin] == 0:
			val = 0
		else:
			val = (weighted["ShowerShape"][ibin]/errors["ShowerShape"][ibin]+weighted["Neutrality"][ibin]/errors["Neutrality"][ibin])/(1/errors["ShowerShape"][ibin] + 1/errors["Neutrality"][ibin])
#		print str(ibin)+"\t"+str(val)
		lister.append(val)
		
#	weights_.Divide(den)
#	outfi = TFile("forShapeMatching/Neutrality.root","RECREATE")
#	weights_.Write()
#	outfi.Close()
	pre = """#include "TH1F.h"
#include "TGraph.h"

"""
#	xaxis = "double xaxis["+str(nbin+1)+"] = {0.0,"
#	yaxis = "double yaxis["+str(nbin+1)+"] = {0.0,"
	yaxis = "double yaxis["+str(nbin)+"] = {"
	for i in lister:
	#for i in range(low+1,high+1,(high-low)/nbin):
#			xaxis+=str(i)+','
		yaxis+=str(i)+','
#                print str(i-1)+"\t"+str(weights_.GetBinContent(i))+"\t"+str(den.GetBinContent(i))
#		yaxis+=str(weights_.GetBinContent(i)/den.GetBinContent(i))+','
	
#		yaxis+=str(weights_.GetBinContent(int(inter)))+','
#		xaxis+='};'
	yaxis+='};'


	mid = """double CFunk(double val){
	TH1F * hist = new TH1F("hist","","""+str(nbin)+","+str(low)+","+str(high)+""");
	for(int i=1;i<"""+str(nbin+1)+""";i++){
		hist->SetBinContent(i,yaxis[i-1]);
	}
	int ibin = hist->GetXaxis()->FindBin(val);
	double outval = hist->GetBinContent(ibin);
	hist->Delete();
	return outval;
}"""

#	mid = """double CFunk(double val){
#	TGraph * grap = new TGraph("""+str(nbin+1)+""",xaxis,yaxis);
#	double outval = grap->Eval(val);
#	grap->Delete();
#	return outval;
#}"""

	output_file = open('CFunk.C','w')
	code_line = pre+yaxis+"\n"+mid
#	code_line = pre+xaxis+"\n"+yaxis+"\n"+mid
	output_file.write(code_line)
#	os.system('root CFunk.C++')
#	fout = TFile("nconst.root","RECREATE")
#	weights_.Write()
#	fout.Close()

def optimizer(var,title,xtitle,nbin,low,high,region,cutter):
	print "*************"+region+"*****************"
	if cutter == "":
		weight = "event_weight"
	else:
		weight = "("+cutter+")*event_weight"
        can = TCanvas("can","can", 1000,700)
        backs = TH1D('a','a',nbin,low,high)
        phojet = TH1D('b','b',nbin,low,high)
        qcd = TH1D('b','b',nbin,low,high)
        Variables = {}
        print "Variable: "+var
        for tre in List:
                histName = var+tre
#                print histName
                Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
		if var == "PhoJet_dPhi":
                        treepj[tre].Draw("TMath::Abs(PhoJet_dPhi)>>"+histName,weight)
                else:
                        treepj[tre].Draw(var+">>"+histName,weight)
#                Variables[tre].Sumw2()
                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
                if tre.startswith("photon_jet"):
                        backs.Add(Variables[tre])
                        phojet.Add(Variables[tre])
                if tre.startswith("qcd"):
                        backs.Add(Variables[tre])
                        qcd.Add(Variables[tre])
		if tre.startswith("dipho"):
                        backs.Add(Variables[tre])
#	backs.Draw("hist")
#        Variables['df1'].Draw("")
#        Variables['data'].Draw("")
#        Variables['df0p1'].Draw("e1same")
#        Variables['df0p001'].Draw("e1same")

        can = TCanvas("can3","can3", 1000,600)
        sm=0
        gr = {}
	updown = ["down"]
        for dec in updown:
                print dec
                sm+=1
                gr[dec] = TMultiGraph()
                gr[dec].SetMaximum(100)
                gr[dec].SetMinimum(0.001)
                gr[dec].SetTitle(";"+xtitle+";Signal/#sqrt{Signal+Backgraound}")
                grap = {}
                ledg = TLegend(0.7, 0.7, 0.9, 0.9)

                for type in darkFactors:
                        binvals = []
                        xaxis = []
                        pint = 0
                        maxval = 0
                        for bin in range(0,nbin):
                                mult = low+float(bin)*float(high-low)/float(nbin)
                                if dec == "down":
                                        if backs.Integral(0,bin) !=0:
                                                xaxis.append(mult)
                                                val = Variables[type].Integral(0,bin)/sqrt(Variables[type].Integral(0,bin)+backs.Integral(0,bin))
                                                print str(val)+"\t"+str(mult)
                                                binvals.append(val)
                                                if(val > pint):
                                                        pint = val
                                                        maxval = bin
                                else:
                                        if backs.Integral(bin,nbin+1) !=0:
                                                xaxis.append(mult)
                                                val = Variables[type].Integral(bin,nbin+1)/sqrt(Variables[type].Integral(bin,nbin+1)+backs.Integral(bin,nbin+1))
                                                print str(val)+"\t"+str(mult)
                                                binvals.append(val)
						if(val >pint):
							pint = val
							maxval = bin
                        print type+"\t"+dec+"\t"+str(len(xaxis))+"\t"+str(maxval)+"\t"+str(low+float(maxval)*float(high-low)/float(nbin))+"\t"+str(pint)
                        grap[type]=ROOT.TGraph(len(xaxis),scipy.array(xaxis),scipy.array(binvals))
                        grap[type].SetLineWidth(2)
                        grap[type].SetLineColor(colors[type])
                        gr[dec].Add(grap[type],"L")
                        ledg.AddEntry(grap[type],"#alpha_{Dark} = #alpha_{EM} x"+str(darks[type]),"l")
        gr[dec].Draw("A");
	gr[dec].GetXaxis().SetLimits(low,high)
        ledg.SetFillColor(0)
#        ledg.Draw()
        gPad.SetLogy(1)
        can.SaveAs("plots/optimized-"+region+"-"+var+".pdf")
