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
                Variables[tre] = TH2F(histName,";"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
                treepj[tre].Draw(yvar+":"+xvar+">>"+histName,weight)
                Variables[tre].Sumw2()
                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
                if tre.startswith("GJets"):
                        backs.Add(Variables[tre])
                if tre.startswith("QCD"):
                        backs.Add(Variables[tre])
                if tre.startswith("DiPhoton"):
                        backs.Add(Variables[tre])
	backs.Draw("COLZ")
	can.SaveAs("plots/Background"+region+"-"+xvar+"-"+yvar+".pdf")

        Variables['Data'] = TH2D('Data',"Data/Background;"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
        treepj['Data'].Draw(yvar+":"+xvar+">>Data",cutter)
        Pull  = TH2D('Pull',title+";"+xtitle+";"+ytitle,xnbin,xlow,xhigh,ynbin,ylow,yhigh)
        Pull = Variables['data'].Clone()
        Pull.Divide(backs)
        Pull.SetMaximum(2.0)
        Pull.SetMinimum(0.0)
        for i in range(xnbin):
                i += 1
		for j in range(ynbin):
			j += 1
			if Variables['data'].GetBinContent(i,j) != 0:
                        	Pull.SetBinContent(i,j,Pull.GetBinContent(i,j))
#				print "x: "+str(i)+" y: "+str(j)+" Data: "+str(Variables['data'].GetBinContent(i,j))+" Backgrounds: "+str(backs.GetBinContent(i,j))+" Data/Backgrounds: "+str(Pull.GetBinContent(i,j))
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
                VariablesNum[tre].Divide(VariablesDen[tre])
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


def optimizer(var,title,xtitle,nbin,low,high,region,cutter):
        darknames = {'df1en0':'Dark Factor = 1E0','df1en1':'Dark Factor = 1E-1','df1en2':'Dark Factor = 1E-2','df1en3':'Dark Factor = 1E-3'}
	print "*************"+region+"*****************"
	weight = cutter
        can = TCanvas("can","can", 1000,700)
        backs = TH1D('a','a',nbin,low,high)
        phojet = TH1D('b','b',nbin,low,high)
        qcd = TH1D('b','b',nbin,low,high)
        Variables = {}
        print "Variable: "+var
        for tre in List:
                histName = var+tre
                Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
                treepj[tre].Draw(var+">>"+histName,weight)
                Variables[tre].Sumw2()
                Variables[tre].Scale(lumi*crossx[tre]*kFact[tre]/Nevents[tre])
                if tre.startswith("GJets"):
                        backs.Add(Variables[tre])
                        phojet.Add(Variables[tre])
                if tre.startswith("QCD"):
                        backs.Add(Variables[tre])
                        qcd.Add(Variables[tre])
		if tre.startswith("DiPhoton"):
                        backs.Add(Variables[tre])
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
                        ledg.AddEntry(grap[type],"f_{D} = "+darknames[type],"l")
        gr[dec].Draw("A");
	gr[dec].GetXaxis().SetLimits(low,high)
        ledg.SetFillColor(0)
        gPad.SetLogy(1)
        can.SaveAs("plots/optimized-"+region+"-"+var+".pdf")
