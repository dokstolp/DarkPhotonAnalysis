from LoadData import *
from plot_fun import *
from plot_variables import *

#vari = [Jet_NchargedHad0,Jet_NConst,Jet_NchargedHad025,Jet_NchargedHad05,Jet_NchargedHad075,Jet_NchargedHad1,Jet_NchargedHad2,Jet_NchargedHad5,Jet_NchargedHad7,Jet_NchargedHad10,Jet_Nelectrons0,Jet_NConst,Jet_Nelectrons025,Jet_Nelectrons05,Jet_Nelectrons075,Jet_Nelectrons1,Jet_Nelectrons2,Jet_Nelectrons5,Jet_Nelectrons7,Jet_Nelectrons10,Jet_Nmuons0,Jet_NConst,Jet_Nmuons025,Jet_Nmuons05,Jet_Nmuons075,Jet_Nmuons1,Jet_Nmuons2,Jet_Nmuons5,Jet_Nmuons7,Jet_Nmuons10,Jet_Ngamma0,Jet_NConst,Jet_Ngamma025,Jet_Ngamma05,Jet_Ngamma075,Jet_Ngamma1,Jet_Ngamma2,Jet_Ngamma5,Jet_Ngamma7,Jet_Ngamma10,Jet_NneutralHad0,Jet_NConst,Jet_NneutralHad025,Jet_NneutralHad05,Jet_NneutralHad075,Jet_NneutralHad1,Jet_NneutralHad2,Jet_NneutralHad5,Jet_NneutralHad7,Jet_NneutralHad10]
#vari = [Jet_NchargedHad0,Jet_NConst,Jet_NchargedHad025,Jet_NchargedHad05,Jet_NchargedHad075,Jet_NchargedHad1,Jet_NchargedHad2,Jet_NchargedHad5,Jet_NchargedHad7,Jet_NchargedHad10,Jet_Nelectrons0,Jet_Ngamma0,Jet_NConst,Jet_Ngamma025,Jet_Ngamma05,Jet_Ngamma075,Jet_Ngamma1,Jet_Ngamma2,Jet_Ngamma5,Jet_Ngamma7,Jet_Ngamma10]
#vari = [Jet_PhoConvLxy]
#vari = [Jet_NConstituents0,Jet_NConst,Jet_NConstituents025,Jet_NConstituents05,Jet_NConstituents075,Jet_NConstituents1,Jet_NConstituents2,Jet_NConstituents5,Jet_NConstituents7,Jet_NConstituents10]
vari = [Jet_NConstituents025]
#vari = [Jet_MinPt,Jet_MinEt]
vari2 = [Jet_NConst,nJets]
cuts = {'noCuts-DiScale':"((neg_weight*event_weight*(isOverlap==0)) && Jet_isPixelPho == 1 && Jet_PhoConvLxy<140)"}
#cuts = {'noCuts-noWeight':"(event_weight)"}


def plot_var(var,title,xtitle,nbin,low,high,region):
	cutter = cuts[region]
	weight = cutter
#	if cutter == "":
#		weight = "event_weight"
#	else:
#		weight = "("+cutter+")*event_weight"
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
#		print tre
		histName = var+tre
		Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
		if tre == "dipho":
			weight = cutter+"&& (JetPho_dR<0.5)"
#			treepj[tre].Draw("TMath::Abs(PhoJet_dPhi)>>"+histName,weight)
#			treepj[tre].Draw(var+">>"+histName,diphoCut)
#		else: 
#		if var == "PhoJet_dPhi":
#			treepj[tre].Draw("TMath::Abs(PhoJet_dPhi)>>"+histName,weight)
#		else:
#			treepj[tre].Draw(var+">>"+histName,diphoCut)

		treepj[tre].Draw(var+">>"+histName,weight)
		Variables[tre].Sumw2()
		scalable = lumi*crossx[tre]*kFact[tre]/Nevents[tre]
#		scalable = 1/Nevents[tre]
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
	if var == "PhoJet_dPhi":
		treepj['data'].Draw("TMath::Abs(PhoJet_dPhi)>>data",cutter,"goff")
	else:
		treepj['data'].Draw(var+">>data",cutter,"goff")

	Variables['data'].SetLineColor(1)
	Variables['data'].SetLineWidth(2)
	stack.Add(Variables['dipho'])
	stack.Add(qcd)
	stack.Add(phojet)
		
	stack.SetMaximum(1000000.0)
#	stack.SetMaximum(3000.0)
	stack.SetMinimum(0.001)
	stack.Draw("hist")
	stack.GetXaxis().SetLabelSize(0)
	backs.Draw("e2same")
	Variables['data'].Draw("e1same")
#	Variables['df1en0'].SetLineColor(4)		
#	Variables['df1en1'].SetLineColor(6)		
#	Variables['df1en2'].SetLineColor(16)		
#	Variables['df1en3'].SetLineColor(2)		
	Variables['df1en0'].Draw("histsame")		
	Variables['df1en1'].Draw("histsame")		
	Variables['df1en2'].Draw("histsame")		
	Variables['df1en3'].Draw("histsame")		
#	print Variables['data'].GetMaximum()
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
#	Pull.Add(backs,-1)
	Pull.Divide(backs)
	Pull.SetMarkerStyle(20)
#	Pull.SetMaximum(10.0)
#	Pull.SetMaximum(20.0)
	Pull.SetMaximum(2.0)
#	Pull.SetMinimum(-10.0)
	Pull.SetMinimum(0.0)
	Pull.SetFillColor(2)
	Pull.GetXaxis().SetTitle(xtitle)
	Pull.GetYaxis().SetTitleSize(0.04)
#	Pull.GetYaxis().SetTitle('#sigma(Data-MC)')
	Pull.GetYaxis().SetTitle('Data/MC')
	Pull.SetMarkerSize(0.7)
	Pull.GetYaxis().SetNdivisions(5)
	line = TF1("line","1",low,high)
	line.SetLineColor(2)
	line.SetLineWidth(2)
	for i in range(nbin):
		i += 1
		if Variables['data'].GetBinContent(i) != 0:
#			Pull.SetBinContent(i,Pull.GetBinContent(i)/Pull.GetBinError(i))
			Pull.SetBinContent(i,Pull.GetBinContent(i))
		else: Pull.SetBinContent(i,0)
#	Pull.Draw("HIST")
	Pull.Draw("e1")
	line.Draw('same')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31) # align right
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	print region
	print "********Yields***********"
        print "Photon + Jet Full         "+str(phojet.Integral(0,nbin+1))
        print "QCD                       "+str(qcd.Integral(0,nbin+1))
        print "DiPhoton         "+str(Variables['dipho'].Integral(0,nbin+1))
	print "Total Background"+str(backs.Integral(0,nbin+1))
        print "Data            "+str(Variables['data'].Integral(0,nbin+1))
        print "df = 10^0       "+str(Variables["df1en0"].Integral(0,nbin+1))
        print "df = 10^-1      "+str(Variables["df1en1"].Integral(0,nbin+1))
        print "df = 10^-2      "+str(Variables["df1en2"].Integral(0,nbin+1))
        print "df = 10^-3      "+str(Variables["df1en3"].Integral(0,nbin+1))
#        print "df = 10^0       "+str(round(Variables["df1en0"].Integral(0,nbin+1)*100,2))
#        print "df = 10^-1      "+str(round(Variables["df1en1"].Integral(0,nbin+1)*100,2))
#        print "df = 10^-2      "+str(round(Variables["df1en2"].Integral(0,nbin+1)*100,2))
#        print "df = 10^-3      "+str(round(Variables["df1en3"].Integral(0,nbin+1)*100,2))

	led = TLegend(0.7, 0.7, 0.94, 0.9)
#	led = TLegend(0.15, 0.7, 0.39, 0.9)
#	led = TLegend(0.6, 0.6, 0.9, 0.9)
	led.AddEntry(phojet,"#gamma + Jet",'f')
	led.AddEntry(qcd,"QCD DiJet",'f')
	led.AddEntry(Variables['dipho'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"#alpha_{Dark} = #alpha_{EM} x 10^{0}")
	led.AddEntry(Variables['df1en1'],"#alpha_{Dark} = #alpha_{EM} x 10^{-1}")
	led.AddEntry(Variables['df1en2'],"#alpha_{Dark} = #alpha_{EM} x 10^{-2}")
	led.AddEntry(Variables['df1en3'],"#alpha_{Dark} = #alpha_{EM} x 10^{-3}")
	led.AddEntry(Variables['data'],"Data")
	led.SetFillColor(0)
	led.Draw("same")
#	can.SaveAs(region+"/"+region+"-"+var+".pdf")
	can.SaveAs("plots/"+region+"-"+var+".pdf")
	can.Close()
	phojet.Delete()
	qcd.Delete()



for var in vari:
	for cal in cuts:
		plot_var(var[0],var[1],var[2],var[3],var[4],var[5],cal)
#		normalized(var[0],var[1],var[2],var[3],var[4],var[5],cal,cuts[cal])
#		optimizer(var[0],var[1],var[2],var[3],var[4],var[5],cal,cuts[cal])
	

#for cal in cuts:
#	aXe(cal,cuts[cal])

#x = vari2[0]
#y = vari2[1]
#for cal in cuts:
#	plot_2D(x[0],y[0],x[1]+y[1],x[2],y[2],x[3],x[4],x[5],y[3],y[4],y[5],cal,cuts[cal])
