from LoadData import *
from plot_fun import *
from plot_variables import *

#List of variables to plot. Find names in plot_variables.py
vari = [Jet_pt]

#Dictionary of cuts and weights to apply. Each weight is set to 1.0 in Data.
#cuts = {'noCuts-NConst80_1':"((event_weight*kfactor*IDScaleFactor)*((isOverlap==0) && (Jet_NConst80==1)))"}
cuts = {'noCuts':"((event_weight*kfactor*IDScaleFactor)*(isOverlap==0))"}

#Main stack plotting function
def plot_var(var,title,xtitle,nbin,low,high,region):
	weight = cuts[region]
	can = TCanvas("can","can", 1000,900)
	can.SetBottomMargin(0.3)
	can.SetRightMargin(0.06)
	can.cd()
	gPad.SetLogy(1)
	gStyle.SetOptStat(0)
	
	#Collecting datasets by process
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
	for tre in List: #Loop through datasets.
		print tre
		histName = var+tre
		Variables[tre] = TH1F(histName,";"+xtitle+";"+"Events",nbin,low,high)
		treepj[tre].Draw(var+">>"+histName,weight)
		Variables[tre].Sumw2()
		scalable = lumi*crossx[tre]*kFact[tre]/Nevents[tre] #Weight each sample according to LoadData.py
		Variables[tre].Scale(scalable)
		if tre.startswith("df"): #Signal Samples
			Variables[tre].SetLineColor(colors[tre])
			Variables[tre].SetLineWidth(3)
		if tre.startswith("GJets"): #Photon + jet samples
			backs.Add(Variables[tre])
			phojet.Add(Variables[tre])
		if tre.startswith("QCD"): #QDC Multijet samples
			backs.Add(Variables[tre])
			qcd.Add(Variables[tre])
		if tre.startswith("DiPhoton"): #DiPhoton samples
			Variables[tre].SetFillColor(3)
			Variables[tre].SetLineColor(3)
			backs.Add(Variables[tre])
	Variables['Data'] = TH1F('Data',";"+xtitle+";"+"Events",nbin,low,high)
	treepj['Data'].Draw(var+">>Data",weight,"goff")
	Variables['Data'].SetLineColor(1)
	Variables['Data'].SetLineWidth(2)

	#Combine backgrounds into stack.
	stack.Add(Variables['DiPhoton'])
	stack.Add(qcd)
	stack.Add(phojet)
	stack.SetMaximum(1000000.0)
	stack.SetMinimum(0.001)
	stack.Draw("hist")
	stack.GetXaxis().SetLabelSize(0)
	backs.Draw("e2same")
	Variables['Data'].Draw("e1same")
	Variables['df1en0'].Draw("histsame")		
	Variables['df1en1'].Draw("histsame")		
	Variables['df1en2'].Draw("histsame")		
	Variables['df1en3'].Draw("histsame")
	
	#Set up and draw Pull plot.
	pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0)
	pad.SetTopMargin(0.7)
	pad.SetFillColor(0)
	pad.SetGridy(1)
	pad.SetFillStyle(0)
	pad.Draw()
	pad.cd(0)
	pad.SetRightMargin(0.06)
	Pull  = TH1D('Pull', 'Pull', nbin, low, high)
	Pull = Variables['Data'].Clone()
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
	for i in range(nbin):
		i += 1
		if Variables['Data'].GetBinContent(i) != 0:
			Pull.SetBinContent(i,Pull.GetBinContent(i))
		else: Pull.SetBinContent(i,0)
	Pull.Draw("e1")
	line.Draw('same')
	latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.035)
        latex2.SetTextAlign(31)
        latex2.DrawLatex(0.87, 0.95, "Work In Progress, "+str(lumi)+" fb^{-1} at #sqrt{s} = 13 TeV");
	print region
	print "********Yields***********"
        print "Photon + Jet Full         "+str(phojet.Integral(0,nbin+1))
        print "QCD                       "+str(qcd.Integral(0,nbin+1))
        print "DiPhoton         "+str(Variables['DiPhoton'].Integral(0,nbin+1))
	print "Total Background"+str(backs.Integral(0,nbin+1))
        print "Data            "+str(Variables['Data'].Integral(0,nbin+1))
        print "df = 10^0       "+str(Variables["df1en0"].Integral(0,nbin+1))
        print "df = 10^-1      "+str(Variables["df1en1"].Integral(0,nbin+1))
        print "df = 10^-2      "+str(Variables["df1en2"].Integral(0,nbin+1))
        print "df = 10^-3      "+str(Variables["df1en3"].Integral(0,nbin+1))

	led = TLegend(0.7, 0.7, 0.94, 0.9)
	led.AddEntry(phojet,"#gamma + Jet",'f')
	led.AddEntry(qcd,"QCD DiJet",'f')
	led.AddEntry(Variables['DiPhoton'],"#gamma + #gamma",'f')
	led.AddEntry(Variables['df1en0'],"#alpha_{Dark} = #alpha_{EM} x 10^{0}")
	led.AddEntry(Variables['df1en1'],"#alpha_{Dark} = #alpha_{EM} x 10^{-1}")
	led.AddEntry(Variables['df1en2'],"#alpha_{Dark} = #alpha_{EM} x 10^{-2}")
	led.AddEntry(Variables['df1en3'],"#alpha_{Dark} = #alpha_{EM} x 10^{-3}")
	led.AddEntry(Variables['Data'],"Data")
	led.SetFillColor(0)
	led.Draw("same")
	can.SaveAs("plots/"+region+"-"+var+".pdf")
	can.Close()
	phojet.Delete()
	qcd.Delete()


#Run plotting function.
for var in vari:
	for cal in cuts:
		plot_var(var[0],var[1],var[2],var[3],var[4],var[5],cal)
#		normalized(var[0],var[1],var[2],var[3],var[4],var[5],cal,cuts[cal])
#		optimizer(var[0],var[1],var[2],var[3],var[4],var[5],cal,cuts[cal])
