import plot_thesis

#For summing the event yields with applied systematic uncertainties
def systematic():
	for tre in List:
		Nominal = 0
		Jet_JEC_u = 0
		Jet_JEC_d = 0
		Pho_Scale_u = 0
		Pho_Scale_d = 0
		Pho_IDScale_u = 0
		Pho_IDScale_d = 0
		for event in treefin[tre]:
			Nominal += event.Nominal
			Jet_JEC_u += event.Jet_JEC_u
			Jet_JEC_d += event.Jet_JEC_d
			Pho_Scale_u += event.Pho_Scale_u
			Pho_Scale_d += event.Pho_Scale_d
			Pho_IDScale_u += event.Pho_IDScale_u
			Pho_IDScale_d += event.Pho_IDScale_d
		Nominal *=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Jet_JEC_u*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Jet_JEC_d*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Pho_Scale_u*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Pho_Scale_d*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Pho_IDScale_u*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
		Pho_IDScale_d*=lumi*crossx[tre]*kFact[tre]/Nevents[tre]
	print "Jet Energy Correction: "+str(Jet_JEC_u/Nominal)+" / "+str(Jet_JEC_d/Nominal)
	print "Photon Energy Correction: "+str(Pho_Scale_u/Nominal)+" / "+str(Pho_Scale_d/Nominal)
	print "Photon ID Scale: "+str(Pho_IDScale_u/Nominal)+" / "+str(Pho_IDScale_d/Nominal)
