#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import math

# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))
import ROOT

directory="XSEC_Tutorial_Outputs_Data/"
if not os.path.exists(directory):
    os.makedirs(directory)
    
#Files: Unnormalized.root, dar_dump.txt
#Unnormalized.root contains the smear matrix and the reco matrix
#dar_dump is the marleydump for Decay at rest

#Also needed: RooUnfold


# In[2]:


mass = 0.000931494061 * 113428.9267 #TODO: Source of these numbers?
m_neg4 = mass ** -4

def PDF(energy):
    return 96 * energy*energy * m_neg4*(mass - (2*energy)) #from Marley Decay At Rest PDF


# In[3]:


'''Made w/ AnalysisScripts/ana_plothits/marEnergyDistrib2.cxx'''


infile = ROOT.TFile("unnormalized.root", "READ")

recoHist = infile.Get("recoMar")
smear = infile.Get("smear") #the 2d histo, but not normalized; In format (true, reco)


binSize = .5


# In[4]:



Response_Matrix = ROOT.RooUnfoldResponse (150, 0.0, 75)

#can assign weight of 1
#Put smear matrix in anticipated Response Matrix form
for x in range(smear.GetNbinsX()):
    for y in range(smear.GetNbinsY()):
        reco = y * binSize
        true = x * binSize
        for entry in range (int(smear.GetBinContent(x, y))):
            if smear.GetBinContent(x, y) != 0:
                Response_Matrix.Fill(reco + .001, true + .001, 1.0) #the .001 is just to get value off of bin edge


# In[5]:


iterations = 4
RooUnfoldBayes_Data    = ROOT.RooUnfoldBayes(Response_Matrix, recoHist, iterations)

error_method = 2

Unfolded_Data_Bayes_h  = RooUnfoldBayes_Data.Hunfold(error_method) 
Unfolded_Data_Bayes_h.SetTitle("Bayes Unfolded Energy")
Unfolded_Data_Bayes_h.SetXTitle("Energy (MeV)")
Unfolded_Data_Bayes_h.SetYTitle("Frequency")


# In[6]:



dumpxs=[]
dumpys=[]

flXscx = []
flXscy = []

dumpGraph = ROOT.TGraph()

'''Dar_dump is made via ./marleydumpxs dar_dump.txt [JSON CONFIG FILE (marley_dar)]'''

#dar dump is the true xsection
with open("dar_dump.txt", 'r') as info:
    for line in info.readlines():
        spaceIndex = line.find(' ')
        xpos = float(line[0:spaceIndex]) #the energy
        ypos = float(line[spaceIndex+1:])
        dumpxs.append(xpos)
        dumpys.append(ypos * (10.0**-4))
        if (xpos <= 52.82):
            flXscx.append(xpos)
            flXscy.append(ypos *(10.0**-42) * PDF(xpos)) #dar dump y is in 10 -42s. Result: trueXSec_i * flux_i
        

i=0
for entry in dumpxs:
    dumpGraph.SetPoint(i, dumpxs[i], dumpys[i])
    i = i+1

dumpGraph.SetLineColor(38)
dumpGraph.Draw("ALP")
dumpGraph.SetName("dumpGraph")
dumpGraph.SetName("MARLEY True XSec")


# In[7]:


xSec = ROOT.TH1D("xSec", "Cross Section", 150, 0, 75)
xGraph = ROOT.TGraphAsymmErrors()
xGraph.SetTitle("Calculated XSec")


numEvents = 10000
liveTime = numEvents/293.83
numArkg = 750
numMols = (numArkg * 1000)/39.948 #39.948 g/mol is the molar mass of Argon
avoConst =  6.02214076e23
numTargets = numMols * avoConst #to get number of nuclei


integral = 0

Rt = np.trapz(flXscy, x=flXscx) #Integrate the flux*Xsec

for binx in range (1, xSec.GetNbinsX() + 1):
    energy = binx * binSize
    if (PDF(energy)<0 and energy > 0):
        break
        
    pdfVal = []
    pdfEn = []
    
    for itr in range (0, 11):
        en = energy - .5 + (itr *.05)
        pdfVal.append(PDF(en))
        pdfEn.append(en)
    intFluxNumer = np.trapz(pdfVal, x = pdfEn) #Numerator: flux over this bin range.
    
    integratedFlux = intFluxNumer/Rt
    
    Efficiency_h   = 1 #TODO: set effeciency
    
    numObserved = Unfolded_Data_Bayes_h.GetBinContent(binx)
    errorHi= Unfolded_Data_Bayes_h.GetBinErrorUp(binx)/(Efficiency_h*integratedFlux*numEvents)
    errorLo=Unfolded_Data_Bayes_h.GetBinErrorLow(binx)/(Efficiency_h*integratedFlux*numEvents)

    if (Efficiency_h*integratedFlux != 0):
        sigma = numObserved/(Efficiency_h*integratedFlux*numEvents)
        xSec.SetBinContent(binx, sigma)
        xGraph.SetPoint(binx, energy, sigma* (10**38)) #Plot, scaling.
        xGraph.SetPointEYlow(binx, errorLo* (10**38))
        xGraph.SetPointEYhigh(binx, errorHi* (10**38))


# In[8]:


multiGraph = ROOT.TMultiGraph()


xGraph.SetLineColor(2)
xGraph.Draw("ALP")
xGraph.SetName("XSecGraph")

multiGraph.Add(xGraph)
multiGraph.Add(dumpGraph)
#You can plot xsection against anything you'd like, just keep scaling in mind

multiGraph.SetName("multigraph")
multiGraph.SetTitle("Cross Sections;MeV;Cross Section (10^-38 cm^2/nucleus)")
#multiGraph.BuildLegend()

outfile = ROOT.TFile("XSEC_Tutorial_Outputs_Data/XSection.root", "RECREATE")
outfile.cd()
Unfolded_Data_Bayes_h.Write()
xSec.Write()
xGraph.Write()
multiGraph.Write()
outfile.Write()


# In[9]:


outfile.Close()


# In[ ]:




