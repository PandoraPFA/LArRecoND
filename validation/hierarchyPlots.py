# Create hierarchy validation plots (eff, purity, completeness, vtx etc.)

import argparse
import array
import math
import os
import ROOT
import sys

class parameters(object):

    def __init__(self, mcFileName, mcTreeName, evtFileName, evtTreeName):
        # MC hierarchy file and tree names
        self.mcFileName = mcFileName
        self.mcTreeName = mcTreeName

        # Event hierarchy file and tree names
        self.evtFileName = evtFileName
        self.evtTreeName = evtTreeName

        # Output histogram filenames
        self.histMCFileName = mcFileName.replace('.root', '_Histos.root')
        self.histEvtFileName = evtFileName.replace('.root', '_Histos.root')
        
        # Quality cuts
        self.minCompleteness = 0.1
        self.minPurity = 0.5
        self.minNSharedHits = 5

        # Max number of hits & momentum (GeV/c)
        self.maxNHits = 1.0e4
        self.maxMtm = 5.0

        # List of particles
        self.pdgList = [(11, 'electron'), (13, 'muon'), (22, 'photon'),
                        (211, 'piplus'), (-211, 'piminus'), (2212, 'proton')]


class histoMCList(object):

    # Object storing the histograms for a given particle type
    def __init__(self, hitsAll, hitsEff, mtmAll, mtmEff, completeness, purity):
        self.hitsAll = hitsAll
        self.hitsEff = hitsEff
        self.mtmAll = mtmAll
        self.mtmEff = mtmEff
        self.completeness = completeness
        self.purity = purity
        

class histoVtxList(object):

    # Object storing the primary vertex dX, dY, dZ and dR histograms
    def __init__(self, hVtxDX, hVtxDY, hVtxDZ, hVtxDR):
        self.hVtxDX = hVtxDX
        self.hVtxDY = hVtxDY
        self.hVtxDZ = hVtxDZ
        self.hVtxDR = hVtxDR
        

class crystalBallFun(object):

    def __call__(self, xArr, pars):

        # Crystal Ball function for vertex residual fits
        x = xArr[0]
        norm = pars[0]
        x0 = pars[1]
        sigmaL = abs(pars[2])
        sigmaR = abs(pars[3])
        alphaL = abs(pars[4])
        nL = abs(pars[5])
        alphaR = abs(pars[6])
        nR = abs(pars[7])

        t = 0.0
        if x < x0:
            t = (x - x0)/sigmaL
        else:
            t = (x - x0)/sigmaR
        
        value = 0.0
        if (t < alphaL):
            value = self.getTail(t, alphaL, nL)
        elif (t <= alphaR):
            value = math.exp(-0.5*t*t)
        else:
            value = self.getTail(-t, alphaR, nR)

        return value*norm

    def getTail(self, t, alpha, n):

        result = 0.0
        if abs(alpha) > 0.0:
            a = math.pow(abs(n/alpha), n) * math.exp(-0.5*alpha*alpha)
            b = (n/alpha) - alpha
            result = a/math.pow(abs(b-t), n)
        return result

 
def getParticleType(mcPDG):

    name = 'Unknown'
    absPDG = abs(mcPDG)
    if absPDG == 13:
        name = 'muon'
    elif absPDG == 11:
        name = 'electron'
    elif absPDG == 2212:
        name = 'proton'
    elif absPDG == 22:
        name = 'photon'
    elif mcPDG == 211:
        name = 'piplus'
    elif mcPDG == -211:
        name = 'piminus'

    return name


def defineMCHistos(pars):

    # Create empty MC histograms for each particle type.
    # Store them in the histogram map, which is returned.
    # Map key = particle type, value = histogram list object
    histMCMap = {}

    # hits binning
    nHitBins = 35
    nHitBinEdges = nHitBins + 1
    hitsBinning = array.array('d', [0.0]*nHitBinEdges)
    for iB in range(nHitBinEdges):
        edge = math.pow(10.0, 1.0 + (iB*1.0 + 2.0)*0.1)
        hitsBinning[iB] = edge
        
    # momentum binning
    nMtmBins = 26
    mtmBinning = array.array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4,
                                   1.6, 2.0, 2.4, 2.8, 3.4, 4.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0])

    # Loop over the particle types
    for iPDG,pLabel in pars.pdgList:
        print('Defining histos for PDG = {0}'.format(pLabel))

        # all hits, hits efficiency
        hHitsAll = ROOT.TH1F('{0}_HitsAll'.format(pLabel), '', nHitBins, hitsBinning)
        hHitsAll.SetDirectory(0)
        hHitsAll.GetXaxis().SetTitle('Number of Hits')
        hHitsAll.GetYaxis().SetTitle('Number of Candidates')
        
        hHitsEff = ROOT.TH1F('{0}_HitsEff'.format(pLabel), '', nHitBins, hitsBinning)
        hHitsEff.SetDirectory(0)
        hHitsEff.GetXaxis().SetTitle('Number of Hits')
        hHitsEff.GetYaxis().SetTitle('Reconstruction Efficiency')

        # all momentum, momentum efficiency
        hMtmAll = ROOT.TH1F('{0}_MtmAll'.format(pLabel), '', nMtmBins, mtmBinning)
        hMtmAll.SetDirectory(0)
        hMtmAll.GetXaxis().SetTitle('True Momentum [GeV]')
        hMtmAll.GetYaxis().SetTitle('Number of Candidates')

        hMtmEff = ROOT.TH1F('{0}_MtmEff'.format(pLabel), '', nMtmBins, mtmBinning)
        hMtmEff.SetDirectory(0)
        hMtmEff.GetXaxis().SetTitle('True Momentum [GeV]')
        hMtmEff.GetYaxis().SetTitle('Reconstruction Efficiency')

        # Completeness
        hCompleteness = ROOT.TH1F('{0}_Completeness'.format(pLabel), '', 51, -0.01, 1.01)
        hCompleteness.SetDirectory(0)
        hCompleteness.GetXaxis().SetTitle('Completeness')
        hCompleteness.GetYaxis().SetTitle('Fraction of Events')

        # Purity
        hPurity = ROOT.TH1F('{0}_Purity'.format(pLabel), '', 51, -0.01, 1.01)
        hPurity.SetDirectory(0)
        hPurity.GetXaxis().SetTitle('Purity')
        hPurity.GetYaxis().SetTitle('Fraction of Events')

        histMCMap[pLabel] = histoMCList(hHitsAll, hHitsEff, hMtmAll, hMtmEff,
                                        hCompleteness, hPurity)

    # Return the histogram map
    return histMCMap


def defineVtxHistos(pars, label = 'all'):

    # Define primary vertex histograms for all interactions
    hVtxDX = ROOT.TH1F('{0}VtxDX'.format(label), '', 100, -5.0, 5.0)
    hVtxDX.SetDirectory(0)
    hVtxDX.GetXaxis().SetTitle('Vertex #DeltaX [cm]')
    hVtxDX.GetYaxis().SetTitle('Number of Candidates')

    hVtxDY = ROOT.TH1F('{0}VtxDY'.format(label), '', 100, -5.0, 5.0)
    hVtxDY.SetDirectory(0)
    hVtxDY.GetXaxis().SetTitle('Vertex #DeltaY [cm]')
    hVtxDY.GetYaxis().SetTitle('Number of Candidates')

    hVtxDZ = ROOT.TH1F('{0}VtxDZ'.format(label), '', 100, -5.0, 5.0)
    hVtxDZ.SetDirectory(0)
    hVtxDZ.GetXaxis().SetTitle('Vertex #DeltaZ [cm]')
    hVtxDZ.GetYaxis().SetTitle('Number of Candidates')

    hVtxDR = ROOT.TH1F('{0}VtxDR'.format(label), '', 100, 0.0, 5.0)
    hVtxDR.SetDirectory(0)
    hVtxDR.GetXaxis().SetTitle('Vertex #DeltaR [cm]')
    hVtxDR.GetYaxis().SetTitle('Number of Candidates')

    # Object storing the list of histograms
    hVtxList = histoVtxList(hVtxDX, hVtxDY, hVtxDZ, hVtxDR)
    return hVtxList
    

def createMCHistos(pars):

    print('mcFileName = {0}, mcTreeName = {1}'.format(pars.mcFileName, pars.mcTreeName))

    # Define performance histograms for each particle type
    histMCMap = defineMCHistos(pars)

    # Open MC hierachy file and its tree
    mcFile = ROOT.TFile.Open(pars.mcFileName, 'read')
    mcTree = mcFile.Get(pars.mcTreeName)

    # Loop over MC hierarchy tree entries
    nMC = mcTree.GetEntries()
    nTotPass = 0
    nTotFail = 0
    nTotal = 0

    for i in range(nMC):

        # Get tree entry
        mcTree.GetEntry(i)

        event = getattr(mcTree, 'event')

        if i%10000 == 0:
            print('MC entry {0}'.format(nMC-i))
            
        # Get MC particle PDG Id
        mcPDG = getattr(mcTree, 'mcPDG')
        # Get particle name type
        mcType = getParticleType(mcPDG)

        # Require known particle type
        if mcType == 'Unknown':
            continue

        # Histogram list for given particle type
        hList = histMCMap[mcType]

        # Number of hits and true momentum
        hHitsAll = hList.hitsAll
        mcNHits = getattr(mcTree, 'mcNHits')
        hHitsAll.Fill(mcNHits)

        mcMtm = getattr(mcTree, 'mcMtm')
        nMatches = getattr(mcTree, 'nMatches')

        hMtmAll = hList.mtmAll
        hMtmAll.Fill(mcMtm)

        # Find best completeness MC-reco match and its corresponding purity
        completeVect = getattr(mcTree, 'completenessVector')
        purityVect = getattr(mcTree, 'purityVector')
        nSharedHitsVect = getattr(mcTree, 'nSharedHitsVector')
            
        complete = 0.0
        purity = 0.0
        nSharedHits = 0
        bestNShared = 0
        for iC,nShared in enumerate(nSharedHitsVect):
            if nShared > bestNShared:
                bestNShared = nShared
                complete = completeVect[iC]
                purity = purityVect[iC]

        # Check for minimum completeness, purity and number of shared hits
        if nMatches > 0 and complete >= pars.minCompleteness and purity >= pars.minPurity \
           and bestNShared >= pars.minNSharedHits:

            # Efficiency numerator (hits & true momentum)
            hHitsEff = hList.hitsEff
            hHitsEff.Fill(mcNHits)

            hMtmEff = hList.mtmEff
            hMtmEff.Fill(mcMtm)

            # Fill completeness and purity histos
            hCompleteness = hList.completeness
            hCompleteness.Fill(complete)

            hPurity = hList.purity
            hPurity.Fill(purity)


    # Process hits and momentum histograms to get their efficiencies,
    # and normalise the completeness and purity histograms.
    # Then, write these all to the histogram output file
    print('Creating {0}'.format(pars.histMCFileName))
    hMCOutFile = ROOT.TFile.Open(pars.histMCFileName, 'recreate')

    # Loop over particle types
    for pdgId, pLabel in pars.pdgList:

        # Histogram list for given particle type
        hList = histMCMap[pLabel]

        # Find efficiencies
        setEffHist(hList.hitsEff, hList.hitsAll)
        setEffHist(hList.mtmEff, hList.mtmAll)

        # Normalise completeness and purity
        nCompleteness = hList.completeness.GetEntries()*1.0
        if nCompleteness > 0.0:
            hList.completeness.Scale(1.0/nCompleteness)

        nPurity = hList.purity.GetEntries()*1.0
        if nPurity > 0.0:
            hList.purity.Scale(1.0/nPurity)

        # Write out the histograms
        hMCOutFile.cd()
        hList.hitsAll.Write()
        hList.hitsEff.Write()
        hList.mtmAll.Write()
        hList.mtmEff.Write()
        hList.completeness.Write()
        hList.purity.Write()
        
    # Close the files
    hMCOutFile.Close()
    mcFile.Close()


def createVtxHistos(pars):
    
    # Event vertex histograms
    hAllVtxList = defineVtxHistos(pars, 'all')
    # Positive and negative drift directions
    hPosVtxList = defineVtxHistos(pars, 'pos')
    hNegVtxList = defineVtxHistos(pars, 'neg')

    # Open event hierarchy file
    evtFile = ROOT.TFile.Open(pars.evtFileName, 'read')
    evtTree = evtFile.Get(pars.evtTreeName)

    # Loop over event entries
    nEvt = evtTree.GetEntries()
    for i in range(nEvt):

        # Get tree entry
        evtTree.GetEntry(i)

        if i%10000 == 0:
            print('Evt entry {0}'.format(nEvt-i))

        # Get vertex residuals
        vtxDx = getattr(evtTree, 'vtxDx')
        vtxDy = getattr(evtTree, 'vtxDy')
        vtxDz = getattr(evtTree, 'vtxDz')
        vtxDr = getattr(evtTree, 'vtxDr')

        # Plot vertex residuals for both drift directions. Assume drift (cathode to anode) module boundaries
        # anode - cathode - cathode - anode etc. are roughly between -350 cm to 350 cm in x every 50cm:
        # e.g neg drift dir -350 to -300, -250 to -200, -150 to -100, -50 to 0, 50 to 100, 150 to 200, 250 to 300
        # positive drift direction will be the other ranges -300 to -250, -200 to -150, etc
        # Neighbouring cathodes are placed next to each other
        trueVtxX = getattr(evtTree, 'trueVtxX')

        # Group drift direction within anode/cathode boundaries
        if (abs(trueVtxX) <= 350) :
            if (trueVtxX <= -300
                or (trueVtxX >= -250 and trueVtxX <= -200)
                or (trueVtxX >= -150 and trueVtxX <= -100)
                or (trueVtxX >= -50 and trueVtxX <= 0)
                or (trueVtxX >= 50 and trueVtxX <= 100)
                or (trueVtxX >= 150 and trueVtxX <= 200)
                or (trueVtxX >= 250 and trueVtxX < 300)):
                # Negative drift
                hNegVtxList.hVtxDX.Fill(vtxDx)
                hNegVtxList.hVtxDY.Fill(vtxDy)
                hNegVtxList.hVtxDZ.Fill(vtxDz)
                hNegVtxList.hVtxDR.Fill(vtxDr)
            else:
                # Positive drift
                hPosVtxList.hVtxDX.Fill(vtxDx)
                hPosVtxList.hVtxDY.Fill(vtxDy)
                hPosVtxList.hVtxDZ.Fill(vtxDz)
                hPosVtxList.hVtxDR.Fill(vtxDr)

        hAllVtxList.hVtxDX.Fill(vtxDx)
        hAllVtxList.hVtxDY.Fill(vtxDy)
        hAllVtxList.hVtxDZ.Fill(vtxDz)
        hAllVtxList.hVtxDR.Fill(vtxDr)

    # Write out histograms
    print('Creating {0}'.format(pars.histEvtFileName))
    hEvtOutFile = ROOT.TFile.Open(pars.histEvtFileName, 'recreate')
    hEvtOutFile.cd()

    hNegVtxList.hVtxDX.Write()
    hNegVtxList.hVtxDY.Write()
    hNegVtxList.hVtxDZ.Write()
    hNegVtxList.hVtxDR.Write()

    hPosVtxList.hVtxDX.Write()
    hPosVtxList.hVtxDY.Write()
    hPosVtxList.hVtxDZ.Write()
    hPosVtxList.hVtxDR.Write()

    hAllVtxList.hVtxDX.Write()
    hAllVtxList.hVtxDY.Write()
    hAllVtxList.hVtxDZ.Write()
    hAllVtxList.hVtxDR.Write()

    # Close files
    hEvtOutFile.Close()
    evtFile.Close()
        
            
def setEffHist(hEff, hAll):

    # Modify the numerator hEff histogram to be the efficiency
    # by dividing by the denominator hAll histogram bin contents.
    # Assumes both have the same binning, which they should do

    # Loop over bins, including under and overflow
    nBins = hEff.GetXaxis().GetNbins()
    for i in range(-1, nBins+1):
        i1 = i + 1
        num = hEff.GetBinContent(i1)
        denom = hAll.GetBinContent(i1)
        # Efficiency and its binomial error
        eff = (num/denom) if (denom > 0.0) else 0.0
        err = math.sqrt(eff*(1.0 - eff)/denom) if (denom > 0.0) else 0.0
        # Update the efficiency bin content
        hEff.SetBinContent(i1, eff)
        hEff.SetBinError(i1, err)


def plotMCHistos(pars):

    # Open MC histogram ROOT file
    hMCFile = ROOT.TFile.Open(pars.histMCFileName, 'read')

    # Define plotting canvas
    theCanvas = ROOT.TCanvas('theCanvas', '', 1350, 700)
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetOptStat(0)
    theCanvas.UseCurrentStyle()

    # For text labelling
    text = ROOT.TLatex()
    text.SetTextSize(0.075)
    text.SetNDC(True)

    # All interactions: muons, protons, piplus and piminus
    theCanvas.Divide(3,2)

    # Hits efficiency
    theCanvas.cd(1)
    muHitsEff = hMCFile.Get('muon_HitsEff')
    muHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    muHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    muHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    pHitsEff = hMCFile.Get('proton_HitsEff')
    pHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    pHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    pHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    eHitsEff = hMCFile.Get('electron_HitsEff')
    eHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    eHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    eHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, 'e')

    theCanvas.cd(4)
    pipHitsEff = hMCFile.Get('piplus_HitsEff')
    pipHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    pipHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    pipHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(5)
    pimHitsEff = hMCFile.Get('piminus_HitsEff')
    pimHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    pimHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    pimHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.cd(6)
    gHitsEff = hMCFile.Get('photon_HitsEff')
    gHitsEff.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    gHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx()
    gHitsEff.Draw()
    text.DrawLatex(0.775, 0.25, '#gamma')

    theCanvas.Print('hierarchy_allHitsEff.png')    

    # Hit populations (efficiency denominators)
    theCanvas.cd(1)
    muHitsAll = hMCFile.Get('muon_HitsAll')
    muHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    muHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, '#mu')

    theCanvas.cd(2)
    pHitsAll = hMCFile.Get('proton_HitsAll')
    pHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    pHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, 'p')

    theCanvas.cd(3)
    eHitsAll = hMCFile.Get('electron_HitsAll')
    eHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    eHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, 'e')

    theCanvas.cd(4)
    pipHitsAll = hMCFile.Get('piplus_HitsAll')
    pipHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    pipHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, '#pi^{+}')

    theCanvas.cd(5)
    pimHitsAll = hMCFile.Get('piminus_HitsAll')
    pimHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    pimHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, '#pi^{-}')

    theCanvas.cd(6)
    gHitsAll = hMCFile.Get('photon_HitsAll')
    gHitsAll.GetXaxis().SetRangeUser(0.0, pars.maxNHits)
    ROOT.gPad.SetLogx()
    gHitsAll.Draw()
    text.DrawLatex(0.5, 0.3, '#gamma')

    theCanvas.Print('hierarchy_allHits.png')

    # Momentum efficiency
    theCanvas.cd(1)
    muMtmEff = hMCFile.Get('muon_MtmEff')
    muMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    muMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    muMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    pMtmEff = hMCFile.Get('proton_MtmEff')
    pMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    pMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    pMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    eMtmEff = hMCFile.Get('electron_MtmEff')
    eMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    eMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    eMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, 'e')

    theCanvas.cd(4)
    pipMtmEff = hMCFile.Get('piplus_MtmEff')
    pipMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    pipMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    pipMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(5)
    pimMtmEff = hMCFile.Get('piminus_MtmEff')
    pimMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    pimMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    pimMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.cd(6)
    gMtmEff = hMCFile.Get('photon_MtmEff')
    gMtmEff.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    gMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    ROOT.gPad.SetLogx(0)
    gMtmEff.Draw()
    text.DrawLatex(0.775, 0.25, '#gamma')

    theCanvas.Print('hierarchy_allMtmEff.png')    

    # Momentum populations (efficiency denominator)
    theCanvas.cd(1)
    muMtmAll = hMCFile.Get('muon_MtmAll')
    muMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    muMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    pMtmAll = hMCFile.Get('proton_MtmAll')
    pMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    pMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    eMtmAll = hMCFile.Get('electron_MtmAll')
    eMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    eMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, 'e')

    theCanvas.cd(4)
    pipMtmAll = hMCFile.Get('piplus_MtmAll')
    pipMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    pipMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(5)
    pimMtmAll = hMCFile.Get('piminus_MtmAll')
    pimMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    pimMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.cd(6)
    gMtmAll = hMCFile.Get('photon_MtmAll')
    gMtmAll.GetXaxis().SetRangeUser(0.0, pars.maxMtm)
    ROOT.gPad.SetLogx(0)
    gMtmAll.Draw()
    text.DrawLatex(0.775, 0.25, '#gamma')

    theCanvas.Print('hierarchy_allMtm.png')

    # Close the histogram MC file
    hMCFile.Close()
    

def plotVtxHistos(pars):

    # Open event histogram ROOT file
    hEvtFile = ROOT.TFile.Open(pars.histEvtFileName, 'read')

    # Define plotting canvas
    theCanvas = ROOT.TCanvas('theCanvas', '', 900, 700)
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(111)
    theCanvas.UseCurrentStyle()

    # All interactions: primary vertex residuals
    theCanvas.Divide(2,2)

    # Define the Crystal Ball fit function. We need to keep this in scope
    # and not be deleted, so we pass it as an argument in setCBFun()
    cbFun = crystalBallFun()

    theCanvas.cd(1)
    vtxDX = hEvtFile.Get('allVtxDX')
    #vtxDXFun = setCBFun(vtxDX, cbFun, 0.0, 0.3)
    #vtxDX.Fit(vtxDXFun)
    vtxDX.Draw()
    ROOT.gPad.Update()
    
    theCanvas.cd(2)
    vtxDY = hEvtFile.Get('allVtxDY')
    #vtxDYFun = setCBFun(vtxDY, cbFun, 0.0, 0.3)
    #vtxDY.Fit(vtxDYFun)
    vtxDY.Draw()
    ROOT.gPad.Update()

    theCanvas.cd(3)
    vtxDZ = hEvtFile.Get('allVtxDZ')
    #vtxDZFun = setCBFun(vtxDZ, cbFun, 0.0, 0.3)
    #vtxDZ.Fit(vtxDZFun)
    vtxDZ.Draw()
    ROOT.gPad.Update()

    theCanvas.cd(4)
    vtxDR = hEvtFile.Get('allVtxDR')
    #vtxDRFun = setCBFun(vtxDR, cbFun, 0.3, 0.3)
    #vtxDR.Fit(vtxDRFun)
    vtxDR.Draw()
    ROOT.gPad.Update()

    theCanvas.Print('hierarchy_allVtx.png')

    # Negative drift direction vertex residuals
    theCanvas.cd(1)
    vtxDX = hEvtFile.Get('negVtxDX')
    vtxDX.Draw()

    theCanvas.cd(2)
    vtxDY = hEvtFile.Get('negVtxDY')
    vtxDY.Draw()

    theCanvas.cd(3)
    vtxDZ = hEvtFile.Get('negVtxDZ')
    vtxDZ.Draw()

    theCanvas.cd(4)
    vtxDR = hEvtFile.Get('negVtxDR')
    vtxDR.Draw()

    theCanvas.Print('hierarchy_negVtx.png')

    # Positive drift direction vertex residuals
    theCanvas.cd(1)
    vtxDX = hEvtFile.Get('posVtxDX')
    vtxDX.Draw()

    theCanvas.cd(2)
    vtxDY = hEvtFile.Get('posVtxDY')
    vtxDY.Draw()

    theCanvas.cd(3)
    vtxDZ = hEvtFile.Get('posVtxDZ')
    vtxDZ.Draw()

    theCanvas.cd(4)
    vtxDR = hEvtFile.Get('posVtxDR')
    vtxDR.Draw()

    theCanvas.Print('hierarchy_posVtx.png')

    # Clost the histogram event file
    hEvtFile.Close()


def setCBFun(hist, cbFun, theMean = -999.0, theSigma = -999.0):

    # Set the Crystal Ball fit function for the given histogram.
    # Note that we pass cbFun = crystalBallFun() as a callable argument
    # so that it remains in scope when the function is returned
    histName = hist.GetName()
    funName = '{0}Fun'.format(histName)
    xAxis = hist.GetXaxis()
    xMin = xAxis.GetXmin()
    xMax = xAxis.GetXmax()
    vtxDXFun = ROOT.TF1(funName, cbFun, xMin, xMax, 8)

    mean = theMean if (theMean > -999.0) else hist.GetMean()
    sigma = theSigma if (theSigma > -999.0) else hist.GetStdDev()
    norm = hist.GetMaximum()

    fun = ROOT.TF1(funName, cbFun, xMin, xMax, 8)
    fun.SetParameters(norm, mean, sigma, sigma, 1.2, 1.0, 1.2, 1.0)
    fun.SetParNames('N', '#mu', '#sigma_{L}', '#sigma_{R}',
                    '#alpha_{L}', 'n_{L}', '#alpha_{R}', 'n_{R}')
    
    return fun
    

def run(args):

    pars = parameters(args.mcFileName, args.mcTreeName, args.evtFileName, args.evtTreeName)

    # Create the histograms. Write them to the ROOT output file with the name
    # "mcFileName_Histos.root", where mcFileName has the .root extension removed
    if args.createHistos == 1:
        createMCHistos(pars)
        createVtxHistos(pars)

    # Plot the histograms
    plotMCHistos(pars)
    plotVtxHistos(pars)

    
def processArgs(parser):

    # Process script arguments
    parser.add_argument('--mcFileName', default='MCHierarchy.root', metavar='fileName',
                        help='MC hierarchy ROOT file [default "MCHierarchy.root"]')

    parser.add_argument('--mcTreeName', default='MC', metavar='treeName',
                        help='MC hierarchy ROOT tree [default "MC"]')

    parser.add_argument('--evtFileName', default='EventHierarchy.root', metavar='fileName',
                        help='Event hierarchy ROOT file [default "EventHierarchy.root"]')

    parser.add_argument('--evtTreeName', default='Events', metavar='treeName',
                        help='Event hierarchy ROOT tree [default "Events"]')

    parser.add_argument('--createHistos', default=1, metavar='int', type=int,
                        help='Recreate histograms [1 = Yes (default), 0 = No]')


if __name__ == '__main__':

    # Process the command line arguments
    # Use "python hierarchyPlots.py --help" to see the full list
    parser = argparse.ArgumentParser(description='List of arguments')
    processArgs(parser)
    args = parser.parse_args()

    run(args)
