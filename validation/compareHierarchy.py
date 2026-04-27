# Compare hierarchy QA histograms made by hierarchyPlots.py
# between two different versions

import argparse
import array
import math
import os
import ROOT
import sys

class parameters(object):

    def __init__(self, hFileName1, eFileName1, hFileName2, eFileName2):
        # Hierarchy histogram files
        self.hFileName1 = hFileName1
        self.eFileName1 = eFileName1
        self.hFileName2 = hFileName2
        self.eFileName2 = eFileName2


def compareHistos(pars):

    # Open MC histogram ROOT file
    print('plotMCHistos using {0}, {1}, {2}, {3}'.format(pars.hFileName1, pars.eFileName1,
                                                         pars.hFileName2, pars.eFileName2))

    # Hierarchy histogram files (MC & event)
    hFile1 = ROOT.TFile.Open(pars.hFileName1, 'read')
    hFile2 = ROOT.TFile.Open(pars.hFileName2, 'read')
    eFile1 = ROOT.TFile.Open(pars.eFileName1, 'read')
    eFile2 = ROOT.TFile.Open(pars.eFileName2, 'read')

    # Define plotting canvas
    theCanvas = ROOT.TCanvas('theCanvas', '', 1350, 700)
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetOptStat(0)
    theCanvas.UseCurrentStyle()

    # For text labelling
    text = ROOT.TLatex()
    text.SetTextSize(0.075)
    text.SetNDC(True)

    # All interactions: muons, protons, piplus, piminus, electrons & photons
    theCanvas.Divide(3,2)

    # Hits efficiency
    maxNHits = 1.0e4

    theCanvas.cd(1)
    h1_muHitsEff = hFile1.Get('muon_HitsEff')
    h1_muHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_muHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_muHitsEff.SetLineColor(ROOT.kBlue)
    h1_muHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_muHitsEff.Draw()
    h2_muHitsEff = hFile2.Get('muon_HitsEff')
    h2_muHitsEff.SetLineColor(ROOT.kRed)
    h2_muHitsEff.SetMarkerColor(ROOT.kRed)
    h2_muHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    h1_pHitsEff = hFile1.Get('proton_HitsEff')
    h1_pHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_pHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pHitsEff.SetLineColor(ROOT.kBlue)
    h1_pHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_pHitsEff.Draw()
    h2_pHitsEff = hFile2.Get('proton_HitsEff')
    h2_pHitsEff.SetLineColor(ROOT.kRed)
    h2_pHitsEff.SetMarkerColor(ROOT.kRed)
    h2_pHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    h1_eHitsEff = hFile1.Get('electron_HitsEff')
    h1_eHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_eHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_eHitsEff.SetLineColor(ROOT.kBlue)
    h1_eHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_eHitsEff.Draw()
    h2_eHitsEff = hFile2.Get('electron_HitsEff')
    h2_eHitsEff.SetLineColor(ROOT.kRed)
    h2_eHitsEff.SetMarkerColor(ROOT.kRed)
    h2_eHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'e')

    theCanvas.cd(4)
    h1_pipHitsEff = hFile1.Get('piplus_HitsEff')
    h1_pipHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_pipHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pipHitsEff.SetLineColor(ROOT.kBlue)
    h1_pipHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_pipHitsEff.Draw()
    h2_pipHitsEff = hFile2.Get('piplus_HitsEff')
    h2_pipHitsEff.SetLineColor(ROOT.kRed)
    h2_pipHitsEff.SetMarkerColor(ROOT.kRed)
    h2_pipHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(5)
    h1_pimHitsEff = hFile1.Get('piminus_HitsEff')
    h1_pimHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_pimHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pimHitsEff.SetLineColor(ROOT.kBlue)
    h1_pimHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_pimHitsEff.Draw()
    h2_pimHitsEff = hFile2.Get('piminus_HitsEff')
    h2_pimHitsEff.SetLineColor(ROOT.kRed)
    h2_pimHitsEff.SetMarkerColor(ROOT.kRed)
    h2_pimHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.cd(6)
    h1_gHitsEff = hFile1.Get('photon_HitsEff')
    h1_gHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h1_gHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_gHitsEff.SetLineColor(ROOT.kBlue)
    h1_gHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h1_gHitsEff.Draw()
    h2_gHitsEff = hFile2.Get('photon_HitsEff')
    h2_gHitsEff.SetLineColor(ROOT.kRed)
    h2_gHitsEff.SetMarkerColor(ROOT.kRed)
    h2_gHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#gamma')

    theCanvas.Print('compare_allHitsEff.png')

    # Momentum efficiency
    maxMtm = 5.0

    theCanvas.cd(1)
    h1_muMtmEff = hFile1.Get('muon_MtmEff')
    h1_muMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_muMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_muMtmEff.SetLineColor(ROOT.kBlue)
    h1_muMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_muMtmEff.Draw()
    h2_muMtmEff = hFile2.Get('muon_MtmEff')
    h2_muMtmEff.SetLineColor(ROOT.kRed)
    h2_muMtmEff.SetMarkerColor(ROOT.kRed)
    h2_muMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    h1_pMtmEff = hFile1.Get('proton_MtmEff')
    h1_pMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_pMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pMtmEff.SetLineColor(ROOT.kBlue)
    h1_pMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_pMtmEff.Draw()
    h2_pMtmEff = hFile2.Get('proton_MtmEff')
    h2_pMtmEff.SetLineColor(ROOT.kRed)
    h2_pMtmEff.SetMarkerColor(ROOT.kRed)
    h2_pMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    h1_eMtmEff = hFile1.Get('electron_MtmEff')
    h1_eMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_eMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_eMtmEff.SetLineColor(ROOT.kBlue)
    h1_eMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_eMtmEff.Draw()
    h2_eMtmEff = hFile2.Get('electron_MtmEff')
    h2_eMtmEff.SetLineColor(ROOT.kRed)
    h2_eMtmEff.SetMarkerColor(ROOT.kRed)
    h2_eMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'e')

    theCanvas.cd(4)
    h1_pipMtmEff = hFile1.Get('piplus_MtmEff')
    h1_pipMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_pipMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pipMtmEff.SetLineColor(ROOT.kBlue)
    h1_pipMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_pipMtmEff.Draw()
    h2_pipMtmEff = hFile2.Get('piplus_MtmEff')
    h2_pipMtmEff.SetLineColor(ROOT.kRed)
    h2_pipMtmEff.SetMarkerColor(ROOT.kRed)
    h2_pipMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(5)
    h1_pimMtmEff = hFile1.Get('piminus_MtmEff')
    h1_pimMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_pimMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_pimMtmEff.SetLineColor(ROOT.kBlue)
    h1_pimMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_pimMtmEff.Draw()
    h2_pimMtmEff = hFile2.Get('piminus_MtmEff')
    h2_pimMtmEff.SetLineColor(ROOT.kRed)
    h2_pimMtmEff.SetMarkerColor(ROOT.kRed)
    h2_pimMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.cd(6)
    h1_gMtmEff = hFile1.Get('photon_MtmEff')
    h1_gMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h1_gMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_gMtmEff.SetLineColor(ROOT.kBlue)
    h1_gMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h1_gMtmEff.Draw()
    h2_gMtmEff = hFile2.Get('photon_MtmEff')
    h2_gMtmEff.SetLineColor(ROOT.kRed)
    h2_gMtmEff.SetMarkerColor(ROOT.kRed)
    h2_gMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#gamma')

    theCanvas.Print('compare_allMtmEff.png')

    # Vertexing: normalise histos
    theCanvas2 = ROOT.TCanvas('theCanvas2', '', 900, 700)
    theCanvas2.UseCurrentStyle()
    theCanvas2.Divide(2,2)

    theCanvas2.cd(1)
    h1_allVtxDX = eFile1.Get('allVtxDX')
    h1_allVtxDX.GetYaxis().SetTitle('')
    h1_allVtxDX.Scale(1.0/h1_allVtxDX.Integral())
    #h1_allVtxDX.GetXaxis().SetRangeUser(0.0, maxMtm)
    #h1_allVtxDX.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_allVtxDX.SetLineColor(ROOT.kBlue)
    h1_allVtxDX.SetMarkerColor(ROOT.kBlue)
    #ROOT.gPad.SetLogx(0)
    h1_allVtxDX.Draw()
    h2_allVtxDX = eFile2.Get('allVtxDX')
    h2_allVtxDX.GetYaxis().SetTitle('')
    h2_allVtxDX.Scale(1.0/h2_allVtxDX.Integral())
    h2_allVtxDX.SetLineColor(ROOT.kRed)
    h2_allVtxDX.SetMarkerColor(ROOT.kRed)
    h2_allVtxDX.Draw('same')
    #text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas2.cd(2)
    h1_allVtxDY = eFile1.Get('allVtxDY')
    h1_allVtxDY.Scale(1.0/h1_allVtxDY.Integral())
    h1_allVtxDY.GetYaxis().SetTitle('')
    #h1_allVtxDY.GetXaxis().SetRangeUser(0.0, maxMtm)
    #h1_allVtxDY.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_allVtxDY.SetLineColor(ROOT.kBlue)
    h1_allVtxDY.SetMarkerColor(ROOT.kBlue)
    #ROOT.gPad.SetLogx(0)
    h1_allVtxDY.Draw()
    h2_allVtxDY = eFile2.Get('allVtxDY')
    h2_allVtxDY.Scale(1.0/h2_allVtxDY.Integral())
    h2_allVtxDY.GetYaxis().SetTitle('')    
    h2_allVtxDY.SetLineColor(ROOT.kRed)
    h2_allVtxDY.SetMarkerColor(ROOT.kRed)
    h2_allVtxDY.Draw('same')
    #text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas2.cd(3)
    h1_allVtxDZ = eFile1.Get('allVtxDZ')
    h1_allVtxDZ.Scale(1.0/h1_allVtxDZ.Integral())
    h1_allVtxDZ.GetYaxis().SetTitle('')
    #h1_allVtxDZ.GetXaxis().SetRangeUser(0.0, maxMtm)
    #h1_allVtxDZ.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_allVtxDZ.SetLineColor(ROOT.kBlue)
    h1_allVtxDZ.SetMarkerColor(ROOT.kBlue)
    #ROOT.gPad.SetLogx(0)
    h1_allVtxDZ.Draw()
    h2_allVtxDZ = eFile2.Get('allVtxDZ')
    h2_allVtxDZ.Scale(1.0/h2_allVtxDZ.Integral())
    h2_allVtxDZ.GetYaxis().SetTitle('')
    h2_allVtxDZ.SetLineColor(ROOT.kRed)
    h2_allVtxDZ.SetMarkerColor(ROOT.kRed)
    h2_allVtxDZ.Draw('same')
    #text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas2.cd(4)
    h1_allVtxDR = eFile1.Get('allVtxDR')
    h1_allVtxDR.Scale(1.0/h1_allVtxDR.Integral())
    h1_allVtxDR.GetYaxis().SetTitle('')
    #h1_allVtxDR.GetXaxis().SetRangeUser(0.0, maxMtm)
    #h1_allVtxDR.GetYaxis().SetRangeUser(0.0, 1.01)
    h1_allVtxDR.SetLineColor(ROOT.kBlue)
    h1_allVtxDR.SetMarkerColor(ROOT.kBlue)
    #ROOT.gPad.SetLogx(0)
    h1_allVtxDR.Draw()
    h2_allVtxDR = eFile2.Get('allVtxDR')
    h2_allVtxDR.Scale(1.0/h2_allVtxDR.Integral())
    h2_allVtxDR.GetYaxis().SetTitle('')
    h2_allVtxDR.SetLineColor(ROOT.kRed)
    h2_allVtxDR.SetMarkerColor(ROOT.kRed)
    h2_allVtxDR.Draw('same')
    #text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas2.Print('compare_allVtx.png')
    
    # Close the histogram files
    hFile1.Close()
    hFile2.Close()
    eFile1.Close()
    eFile2.Close()


def run(args):

    pars = parameters(args.hFileName1, args.eFileName1, args.hFileName2, args.eFileName2)

    # Plot the comparisons
    compareHistos(pars)


def processArgs(parser):

    # Process script arguments
    parser.add_argument('--hFileName1', default='MCHierarchy_Histos1.root', metavar='fileName',
                        help='1st MC hierarchy histogram ROOT file [default "MCHierarchy_Histos1.root"]')

    parser.add_argument('--eFileName1', default='EventHierarchy_Histos1.root', metavar='fileName',
                        help='1st event hierarchy histogram ROOT file [default "EventHierarchy_Histos1.root"]')

    parser.add_argument('--hFileName2', default='MCHierarchy_Histos2.root', metavar='fileName',
                        help='2nd MC hierarchy histogram ROOT file [default "MCHierarchy_Histos2.root"]')

    parser.add_argument('--eFileName2', default='EventHierarchy_Histos2.root', metavar='fileName',
                        help='2nd event hierarchy histogram ROOT file [default "EventHierarchy_Histos2.root"]')


if __name__ == '__main__':

    # Process the command line arguments
    # Use "python hierarchyPlots.py --help" to see the full list
    parser = argparse.ArgumentParser(description='List of arguments')
    processArgs(parser)
    args = parser.parse_args()

    run(args)
