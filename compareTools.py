import os
import errno
import pprint

from ROOT import TH1F, TFile, TCanvas, TPad, TLegend, \
    TGraphAsymmErrors, Double, TLatex, TMath, TPaveStats, \
    gStyle, gPad

pp = pprint.PrettyPrinter(indent=4)


def ensureDir(file_path):
    if '/' not in file_path:
        return
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                print "Bad name for directory:", file_path
                raise


def save(canvas, name, extension='.png'):
    name = name.replace(' ', '').replace('&&', '')
    ensureDir(name)
    canvas.SaveAs(name + extension)


def configureLegend(leg, ncolumn):
    leg.SetNColumns(ncolumn)
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(42)


def overlay(graphs, header, addon, runtype,
            tlabel, comparePerReleaseSuffix=""):
    dir_translator = {
        "1p": "1prong",
        "1ppi0": "1prongpizero",
        "2p": "2prong",
        "3p": "3prong",
        "3ppi0": "3prongpizero",
        "oldDM": "oldDM",
        "newDM": "newDM",
        "newDMwo2p": "newDMwithout2prong",
    }

    ymin = min(TMath.MinElement(g.GetN(), g.GetY()) for g in graphs)
    ymax = max(TMath.MaxElement(g.GetN(), g.GetY()) for g in graphs)

    c = TCanvas(tlabel)

    # Upper plot will be in pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)  # Upper and lower plot are joined
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()              # pad1 becomes the current pad
    pad1.SetLogy(0)

    leg = TLegend(0.2, 0.7, 0.5, 0.9)
    configureLegend(leg, 1)

    ratio_graphs = []
    for i_graph, graph in enumerate(graphs):

        graph.GetYaxis().SetTitle('efficiency')
        graph.SetLineWidth(2)
        graph.SetMarkerSize(1)
        graph.SetMaximum(ymax * 1.4)
        graph.SetMinimum(ymin * 0.80)
        graph.GetXaxis().SetLabelSize(0.0)
        graph.GetYaxis().SetTitleSize(0.06)
        graph.GetYaxis().SetLabelSize(0.06)
        # hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin()+(3*(ii-3)),
        #                           hist.GetXaxis().GetXmax()+(3*(ii-3)))
        graph.Draw('ap' if i_graph == 0 else 'psame')

        legname = graph.GetName()
        # graph.GetPoint(0, x, y)

        if i_graph == 0:
            refg = graph.Clone()
        else:
            rp = graph.Clone()
            for i_points in range(rp.GetN()):
                if refg.GetPointY(i_points) != 0.0:
                    rp.SetPoint(i_points, graph.GetPointX(i_points), graph.GetPointY(i_points)/refg.GetPointY(i_points))
                else:
                    rp.SetPoint(i_points, graph.GetPointX(i_points), 0.0)
            ratio_graphs.append(rp)

        leg.AddEntry(graph, legname, 'lep')

    leg.Draw()

    tex = TLatex(
        (graphs[-1].GetXaxis().GetXmin() +
         0.01 * (graphs[-1].GetXaxis().GetXmax() -
                 graphs[-1].GetXaxis().GetXmin())),
        ymax * 1.4,
        addon.replace('tau_', '')
    )
    tex.SetTextAlign(10)
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.Draw()

    # xshift = 0.0
    xshift = 0.95
    # yshift = ymax*1.4
    yshift = 0.95
    if tlabel.find('QCD') != -1:
        pad1.SetTopMargin(0.1)
        # xshift = 0.0
        # yshift = ymax*1.45
        yshift = 0.95
    if runtype.find('TTbarTau') != -1:
        xshift = 0.78
    tex2 = TLatex(
        # (graphs[-1].GetXaxis().GetXmax() +
        #  xshift * (graphs[-1].GetXaxis().GetXmax() -
        #            graphs[-1].GetXaxis().GetXmin())),
        xshift,
        yshift,
        tlabel
    )
    tex2.SetNDC()
    tex2.SetTextAlign(32)
    # tex2.SetTextAlign(30)
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.Draw()


    # lower plot will be in pad
    c.cd()          # Go back to the main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.0, 1, 0.27)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.Draw()
    pad2.cd()
    for ii, ratio in enumerate(ratio_graphs):
        ratio.SetTitle("")
        ratio.GetXaxis().SetRangeUser(graphs[-1].GetXaxis().GetXmin(), graphs[-1].GetXaxis().GetXmax())
        # ratio.GetXaxis().SetTitle('#epsilon_{s}')
        ratio.GetXaxis().SetTitleOffset(0.83)
        ratio.GetXaxis().SetNdivisions(507)
        ratio.GetYaxis().SetTitle('ratio')
        ratio.GetYaxis().SetNdivisions(305)
        ratio.SetMinimum(0.82)
        ratio.SetMaximum(1.18)
        ratio.GetYaxis().SetTitleOffset(0.5)
        ratio.GetYaxis().SetTitleSize(0.16)
        ratio.GetYaxis().SetLabelSize(0.16)
        ratio.GetXaxis().SetTitleSize(0.16)
        ratio.GetXaxis().SetLabelSize(0.16)
        ratio.Draw('APsame')

    c.cd()          # Go back to the main canvas

    eta = '_eta' if '_eta' in header else ''

    dir_name = header.split('_')[0]
    save(
        c,
        'compare_' + runtype + comparePerReleaseSuffix +
        '/' + dir_name + eta + '/' + header
    )

    try:
        directory2 = dir_translator[header.split('_')[1]]
        save(
            c,
            'compare_' + runtype + comparePerReleaseSuffix +
            '/' + directory2 + eta + '/' + header
        )
    except (IndexError, KeyError):
        pass

    save(
        c,
        'compare_' + runtype +
        comparePerReleaseSuffix + '/all' +
        eta + '/' + header
    )
    pad1 = None
    pad2 = None
    c = None


def hoverlay(hists, xtitle, ytitle,
             name, runtype, tlabel,
             xlabel, xlabel_eta, comparePerReleaseSuffix=""):
    c = TCanvas(xlabel)

    # Upper plot will be in pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)  # Upper and lower plot are joined
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()              # pad1 becomes the current pad
    pad1.SetLogy(0)
    if any(subname in name for subname in ['isoPt', 'outOfConePt', 'IsoRaw', 'deepTau']):
        pad1.SetLogy()

    ymax = max([hist.GetMaximum() for hist in hists])
    leg = TLegend(0.65, 0.65, 0.85, 0.9)
    configureLegend(leg, 1)

    hratios = []
    for i_hist, hist in enumerate(hists):

        hist.SetMaximum(ymax * 1.2)
        hist.SetMinimum(0.)
        # if c.GetLogy > 0:
        if pad1.GetLogy > 0:
            hist.SetMinimum(0.001)
        hist.SetMarkerSize(0.)
        hist.GetXaxis().SetLabelSize(0.0)
        hist.GetYaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetLabelSize(0.06)
        hist.GetXaxis().SetTitle(xtitle)
        hist.GetYaxis().SetTitle(ytitle)

        if i_hist == 0:
            hist.Draw('h')
            hist_TAR = hist.Clone()
        else:
            hist.Draw('hsame')
            # ihr = hist.Clone()
            # ihr.Sumw2()
            #ihr.Divide(hists[0])

            # ihr = hist_TAR.Clone()
            # ihr.Divide(hist)
            ihr = hist.Clone()
            ihr.Divide(hist_TAR)

            ihr.SetStats(0)
            ihr.SetLineColor(hist.GetLineColor())
            ihr.SetMarkerColor(hist.GetMarkerColor())
            ihr.SetMarkerStyle(hist.GetMarkerStyle())
            hratios.append(ihr)

        leg.AddEntry(hist, hist.GetName(), "l")

    leg.Draw()

    xshift = 0.0
    yshift = ymax*1.4
    if tlabel.find('QCD') != -1:
        # pad1.SetTopMargin(0.1)
        xshift = 0.0
        # yshift = ymax*1.45
    if runtype.find('TTbarTau') != -1:
        xshift = 0.78
    tex2 = TLatex(
        0.95,
        0.95,
        tlabel
    )
    tex2.SetNDC()
    tex2.SetTextAlign(32)
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.043)
    tex2.Draw()

    # lower plot will be in pad
    c.cd()          # Go back to the main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.0, 1, 0.27)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.Draw()
    pad2.cd()
    for ii, hist in enumerate(hratios):
        hist.GetXaxis().SetTitle(xtitle)
        hist.GetXaxis().SetTitleOffset(0.83)
        hist.GetXaxis().SetNdivisions(507)
        hist.GetYaxis().SetTitle('ratio')
        hist.GetYaxis().SetNdivisions(503)
        hist.SetMinimum(0.75)
        hist.SetMaximum(1.25)
        hist.GetYaxis().SetTitleOffset(0.45)
        hist.GetYaxis().SetTitleSize(0.16)
        hist.GetYaxis().SetLabelSize(0.16)
        hist.GetXaxis().SetTitleSize(0.16)
        hist.GetXaxis().SetLabelSize(0.16)
        if ii == 0:
            hist.Draw('ep')
        else:
            hist.Draw('epsame')

    c.cd()          # Go back to the main canvas
    save(
        c,
        ('compare_' +
         runtype +
         comparePerReleaseSuffix +
         '/histograms/hist_' +
         name)
    )
    pad1 = None
    pad2 = None
    c = None

def coverlay(hists, xtitle, ytitle,
             name, runtype, tlabel,
             xlabel, xlabel_eta, sellabel, comparePerReleaseSuffix="", norm="", ndim=1):

    gStyle.SetOptStat(111110)
    c = TCanvas()
    if ndim == 1:
        c.SetRightMargin(0.10)
    elif ndim == 2:
        c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.10)
    # Upper plot will be in pad1
    # pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    # pad1.SetBottomMargin(0.03)  # Upper and lower plot are joined
    # pad1.Draw()             # Draw the upper pad: pad1
    # pad1.cd()              # pad1 becomes the current pad
    # c.SetLogy(0)
    # if any(subname in name for subname in ['isoPt', 'outOfConePt', 'IsoRaw', 'deepTau']):
    #     c.SetLogy()
    if ndim == 1:
        c.SetLogy()
    if ndim == 2:
        c.SetLogz()

    ymax = max([hist.GetMaximum() for hist in hists])
    # leg = TLegend(0.65, 0.65, 0.85, 0.9)
    # configureLegend(leg, 1)

    hratios = []
    for i_hist, hist in enumerate(hists):

        hist.SetMaximum(ymax * 1.2)
        hist.SetMinimum(0.)
        # if c.GetLogy > 0:
        if c.GetLogy > 0:
            hist.SetMinimum(0.001)
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.GetXaxis().SetTitle(xtitle)
        hist.GetYaxis().SetTitle(ytitle)

        if i_hist == 0:
            if ndim == 1:
                hist.SetMarkerSize(0.)
                hist.Draw('h')
                # gPad.Update()
                hist.SetStats(1)
                # st = hist.GetListOfFunctions().FindObject("stats")
                # # st.SetOptStat(1110)
                # # st = c.GetPrimitive("stats")
                # # st.SetNDC()
                # st.SetX1NDC(0.7)
                # st.SetX2NDC(0.85)
            elif ndim == 2:
                hist.SetMarkerSize(1.2)
                min_bincontent = 0.005
                if norm == "abs":
                    gStyle.SetPaintTextFormat("g")
                else:
                    gStyle.SetPaintTextFormat("4.2f")
                    hist.SetMaximum(100)
                    hist.SetMinimum(0.0049)
                if norm == "global":
                    integral = hist.Integral(1,hist.GetNbinsX(),1,hist.GetNbinsY())
                    for i in range(1, hist.GetNbinsX() + 1):
                        for j in range(1, hist.GetNbinsY() + 1):
                            bincontent = hist.GetBinContent(i,j)/integral*100 if integral > 0 else 0
                            bincontent = min_bincontent if (bincontent > 0 and bincontent < min_bincontent) else bincontent
                            hist.SetBinContent(i,j,bincontent)
                elif norm == "col":
                    for i in range(1, hist.GetNbinsX() + 1):
                        integral = hist.Integral(i,i,1,hist.GetNbinsY())
                        for j in range(1, hist.GetNbinsY() + 1):
                            bincontent = hist.GetBinContent(i,j)/integral*100 if integral > 0 else 0
                            bincontent = min_bincontent if (bincontent > 0 and bincontent < min_bincontent) else bincontent
                            hist.SetBinContent(i,j,bincontent)
                elif norm == "row":
                    for j in range(1, hist.GetNbinsY() + 1):
                        integral = hist.Integral(1,hist.GetNbinsX(),j,j)
                        for i in range(1, hist.GetNbinsX() + 1):
                            bincontent = hist.GetBinContent(i,j)/integral*100 if integral > 0 else 0
                            bincontent = min_bincontent if (bincontent > 0 and bincontent < min_bincontent) else bincontent
                            hist.SetBinContent(i,j,bincontent)
                hist.Draw("COLZ TEXT")
                gPad.RedrawAxis()
        #     hist_TAR = hist.Clone()
        else:
            if ndim == 1:
                hist.Draw('hsame')
            elif ndim == 2:
                hist.Draw("TEXTCOLZsame")
        #     # ihr = hist.Clone()
        #     # ihr.Sumw2()
        #     #ihr.Divide(hists[0])
        #
        #     # ihr = hist_TAR.Clone()
        #     # ihr.Divide(hist)
        #     ihr = hist.Clone()
        #     ihr.Divide(hist_TAR)
        #
        #     ihr.SetStats(0)
        #     ihr.SetLineColor(hist.GetLineColor())
        #     ihr.SetMarkerColor(hist.GetMarkerColor())
        #     ihr.SetMarkerStyle(hist.GetMarkerStyle())
        #     hratios.append(ihr)
        #
        # leg.AddEntry(hist, hist.GetName(), "l")

    # leg.Draw()
    #
    # xshift = 0.0
    # yshift = ymax*1.4
    # if tlabel.find('QCD') != -1:
    #     # pad1.SetTopMargin(0.1)
    #     xshift = 0.0
    #     # yshift = ymax*1.45
    # if runtype.find('TTbarTau') != -1:
    #     xshift = 0.78
    tex2 = TLatex(
        0.95,
        0.96,
        tlabel
    )
    tex2.SetNDC()
    tex2.SetTextAlign(32)
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.043)
    tex2.Draw()

    seltex = TLatex(
        0.02,
        0.96,
        sellabel
    )
    seltex.SetNDC()
    seltex.SetTextAlign(12)
    seltex.SetTextFont(42)
    seltex.SetTextSize(0.023)
    seltex.Draw()

    #
    # # lower plot will be in pad
    # c.cd()          # Go back to the main canvas before defining pad2
    # pad2 = TPad("pad2", "pad2", 0, 0.0, 1, 0.27)
    # pad2.SetTopMargin(0.03)
    # pad2.SetBottomMargin(0.35)
    # pad2.Draw()
    # pad2.cd()
    # for ii, hist in enumerate(hratios):
    #     hist.GetXaxis().SetTitle(xtitle)
    #     hist.GetXaxis().SetTitleOffset(0.83)
    #     hist.GetXaxis().SetNdivisions(507)
    #     hist.GetYaxis().SetTitle('ratio')
    #     hist.GetYaxis().SetNdivisions(503)
    #     hist.SetMinimum(0.75)
    #     hist.SetMaximum(1.25)
    #     hist.GetYaxis().SetTitleOffset(0.45)
    #     hist.GetYaxis().SetTitleSize(0.16)
    #     hist.GetYaxis().SetLabelSize(0.16)
    #     hist.GetXaxis().SetTitleSize(0.16)
    #     hist.GetXaxis().SetLabelSize(0.16)
    #     if ii == 0:
    #         hist.Draw('ep')
    #     else:
    #         hist.Draw('epsame')

    c.cd()          # Go back to the main canvas
    save(
        c,
        ('compare_' +
         runtype +
         comparePerReleaseSuffix +
         '/taumatching/hist_' +
         name + '_' +
         norm + 'norm_' +
         str(ndim) +'D')
    )
    c = None
    gStyle.SetOptStat(0)


def findLooseId(hname):
    looseIdDict = {
        'tau_byLooseIsolationMVArun2v1DBoldDMwLT': [
            'byLooseIsolationMVArun2v1DBoldDMwLT',
            'byMediumIsolationMVArun2v1DBoldDMwLT',
            'byTightIsolationMVArun2v1DBoldDMwLT'
        ],
        'tau_byLooseIsolationMVArun2v1PWoldDMwLT': [
            'byLooseIsolationMVArun2v1PWoldDMwLT',
            'byMediumIsolationMVArun2v1PWoldDMwLT',
            'byTightIsolationMVArun2v1PWoldDMwLT'
        ],

        'tau_byLooseIsolationMVArun2017v1DBoldDMwLT2017': [
            'byLooseIsolationMVArun2017v1DBoldDMwLT2017',
            'byMediumIsolationMVArun2017v1DBoldDMwLT2017',
            'byTightIsolationMVArun2017v1DBoldDMwLT2017'
        ],
        'tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017': [
            'byLooseIsolationMVArun2017v2DBoldDMwLT2017',
            'byMediumIsolationMVArun2017v2DBoldDMwLT2017',
            'byTightIsolationMVArun2017v2DBoldDMwLT2017'
        ],
        'tau_byLooseIsolationMVArun2v1DBoldDMwLT2016': [
            'byLooseIsolationMVArun2v1DBoldDMwLT2016',
            'byMediumIsolationMVArun2v1DBoldDMwLT2016',
            'byTightIsolationMVArun2v1DBoldDMwLT2016'
        ],
        'tau_byLooseIsolationMVArun2v1DBnewDMwLT2016': [
            'byLooseIsolationMVArun2v1DBnewDMwLT2016',
            'byMediumIsolationMVArun2v1DBnewDMwLT2016',
            'byTightIsolationMVArun2v1DBnewDMwLT2016'
        ],
        'tau_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017': [
            'byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017',
            'byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017',
            'byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017'
        ]
    }
    for key, value in looseIdDict.items():
        if hname in value:
            return key

    return 'tau_byLooseIsolationMVArun2v1DBoldDMwLT'


def shiftAlongX(tGraph, numberOfGraphs, index):
    for binNumber in xrange(tGraph.GetN()):
        x = Double(-1)
        y = Double(-1)
        tGraph.GetPoint(binNumber, x, y)
        shift = (tGraph.GetErrorXhigh(binNumber)) / (numberOfGraphs + 1)
        x = x + shift * index
        tGraph.SetPoint(binNumber, x, y)


def makeEffPlotsVars(tree,
                     varx,
                     numeratorAddSelection,
                     baseSelection,
                     binning,
                     xtitle='', header='', addon='', marker=20, col=1):

    _denomHist_ = TH1F('h_effp_' + addon,
                       'h_effp' + addon,
                       len(binning) - 1,
                       binning)
    _nominatorHist_ = TH1F('ah_effp_' + addon, 'ah_effp' + addon,
                           len(binning) - 1,
                           binning)

    tree.Draw(varx + ' >> ' + _denomHist_.GetName(), baseSelection)
    tree.Draw(varx + ' >> ' + _nominatorHist_.GetName(),
              baseSelection + ' && ' + numeratorAddSelection)

    g_eff = TGraphAsymmErrors()
    g_eff.Divide(_nominatorHist_, _denomHist_, "cl=0.683 b(1,1) mode")
    g_eff.GetXaxis().SetTitle(xtitle)
    g_eff.GetYaxis().SetTitle('efficiency')
    g_eff.GetYaxis().SetNdivisions(507)
    g_eff.SetLineWidth(3)
    g_eff.SetName(header)
    g_eff.SetMinimum(0.)
    g_eff.GetYaxis().SetTitleOffset(1.3)
    g_eff.SetMarkerStyle(marker)
    g_eff.SetMarkerSize(1)
    g_eff.SetMarkerColor(col)
    g_eff.SetLineColor(col)
    g_eff.Draw('ap')

    return g_eff


def fillSampledic(globaltags, releases, runtype, inputfiles=None, trees=None):
    sampledict = {}
    styles = [
        {'col': 1, 'marker': 26, 'width': 2},
        {'col': 2, 'marker': 25, 'width': 2},
        {'col': 4, 'marker': 21, 'width': 2},
        {'col': 2, 'marker': 21, 'width': 2},
        {'col': 7, 'marker': 24, 'width': 2},
        {'col': 41, 'marker': 20, 'width': 2},
        {'col': 6, 'marker': 22, 'width': 2},
    ]

    for index, globalTag in enumerate(globaltags):
        name = releases[index]+"_"+globalTag
        sampledict[name] = styles[index]

        jet_run_types = ['QCD', 'TTbar']
        muon_run_types = ['ZMM', 'ZpMM']
        ele_run_types = ['ZEE']
        if runtype in jet_run_types:
            runtype = runtype + "_genJets"
        if runtype in muon_run_types:
            runtype = runtype + "_genMuon"
        if runtype in ele_run_types:
            runtype = runtype + "_genEle"

        if not inputfiles:
            sampledict[name]['file'] = TFile('Myroot_{}_{}_{}.root'.format(
                releases[index],
                globalTag,
                runtype))
        else:
            sampledict[name]['file'] = TFile(inputfiles[index])
        if not trees:
            sampledict[name]['tree'] = sampledict[name]['file'].Get('per_tau')
        else:
            print trees[index]
            sampledict[name]['tree'] = sampledict[name]['file'].Get(trees[index])

        # adding the index such that we can sort the dictionary later to have correct ratio plots
        sampledict[name]['index'] = index

    return sampledict
