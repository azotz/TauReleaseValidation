import ROOT
from ROOT import gROOT, gStyle, TEfficiency, TRatioPlot, TPad, TLine
from officialStyle import officialStyle

gROOT.SetBatch(True)
officialStyle(gStyle)

colours = [1, 2, 3, 4, 6, 7, 8, 9, 47, 46, 44, 43, 42, 41, 40]
markers = [20, 21, 22, 23, 24, 25, 26, 27]


def histsToRoc(hsig, hbg, w_error=False):
    '''Produce ROC curve from 2 input histograms.
    Partly adapted from Giovanni's ttH code.
    '''
    nbins = hsig.GetNbinsX() + 2 - 1 # include under/overflow; remove events not passing selection
    si = [hsig.GetBinContent(i) for i in xrange(nbins+1)]
    bi = [hbg.GetBinContent(i) for i in xrange(nbins+1)]
    del si[1]
    del bi[1]

    if hsig.GetMean() > hbg.GetMean():
        si.reverse()
        bi.reverse()

    sums, sumb = sum(si), sum(bi)
    if sums == 0 or sumb == 0:
        print 'WARNING: Either signal or background histogram empty', sums, sumb
        return None

    for i in xrange(1, nbins):
        si[i] += si[i - 1]
        bi[i] += bi[i - 1]
    fullsi, fullbi = si[:], bi[:]
    si, bi = [], []
    for i in xrange(1, nbins):
        # skip negative weights
        if si and (fullsi[i] < si[-1] or fullbi[i] < bi[-1]):
            continue
        # skip repetitions
        if fullsi[i] != fullsi[i - 1] or fullbi[i] != fullbi[i - 1]:
            si.append(fullsi[i])
            bi.append(fullbi[i])

    if len(si) == 2:
        si = [si[0]]
        bi = [bi[0]]

    bins = len(si)

    if not w_error:
        roc = ROOT.TGraph(bins)
        for i in xrange(bins):
            roc.SetPoint(i, si[i] / sums, bi[i] / sumb)

        return roc

    roc = ROOT.TGraphAsymmErrors(bins)
    for i in xrange(bins):
        interval = 0.683

        e_s_low = si[i] / sums - TEfficiency.ClopperPearson(sums, si[i], interval, False)
        e_s_up = TEfficiency.ClopperPearson(sums, si[i], interval, True) - si[i] / sums
        e_b_low = bi[i] / sumb - TEfficiency.ClopperPearson(sumb, bi[i], interval, False)
        e_b_up = TEfficiency.ClopperPearson(sumb, bi[i], interval, True) - bi[i] / sumb

        roc.SetPoint(i, si[i] / sums, bi[i] / sumb)
        roc.SetPointError(i, e_s_low, e_s_up, e_b_low, e_b_up)

    return roc


def makeLegend(rocs, textSize=0.035, left=True):
    (x1, y1, x2, y2) = (.18 if left else .68, .76 - textSize * max(len(rocs) - 3, 0), .45 if left else .95, .88)
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetFillColor(0)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetLineWidth(0)
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for key, roc in rocs:
        leg.AddEntry(roc, key, 'lp')
    leg.Draw()

    return leg


def makeROCPlot(rocs, set_name, ymin=0., ymax=1., xmin=0., xmax=1., logy=False):
    print "makeROCPlot"
    allrocs = ROOT.TMultiGraph(set_name, '')
    point_graphs = []
    ratio_graphs = []
    i_marker = 0
    c = ROOT.TCanvas()
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)  # Upper and lower plot are joined
    pad1.Draw()             # Draw the upper pad: pad1
    pad1.cd()              # pad1 becomes the current pad
    pad1.SetLogy(0)
    # pad1.SetLogx()
    if ymin > 0. and logy:
        pad1.SetLogy()

    for i_col, graph in enumerate(rocs):
        col = colours[i_col]
        graph.SetLineColor(col)
        graph.SetMarkerColor(col)
        graph.SetLineWidth(3)
        graph.SetMarkerStyle(9)
        graph.SetMarkerSize(0)
        if i_col == 0:
            refg = graph.Clone()
        else:
            n_points = refg.GetN() if refg.GetN() > graph.GetN() else graph.GetN()
            # rp = ROOT.TGraphAsymmErrors(n_points)
            rp = graph.Clone()
            rp.Set(n_points)
            # n_points = n_points - 1
            for i_points in range(n_points):
                x_refg = refg.GetPointX(i_points)
                print refg.GetPointX(i_points), refg.GetPointY(i_points), "|", graph.GetPointX(i_points), graph.Eval(x_refg)
                if refg.GetPointY(i_points) != 0.0:
                    rp.SetPoint(i_points, x_refg, graph.Eval(x_refg)/refg.GetPointY(i_points))
                else:
                    rp.SetPoint(i_points, x_refg, 0.0)
            ratio_graphs.append(rp)
        if graph.GetN() > 10:
            allrocs.Add(graph)
        else:
            graph.SetMarkerStyle(markers[i_marker])
            i_marker += 1
            graph.SetMarkerSize(1)
            point_graphs.append(graph)


    # allrocs.Draw('APL')
    # allrocs.Draw('AL')

    allrocs.GetXaxis().SetTitle('#epsilon_{s}')
    allrocs.GetYaxis().SetTitle('#epsilon_{b}')
    allrocs.GetYaxis().SetDecimals(True)

    allrocs.GetXaxis().SetLabelSize(0.0)
    allrocs.GetYaxis().SetTitleSize(0.06)
    allrocs.GetYaxis().SetLabelSize(0.06)

    allrocs.GetYaxis().SetRangeUser(ymin, ymax)
    allrocs.GetXaxis().SetRangeUser(xmin, xmax)


    allrocs.Draw('APL')
    # allrocs.Draw('AL')

    for graph in point_graphs:
        graph.Draw('Psame')

    allrocs.leg = makeLegend(zip([r.title for r in rocs], rocs))

    c.cd()          # Go back to the main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.0, 1, 0.27)
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    # pad2.SetLogx()
    pad2.Draw()
    pad2.cd()
    # axis_graph = rocs[0].Clone()
    # axis_graph.Set(2)
    # axis_graph = ROOT.TGraph(2)
    # axis_graph.SetPoint(0, allrocs.GetXaxis().GetXmin(), 1.0); axis_graph.SetPoint(1, 0.999, 1.0);
    # axis_graph.SetTitle("")
    # axis_graph.SetMarkerSize(0)
    # axis_graph.SetLineWidth(0)
    # # axis_graph.GetYaxis().SetRangeUser(0.75,1.25)
    # axis_graph.SetMinimum(0.75)
    # axis_graph.SetMaximum(1.25)
    # axis_graph.GetXaxis().SetRangeUser(allrocs.GetXaxis().GetXmin(), allrocs.GetXaxis().GetXmax())
    # pad2.Update()
    # pad2.RedrawAxis()
    # axis_graph.Draw("APL")
    # line = TLine(allrocs.GetXaxis().GetXmin(),1.,allrocs.GetXaxis().GetXmin(),1.)
    # line.SetLineColor(1)
    # line.Draw()
    # print pad2.GetUxmin(),pad2.GetUxmax()
    for ii, ratio in enumerate(ratio_graphs):
        ratio.SetTitle("")
        print allrocs.GetXaxis().GetXmin(), allrocs.GetXaxis().GetXmax()
        ratio.GetXaxis().SetRangeUser(allrocs.GetXaxis().GetXmin(), allrocs.GetXaxis().GetXmax())
        # ratio.GetXaxis().SetRangeUser(xmin, allrocs.GetXaxis().GetXmax())
        ratio.GetXaxis().SetTitle('#epsilon_{s}')
        ratio.GetXaxis().SetTitleOffset(0.83)
        ratio.GetXaxis().SetNdivisions(507)
        ratio.GetYaxis().SetTitle('ratio')
        ratio.GetYaxis().SetNdivisions(305)
        ratio.SetMinimum(0.75)
        ratio.SetMaximum(1.25)
        ratio.GetYaxis().SetTitleOffset(0.5)
        ratio.GetYaxis().SetTitleSize(0.16)
        ratio.GetYaxis().SetLabelSize(0.16)
        ratio.GetXaxis().SetTitleSize(0.16)
        ratio.GetXaxis().SetLabelSize(0.16)
        ratio.Draw('APL')
        # line.Draw("")
        # ratio.Draw('APLsame')
    pad2.Update()
    pad2.RedrawAxis()
    # c.cd()
    # c.Update()
    c.Print(set_name + '.png')

    return allrocs
