''' Produces plots for tau release/data validation using the trees produced by
produceTauValTree.py
Authors: Yuta Takahashi, Michal Bluj, Jan Steggemann.
'''

import re
import warnings
from array import array
from collections import namedtuple

# The following needs to come before any other ROOT import and before argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from officialStyle import officialStyle
from variables import vardict, hvardict, cvardict
from compareTools import overlay, hoverlay, coverlay, makeEffPlotsVars, fillSampledic, findLooseId, shiftAlongX

from ROOT import gROOT, gStyle, TH1F, TH2F

import argparse
from relValTools import addArguments, dprint

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

RuntypeOptions = namedtuple("RuntypeOptions", "tlabel xlabel xlabel_eta")
options_dict = {
    'DYToLL': RuntypeOptions(tlabel='Z #rightarrow ll', xlabel='gen. lepton p_{T}^{vis} (GeV)', xlabel_eta='gen. lepton #eta^{vis}'),
    'ZTT': RuntypeOptions(tlabel='Z #rightarrow #tau#tau', xlabel='gen. tau p_{T}^{vis} (GeV)', xlabel_eta='gen. tau #eta^{vis}'),
    'ZEE': RuntypeOptions(tlabel='Z #rightarrow ee', xlabel='electron p_{T} (GeV)', xlabel_eta='electron #eta'),
    'ZMM': RuntypeOptions(tlabel='Z #rightarrow #mu#mu', xlabel='muon p_{T} (GeV)', xlabel_eta='muon #eta'),
    'QCD': RuntypeOptions(tlabel='QCD, flat #hat{p}_{T} 15-3000GeV', xlabel='jet p_{T} (GeV)', xlabel_eta='jet #eta'),
    'TTbar': RuntypeOptions(tlabel='TTbar', xlabel='jet p_{T} (GeV)', xlabel_eta='jet #eta'),
    'TTbarTau': RuntypeOptions(tlabel='TTbar #rightarrow #tau+X', xlabel='gen. tau p_{T}^{vis} (GeV)', xlabel_eta='gen. tau #eta^{vis}'),
    'TenTaus': RuntypeOptions(tlabel='Ten taus', xlabel='gen. tau p_{T}^{vis} (GeV)', xlabel_eta='gen. tau #eta^{vis}'),
    'truetauDY': RuntypeOptions(tlabel='true #tau_{h} (DY)', xlabel='gen. tau p_{T}^{vis} (GeV)', xlabel_eta='gen. tau #eta^{vis}'),
    'efakeDY': RuntypeOptions(tlabel='e #rightarrow #tau_{h} (DY)', xlabel='gen. e p_{T} (GeV)', xlabel_eta='gen. e #eta'),
    'mufakeDY': RuntypeOptions(tlabel='#mu #rightarrow #tau_{h} (DY)', xlabel='gen. #mu p_{T} (GeV)', xlabel_eta='gen. #mu #eta'),
    'jfakeDY': RuntypeOptions(tlabel='j #rightarrow #tau_{h} (DY)', xlabel='gen. jet p_{T} (GeV)', xlabel_eta='gen. jet #eta'),
    'jfakeQCD': RuntypeOptions(tlabel='j #rightarrow #tau_{h} (QCD, flat #hat{p}_{T} 15-3000GeV)', xlabel='gen. jet p_{T} (GeV)', xlabel_eta='gen. jet #eta}'),
}


def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def word_finder(expr):
    words = re.compile(r'\w+').findall(expr)
    return [w for w in words if not is_number(w) and w not in ['min', 'max']]


def efficiency_plots(d_sample, var_name, hdict):
    graphs = []
    graphs_eta = []

    for rel, rdict in sorted(d_sample.items(), key=lambda item: item[1]["index"]):
        tree = rdict['tree']
        if 'leaves' not in rdict:
            rdict['leaves'] = [leaf.GetName() for leaf in tree.GetListOfLeaves()]
        used_vars = word_finder(hdict['var'])
        if not set(used_vars).issubset(rdict['leaves']):
            with open('missing_leaves.txt', 'a+') as f:
              print >> f, var_name + ' is missing in input file ' + rdict['file'].GetName()
            warnings.warn(
                var_name + ' is missing in input file ' + rdict['file'].GetName())
            return
        num_sel = reco_cut
        den_sel = '1'
        discriminators = {"": den_sel}
        if 'against' in var_name:
            den_sel = gen_cut + ' && ' + loose_id

        rel = "tauReco @ miniAOD" if rel=="slimmedTaus_slimmedTaus" else "tauReco @ AOD"

        for mvaIDname, sel in discriminators.items():
            graphs.append(makeEffPlotsVars(tree=tree,
                                           varx='tau_genpt',
                                           numeratorAddSelection=num_sel +
                                           '&&' + hdict['var'],
                                           baseSelection=sel,
                                           binning=ptPlotsBinning,
                                           xtitle=options_dict[runtype].xlabel,
                                           header=rel + mvaIDname, addon=rel + mvaIDname,
                                           marker=rdict['marker'],
                                           col=rdict['col']))

            graphs_eta.append(makeEffPlotsVars(tree=tree,
                                               varx='tau_geneta',
                                               numeratorAddSelection=num_sel +
                                               '&&' + hdict['var'],
                                               baseSelection=sel,
                                               binning=etaPlotsBinning,
                                               xtitle=options_dict[runtype].xlabel_eta,
                                               header=rel + mvaIDname, addon=rel + mvaIDname,
                                               marker=rdict['marker'],
                                               col=rdict['col']))

    overlay(graphs=graphs,
            header=var_name,
            addon=hdict['title'],
            runtype=runtype,
            tlabel=options_dict[runtype].tlabel)

    overlay(graphs=graphs_eta,
            header=var_name + '_eta',
            addon=hdict['title'] + '_eta',
            runtype=runtype,
            tlabel=options_dict[runtype].tlabel)


def eff_plots_single(d_sample, vars_to_compare, var_dict):
    '''Adapted from Olena's code - can possibly merge it with efficiency_plots
    '''
    if not vars_to_compare:
        return

    hists = []
    histseta = []

    for index, var_name in enumerate(vars_to_compare):
        hdict = var_dict[var_name]
        if varyLooseId and 'IsolationMVA' in var_name:
            loose_id = 'tau_decayModeFindingOldDMs > 0.5 && ' + findLooseId(var_name)

        for _, rdict in sorted(d_sample.items(), key=lambda item: item[1]["index"]):
            tree = rdict['tree']
            if 'leaves' not in rdict:
                rdict['leaves'] = [leaf.GetName() for leaf in tree.GetListOfLeaves()]
            used_vars = word_finder(hdict['var'])
            if not set(used_vars).issubset(rdict['leaves']):
                warnings.warn(
                    var_name + ' is missing in input file ' + rdict['file'].GetName())
                return
            num_sel = reco_cut
            den_sel = '1'
            discriminators = {"loose_id": den_sel}
            if 'against' in var_name:
                den_sel = gen_cut + ' && ' + loose_id

            for mvaIDname, sel in discriminators.items():
                dprint("\n\tmvaIDname:", mvaIDname, "hdict['var']:", hdict['var'])

                hists.append(makeEffPlotsVars(tree=tree,
                                              varx='tau_genpt',
                                              numeratorAddSelection=num_sel + '&&' + hdict['var'],
                                              baseSelection=sel,
                                              binning=ptPlotsBinning,
                                              xtitle=options_dict[runtype].xlabel,
                                              header=var_name + mvaIDname, addon=var_name + mvaIDname,
                                              marker=rdict['marker'],
                                              col=int(colors[index])))

                shiftAlongX(hists[-1], len(vars_to_compare), index)

                histseta.append(makeEffPlotsVars(tree=tree,
                                                 varx='tau_geneta',
                                                 numeratorAddSelection=num_sel + '&&' + hdict['var'],
                                                 baseSelection=sel,
                                                 binning=etaPlotsBinning,
                                                 xtitle=options_dict[runtype].xlabel_eta,
                                                 header=var_name + mvaIDname, addon=var_name + mvaIDname,
                                                 marker=rdict['marker'],
                                                 col=int(colors[index])))

                shiftAlongX(histseta[-1], len(vars_to_compare), index)

    overlay(graphs=hists,
            header=vars_to_compare[0],
            addon=hdict['title'],
            runtype=runtype,
            tlabel=options_dict[runtype].tlabel,
            comparePerReleaseSuffix="_comparePerRelease")

    overlay(graphs=histseta,
            header=vars_to_compare[0] + '_eta',
            addon=hdict['title'] + '_eta',
            runtype=runtype,
            tlabel=options_dict[runtype].tlabel,
            comparePerReleaseSuffix="_comparePerRelease")


def var_plots(d_sample, var_name, hdict):
    hists = []
    trees = []

    for rel, rdict in sorted(d_sample.items(), key=lambda item: item[1]["index"]):

        tree = rdict['tree']
        trees.append(tree)
        if 'leaves' not in rdict:
            rdict['leaves'] = [leaf.GetName() for leaf in tree.GetListOfLeaves()]
        used_vars = word_finder(hdict['var'])
        if not set(used_vars).issubset(rdict['leaves']):
            warnings.warn(
                var_name + ' is missing in input file ' + rdict['file'].GetName())
            return
        hist = TH1F('h_' + var_name + '_' + rel, 'h_' + var_name +
                    '_' + rel, hdict['nbin'], hdict['min'], hdict['max'])

        if rel=="slimmedTaus_slimmedTaus":
            rel = "tauReco @ AOD"
        elif rel=="selectedPatTaus_selectedPatTaus":
            rel = "tauReco @ miniAOD"

        hist.GetYaxis().SetNdivisions(507)
        hist.SetLineColor(rdict['col'])
        hist.SetLineWidth(rdict['width'])
        hist.SetMinimum(0)
        hist.SetName(rel)
        hist.Sumw2()
        hist.GetXaxis().SetTitle(hdict['title'])

        # tree.Project(hist.GetName(), hdict['var'], hdict['sel'])
        #
        # if hist.Integral(0, hist.GetNbinsX() + 1) > 0:
        #     hist.Scale(1. / hist.Integral(0, hist.GetNbinsX() + 1))

        hists.append(hist)

    for i, tree in enumerate(trees):
        if args.tau_matching:
            if i == 0:
                trees[0].AddFriend(trees[1], "ft")
            elif i == 1:
                trees[1].AddFriend(trees[0], "ft")

        if additional_selection != "":
            hdict['sel'] = hdict['sel'] + '&&' + additional_selection
        # hdict['sel'] = hdict['sel'] + '&&tau_dm==0'

        if hists[i].Integral(0, hists[i].GetNbinsX() + 1) > 0:
            hists[i].Scale(1. / hists[i].Integral(0, hists[i].GetNbinsX() + 1))
        tree.Project(hists[i].GetName(), hdict['var'], hdict['sel'])

    hoverlay(hists=hists,
             xtitle=hdict['title'],
             ytitle='a.u.',
             name=var_name,
             runtype=runtype,
             tlabel=options_dict[runtype].tlabel,
             xlabel=options_dict[runtype].xlabel,
             xlabel_eta=options_dict[runtype].xlabel_eta)

def cvar_plots(d_sample, var_name, hdict):
    hists = []
    trees = []
    rels = []
    rdicts = []
    for rel, rdict in sorted(d_sample.items(), key=lambda item: item[1]["index"]):

        tree = rdict['tree']
        print rel
        print rdict['tree']
        trees.append(tree)
        if rel=="slimmedTaus_slimmedTaus":
            rel = "tauReco @ AOD"
        elif rel=="selectedPatTaus_selectedPatTaus":
            rel = "tauReco @ miniAOD"
        rels.append(rel)
        if 'leaves' not in rdict:
            rdict['leaves'] = [leaf.GetName() for leaf in tree.GetListOfLeaves()]
        rdicts.append(rdict)
        used_vars = word_finder(hdict['var'])
        if not set(used_vars).issubset(rdict['leaves']):
            warnings.warn(
                var_name + ' is missing in input file ' + rdict['file'].GetName())
            return
    if hdict['dim'] == 1:
        hist = TH1F('h_' + var_name + '_' + rels[0], 'h_' + var_name +
                    '_' + rels[0], hdict['nbin'], hdict['min'], hdict['max'])
    elif hdict['dim'] == 2:
        hist = TH2F('h_' + var_name + '_' + rels[0], 'h_' + var_name +
                    '_' + rels[0], hdict['nbin'], hdict['min'], hdict['max'], hdict['nbin'], hdict['min'], hdict['max'])


    aodtree = None
    miniaodtree = None
    if rels[0] == "tauReco @ AOD":
        aodtree = trees[0]
        miniaodtree = trees[1]
    elif rels[1] == "tauReco @ AOD":
        aodtree = trees[1]
        miniaodtree = trees[0]

    hist.GetYaxis().SetNdivisions(507)
    hist.SetLineColor(rdicts[0]['col'])
    hist.SetLineWidth(rdicts[0]['width'])
    hist.SetMinimum(0)
    hist.SetName(rels[0])
    hist.Sumw2()
    hist.GetXaxis().SetTitle(hdict['title'])


    aodtree.AddFriend(miniaodtree, "ft")
    xtitle = hdict['title']
    ytitle = 'a.u.'

    if additional_selection != "":
        hdict['sel'] = hdict['sel'] + '&&' + additional_selection
    # 'tau_dm==ft.tau_dm'
    # hdict['sel'] = hdict['sel'] + '&&tau_dm==0'
    # hdict['sel'] = hdict['sel'] + "&&(tau_hcalEnergyLeadChargedHadrCand-ft.tau_hcalEnergyLeadChargedHadrCand>0.001)"

    if hdict['dim'] == 1:
        if hdict['norm'] == 'abs':
            xtitle = 'AOD ' + hdict['title'] + ' - miniAOD ' + hdict['title']
            # aodtree.Project(hist.GetName(), hdict['var'] + "-ft." + hdict['var'], hdict['sel'])
            aodtree.Project(hist.GetName(), hdict['var'] + '-' + hdict['var'].replace('tau_', 'ft.tau_'), hdict['sel'])
        elif hdict['norm'] == 'rel':
            xtitle = '(AOD ' + hdict['title'] + ' - miniAOD ' + hdict['title'] + ')/AOD ' + hdict['title']
            aodtree.Project(hist.GetName(), '(' + hdict['var'] + '-' + hdict['var'].replace('tau_', 'ft.tau_') + ')/' + hdict['var'], hdict['sel'])
        if hist.Integral(0, hist.GetNbinsX() + 1) > 0:
            hist.Scale(1. / hist.Integral(0, hist.GetNbinsX() + 1))
    elif hdict['dim'] == 2:
        xtitle = 'AOD ' + hdict['title']
        ytitle = 'miniAOD ' + hdict['title']
        # if hdict['norm'] == 'abs':
        #     # aodtree.Project(hist.GetName(), hdict['var'] + ":ft." + hdict['var'], hdict['sel'])
        aodtree.Project(hist.GetName(), hdict['var'].replace('tau_', 'ft.tau_') + ':' + hdict['var'], hdict['sel'])


    hists.append(hist)

    coverlay(hists=hists,
             xtitle=xtitle,
             ytitle=ytitle,
             name=var_name,
             runtype=runtype,
             tlabel=options_dict[runtype].tlabel,
             xlabel=options_dict[runtype].xlabel,
             xlabel_eta=options_dict[runtype].xlabel_eta,
             sellabel=additional_selection,
             norm=hdict['norm'],
             ndim=hdict['dim'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    addArguments(parser, produce=False, compare=True)
    args = parser.parse_args()
    part = args.part
    totalparts = args.totalparts
    inputfiles = args.inputfiles

    runtype = args.runtype
    releases = args.releases
    globaltags = args.globalTags
    folders = args.folders
    # The following three are for Olena's variable comparison
    variables = args.variables
    varyLooseId = args.varyLooseId
    colors = args.colors
    additional_selection = args.selection

    sampledict = fillSampledic(
        globaltags, releases, runtype, inputfiles, folders)


    ptPlotsBinning = array('d', [20, 200]) if args.onebin else array(
        'd', [20, 30, 40, 50, 60, 70, 80, 100, 150, 200])
    etaPlotsBinning = array('d', [-2.4, 2.4]) if args.onebin else array(
        'd', [round(-2.4 + i * 0.4, 1) for i in range(13)])
    reco_cut = 'tau_pt > 20 && abs(tau_eta) < 2.3'
    gen_cut = 'tau_genpt > 20 && abs(tau_geneta) < 2.3'
    # loose_id = 'tau_decayModeFinding > 0.5 && tau_byLooseCombinedIsolationDeltaBetaCorr3Hits > 0.5'
    # loose_id = 'tau_decayModeFinding > 0.5 && tau_byLooseIsolationMVArun2v1DBoldDMwLT > 0.5'
    loose_id = '1.0'

    if part in [0, 1]:
        print "First part of plots"
        for h_name, h_dict in vardict.items():
            efficiency_plots(sampledict, h_name, h_dict)

        # Add Olena's per-release/GT plots into this script
        if variables and len(releases) == 1 and len(globaltags) == 1:
            eff_plots_single(sampledict, variables, vardict)

        print "End first part of plots"
    if part == 1:
        exit()
    elif part != 0:
        print str(part)+". part of plots"

    print "Total plots that should be made: "+str(len(hvardict.items()))
    if part == 2:
        for index, (h_name, h_dict) in enumerate(hvardict.iteritems()):
            # if index >= float(len(hvardict.items())) / (totalparts-1) * (part-1): break
            # if index < float(len(hvardict.items())) / (totalparts-1) * (part-2): continue

            if runtype not in ['ZTT', 'TTbarTau', 'TenTaus', 'truetauDY'] and h_name.find('pt_resolution') != -1:
                continue

            print "Doing",index+1, ":", h_name
            var_plots(sampledict, h_name, h_dict)

    if part == 3:
        if args.tau_matching:
            print "Total plots that should be made: "+str(len(cvardict.items()))
            for index, (c_name, c_dict) in enumerate(cvardict.iteritems()):
                # if part != 0:
                #     if index >= float(len(cvardict.items())) / (totalparts-1) * (part-1): break
                #     if index < float(len(cvardict.items())) / (totalparts-1) * (part-2): continue
                #
                # if runtype not in ['ZTT', 'TTbarTau', 'TenTaus', 'truetauDY'] and c_name.find('pt_resolution') != -1:
                #     continue

                print "Doing",index+1, ":", c_name
                cvar_plots(sampledict, c_name, c_dict)

    print "Finished"
