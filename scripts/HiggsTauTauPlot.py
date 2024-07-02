# HiggsTauTauPlot.py

import os
import sys
import argparse
from collections import OrderedDict
from prettytable import PrettyTable
import copy
from uncertainties import ufloat

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from python import Analysis  # noqa: F401
from python import root_loader  # noqa: F401
import ROOT  # noqa: F401
from python import plotting  # noqa: F401

ROOT.TH1.SetDefaultSumw2(True)

python_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python')
sys.path.append(python_dir)
parser = argparse.ArgumentParser()
parser.add_argument('--input_folder', type=str, help='Input folder')
parser.add_argument('--parameter_file', type=str, help='Parameter file')
parser.add_argument('--outputfolder', default='output', help='Output folder')
parser.add_argument('--channel', default='mm', help='Channel to run on')
parser.add_argument('--era', default='2016', help='Era to run on')
parser.add_argument('--method', default=1, help='Method to run on')
parser.add_argument('--sel', type=str, help='Additional Selection to apply')
parser.add_argument('--var', type=str, help='Variable to plot')
parser.add_argument('--do_ss', action='store_true', help='Do SS')
args = parser.parse_args()

available_channels = ['mm', 'em', 'mt', 'et', 'tt']
if args.channel not in available_channels:
    raise ValueError("Invalid channel. Please choose from: {}".format(available_channels))
available_methods = ["1"]
print(type(args.method))
if args.method not in available_methods:
    raise ValueError("Invalid method. Please choose from: {}".format(available_methods))

table = PrettyTable()
table.field_names = ['Details', 'Choices']
table.add_row(['Output Folder', args.outputfolder])
table.add_row(['Channel', args.channel])
table.add_row(['Era', args.era])
table.add_row(['Method', args.method])
table.add_row(['Selection', args.sel])
table.add_row(['Variable', args.var])

method = args.method

categories = {}
if args.era == "Run3_2022":
    if args.channel == "mm":
        categories['baseline'] = '(pt_1>25 && iso_1<0.15 && iso_2<0.15 && trg_singlemuon)'

categories['inclusive'] = '(1)'
categories['nobtag'] = '(n_bjets==0)'
categories['btag'] = '(n_bjets>=1)'
categories['w_sdb'] = 'mt_1>70.'
categories['w_shape'] = ''
qcd_os_ss_ratio = 1.07

# Samples
if args.channel == "mm":
    if args.era == "Run3_2022":
        data_samples = [
            'SingleMuon_Run2022C','Muon_Run2022C','Muon_Run2022D'
        ]
        ztt_samples = [
            'DYto2L_M-10to50_madgraphMLM','DYto2L_M-50_madgraphMLM','DYto2L_M-50_madgraphMLM_ext1','DYto2L_M-50_1J_madgraphMLM','DYto2L_M-50_2J_madgraphMLM','DYto2L_M-50_3J_madgraphMLM','DYto2L_M-50_4J_madgraphMLM'
        ]
        top_samples = [
            'TTto2L2Nu','TTto2L2Nu_ext1','TTtoLNu2Q','TTtoLNu2Q_ext1','TTto4Q','TTto4Q_ext1','TT','TT_ext1'
        ]
        vv_samples = [
            'WW','WZ','ZZ','ST_t-channel_top_4f_InclusiveDecays','ST_t-channel_antitop_4f_InclusiveDecays','ST_tW_top_2L2Nu','ST_tW_top_2L2Nu_ext1','ST_tW_antitop_2L2Nu','ST_tW_antitop_2L2Nu_ext1','ST_tW_top_4Q','ST_tW_antitop_4Q','ST_tW_antitop_4Q_ext1','ST_tW_top_LNu2Q','ST_tW_top_LNu2Q_ext1','ST_tW_antitop_LNu2Q','ST_tW_antitop_LNu2Q_ext1'
        ]
        wjets_samples = [
            'WtoLNu_madgraphMLM','WtoLNu_madgraphMLM_ext1','WtoLNu_1J_madgraphMLM','WtoLNu_2J_madgraphMLM','WtoLNu_3J_madgraphMLM','WtoLNu_4J_madgraphMLM'
        ]
        signal_samples = []

# Additional selections to separate MC samples by Gen Flags
gen_sels = {}
if args.channel == "mm":
    # 1 = prompt muon
    # 15 = muon from prompt tau
    # 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched
    gen_sels['ll_sel'] = '(genPartFlav_1==1 & genPartFlav_2==1)'
    gen_sels['tt_sel'] = '(genPartFlav_1==15 & genPartFlav_2==15)'
    gen_sels['j_sel'] = '(!(' + gen_sels['ll_sel'] + ') && !(' + gen_sels['tt_sel'] + '))'

    z_sels = {}
    z_sels['ztt_sel'] = gen_sels['tt_sel']
    z_sels['zl_sel'] = gen_sels['ll_sel']
    z_sels['zj_sel'] = gen_sels['j_sel']
    top_sels = {}
    top_sels['ttt_sel'] = gen_sels['ll_sel']
    top_sels['ttj_sel'] = '!(' + gen_sels['ll_sel'] + ')'
    vv_sels = {}
    vv_sels['vvt_sel'] = gen_sels['ll_sel']
    vv_sels['vvj_sel'] = '!(' + gen_sels['ll_sel'] + ')'

# ------------------------------------------------------------------------------------------------------------------------


def BuildCutString(wt='', sel='', cat='', sign='os',bkg_sel=''):
    full_selection = '(1)'
    if wt != '':
        full_selection = '(' + wt + ')'
    if sel != '':
        full_selection += '* (' + sel + ')'
    if sign != '':
        full_selection += '* (' + sign + ')'
    if bkg_sel != '':
        full_selection += '* (' + bkg_sel + ')'
    if cat != '':
        full_selection += '* (' + cat + ')'
    return full_selection


def GetZTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['ztt_sel'])
    return ana.SummedFactory('ZTT'+add_name, samples, plot, full_selection)

def GetZLLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zll_sel'])
    return ana.SummedFactory('ZLL'+add_name, samples, plot, full_selection)

def GetZLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zl_sel'])
    return ana.SummedFactory('ZL'+add_name, samples, plot, full_selection)

def GetZJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zj_sel'])
    return ana.SummedFactory('ZJ'+add_name, samples, plot, full_selection)

def GenerateZLL(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True, doZL=True, doZJ=True):
    if args.channel == 'em':
        zll_node = GetZLLNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
        ana.nodes[nodename].AddNode(zll_node)
    else:
        if doZL:
            zl_node = GetZLNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
            ana.nodes[nodename].AddNode(zl_node)
        if doZJ:
            zj_node = GetZJNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
            ana.nodes[nodename].AddNode(zj_node)  

def GenerateZTT(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    ztt_node = GetZTTNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)  
    ana.nodes[nodename].AddNode(ztt_node)
    split_taus = False
    if split_taus:
      ztt_node_rho = GetZTTNode(ana, '_rho'+add_name, samples, plot, wt, sel+'&&tauFlag_2==1', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_rho) 
      ztt_node_a1 = GetZTTNode(ana, '_a1'+add_name, samples, plot, wt, sel+'&&tauFlag_2==2', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_a1)
      ztt_node_pi = GetZTTNode(ana, '_pi'+add_name, samples, plot, wt, sel+'&&tauFlag_2==0', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_pi)
      ztt_node_other = GetZTTNode(ana, '_other'+add_name, samples, plot, wt, sel+'&&(tauFlag_2<0 || tauFlag_2>2)', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_other)  

def GetTTJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttj_sel'])
  return ana.SummedFactory('TTJ'+add_name, samples, plot, full_selection)


def GetTTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttt_sel'])
  return ana.SummedFactory('TTT'+add_name, samples, plot, full_selection)

def GenerateTop(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True, doTTT=True, doTTJ=True):
  wt_=wt#+"*wt_tquark_down"
  if doTTT:
      ttt_node = GetTTTNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttt_node)

  if doTTJ:
      ttj_node = GetTTJNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttj_node)

def GetVVTNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True): 
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvt_sel'])
  return ana.SummedFactory('VVT'+add_name, samples, plot, full_selection)

def GetVVJNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True): 
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvj_sel'])
  return ana.SummedFactory('VVJ'+add_name, samples, plot, full_selection)

def GenerateVV(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True, doVVT=True, doVVJ=True): 
  if doVVT:
      vvt_node = GetVVTNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvt_node)

  if doVVJ:
      vvj_node = GetVVJNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvj_node)

def GetWNode(ana, name='W', samples=[], data=[], plot='',plot_unmodified='', wt='', sel='', cat='', cat_data='', method=8, qcd_factor=qcd_os_ss_ratio, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, '')
    if categories['w_shape'] != '': shape_cat = categories['w_shape']
    else: shape_cat = cat
    shape_selection = BuildCutString(wt, sel, shape_cat, OSSS, '')
    full_selection = BuildCutString(wt, sel, cat, OSSS)
    ss_selection = BuildCutString(wt, '', cat, '!os', '')
    os_selection = BuildCutString(wt, '', cat, 'os', '')
    control_sel = categories['w_sdb']
    w_control_full_selection = BuildCutString(wt, control_sel, cat, OSSS)
    w_control_full_selection_os = BuildCutString(wt, control_sel, cat)
    w_control_full_selection_ss = BuildCutString(wt, control_sel, cat, '!os')
    w_control_full_selection_os_data = BuildCutString("weight", control_sel, cat_data)
    w_control_full_selection_ss_data = BuildCutString("weight", control_sel, cat_data, '!os')
    btag_extrap_num_node = None
    btag_extrap_den_node = None
    subtract_node_os = GetSubtractNode(ana,'_os',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,True,False) 
    subtract_node_ss = GetSubtractNode(ana,'_ss',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,False,False)

    if shape_selection == full_selection:
        w_shape = None
    else:
        w_shape = ana.SummedFactory('w_shape', samples, plot, shape_selection)
    w_node = Analysis.HttWOSSSNode(name,
      ana.SummedFactory('data_os', data, plot_unmodified, w_control_full_selection_os_data),
      subtract_node_os,
      ana.SummedFactory('data_ss', data, plot_unmodified, w_control_full_selection_ss_data),
      subtract_node_ss,
      ana.SummedFactory('W_cr', samples, plot, w_control_full_selection),
      ana.SummedFactory('W_sr', samples, plot, full_selection),
      ana.SummedFactory('W_os', samples, plot, os_selection),
      ana.SummedFactory('W_ss', samples, plot, ss_selection),
      w_shape,
      qcd_factor,
      get_os,
      btag_extrap_num_node,
      btag_extrap_den_node)

    return w_node

def GenerateW(ana, add_name='', samples=[], data=[], wg_samples=[], plot='', plot_unmodified='', wt='', sel='', cat='', cat_data='', method=12, qcd_factor=qcd_os_ss_ratio, get_os=True):
  w_node_name = 'W'
  ana.nodes[nodename].AddNode(GetWNode(ana, w_node_name+add_name, samples, data, plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_factor, get_os))

def GetSubtractNode(ana,add_name,plot,plot_unmodified,wt,sel,cat,cat_data,method,qcd_os_ss_ratio,OSSS,includeW=False,w_shift=None):
    subtract_node = Analysis.SummedNode('total_bkg'+add_name)
    if includeW:
        if w_shift is not None: w_wt = '%s*%f' %(wt,w_shift)
        else: w_wt = wt
        w_node = GetWNode(ana, 'W', wjets_samples, data_samples, plot, plot_unmodified, w_wt, sel, cat, cat_data,method, qcd_os_ss_ratio, OSSS)
        subtract_node.AddNode(w_node)
    ttt_node = GetTTTNode(ana, "", top_samples, plot, wt, sel, cat, top_sels, OSSS)
    ttj_node = GetTTJNode(ana, "", top_samples, plot, wt, sel, cat, top_sels, OSSS)
    vvt_node = GetVVTNode(ana, "", vv_samples, plot, wt, sel, cat, vv_sels, OSSS)
    vvj_node = GetVVJNode(ana, "", vv_samples, plot, wt, sel, cat, vv_sels, OSSS)
    subtract_node.AddNode(ttt_node)
    subtract_node.AddNode(ttj_node)
    subtract_node.AddNode(vvt_node)
    subtract_node.AddNode(vvj_node)

    ztt_node = GetZTTNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
    subtract_node.AddNode(ztt_node)

    zl_node = GetZLNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
    zj_node = GetZJNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
    subtract_node.AddNode(zl_node)
    subtract_node.AddNode(zj_node)

    return subtract_node

def GenerateQCD(ana, add_name='', data=[], plot='', plot_unmodified='', wt='', sel='', cat='', cat_data='', method=8, qcd_factor=qcd_os_ss_ratio, get_os=True,w_shift=None):
    shape_node = None
    OSSS = "!os"
    if get_os: OSSS = "os"

    if args.channel != 'tt':
        sub_shift='*1.0'
        if 'qcd_sub_up' in systematic: sub_shift = '*1.1'
        if 'qcd_sub_down' in systematic: sub_shift = '*0.9'

        qcd_os_ss_factor = qcd_factor
        weight = wt

        shape_node = None

        full_selection = BuildCutString("weight", sel, cat_data, '!os')
        subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight+sub_shift,sel,cat,cat_data,method,qcd_os_ss_ratio,False,True)
        if get_os: qcd_ratio = qcd_os_ss_factor
        else: qcd_ratio = 1.0
        ana.nodes[nodename].AddNode(Analysis.HttQCDNode('QCD'+add_name,
          ana.SummedFactory('data_ss', data, plot_unmodified, full_selection),
          subtract_node,
          qcd_ratio,
          shape_node))

def RunPlotting(ana, cat='',cat_data='', sel='', add_name='', wt='wt', do_data=True, samples_to_skip=[], outfile='output.root',ff_syst_weight=None):
    # produce template for observed data
    zll_samples = ztt_samples
    doZL = True
    doZJ = True
    doTTT = True
    doTTJ = True
    doVVT = True
    doVVJ = True

    if do_data:
        if args.do_ss:
          OSSS = '!os'
        else:
          OSSS = 'os'
        weight='weight'
        full_selection = BuildCutString(weight, sel, cat_data, OSSS)
        ana.nodes[nodename].AddNode(ana.SummedFactory('data_obs', data_samples, plot_unmodified, full_selection))
    GenerateZTT(ana, add_name, ztt_samples, plot, wt, sel, cat, z_sels, not args.do_ss)
    GenerateZLL(ana, add_name, zll_samples, plot, wt, sel, cat, z_sels, not args.do_ss,doZL,doZJ)
    GenerateTop(ana, add_name, top_samples, plot, wt, sel, cat, top_sels, not args.do_ss, doTTT, doTTJ)
    GenerateVV(ana, add_name, vv_samples, plot, wt, sel, cat, vv_sels, not args.do_ss, doVVT, doVVJ)
    GenerateW(ana, add_name, wjets_samples, data_samples, [], plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_os_ss_ratio, not args.do_ss)
    GenerateQCD(ana, add_name, data_samples, plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_os_ss_ratio, not args.do_ss)

def PrintSummary(nodename='', data_strings=['data_obs'], add_names=''):
    print('')
    print('################### Summary ###################')
    nodes = analysis.nodes[nodename].SubNodes()
    bkg_total = ufloat(0.000000001,0.000000001)
    sig_total = ufloat(0.000000001,0.000000001)
    for node in nodes:
        if args.channel == 'em' and node.name == 'W': continue
        if node.shape.rate.n == 0: per_err = 0
        else: per_err = node.shape.rate.s/node.shape.rate.n
        print(node.name.ljust(10) , ("%.2f" % node.shape.rate.n).ljust(10), '+/-'.ljust(5), ("%.2f" % node.shape.rate.s).ljust(7), "(%.4f)" % per_err)
        if True in [node.name.find(add_name) != -1 and add_name != '' for add_name in add_names]: continue
        if len(signal_samples) != 0: sig_samp_cond = [node.name.find(sig) != -1 for sig in signal_samples.keys()]
        else: sig_samp_cond = []
        if True in sig_samp_cond and node.name.find("_SM"+args.add_sm_background) ==-1:
            sig_total += node.shape.rate
        elif node.name not in data_strings or (args.method == 18 and 'jetFakes' == node.name):
            bkg_total += node.shape.rate
    if bkg_total.n == 0: per_err = 0
    else: per_err = bkg_total.s/bkg_total.n
    print('Total bkg'.ljust(10) , ("%.2f" % bkg_total.n).ljust(10), '+/-'.ljust(5), ("%.2f" % bkg_total.s).ljust(7), "(%.4f)" % per_err)
    if sig_total.n == 0: per_err = 0
    else: per_err = sig_total.s/sig_total.n
    print('Total sig'.ljust(10) , ("%.2f" % sig_total.n).ljust(10), '+/-'.ljust(5), ("%.2f" % sig_total.s).ljust(7), "(%.4f)" % per_err)
    print('###############################################')
    print('')

def GetTotals(ana,add_name="",outfile='outfile.root'):
    # add histograms to get totals for backgrounds split into real/fake taus and make a total backgrounds histogram
    from itertools import chain
    outfile.cd(nodename)
    nodes = ana.nodes[nodename].SubNodes()
    first_hist=True
    for node in nodes:
      if add_name not in node.name: continue
      if node.name in list(chain.from_iterable(signal_samples)): continue
      # if node.name == "data_obs": continue
      if node.name.endswith("Up"): continue
      if node.name.endswith("Down"): continue
      if node.name.endswith("extrap"): continue
      if first_hist:
          total_bkg = ana.nodes[nodename].nodes[node.name].shape.hist.Clone()
          first_hist=False
      else: total_bkg.Add(ana.nodes[nodename].nodes[node.name].shape.hist.Clone())
    if not first_hist:
      total_bkg.SetName('total_bkg'+add_name)
      total_bkg.Write()
    outfile.cd()

def FixBins(ana,outfile='output.root'):
    #Fix empty histograms
    nodes = ana.nodes[nodename].SubNodes()
    for node in nodes:
        if 'data_obs' in node.name: continue
        hist = outfile.Get(nodename+'/'+node.name)
        outfile.cd(nodename)
        #Fix empty histogram
        if hist.Integral() == 0.0:
            hist.SetBinContent(hist.GetNbinsX()//2, 0.00001)
            hist.SetBinError(hist.GetNbinsX()//2, 0.00001)
            hist.Write(hist.GetName(),ROOT.TObject.kWriteDelete)
        outfile.cd()

def Get1DBinNumFrom2D(h2d,xbin,ybin):
    Nxbins = h2d.GetNbinsX()
    return (ybin-1)*Nxbins + xbin -1

def Get1DBinNumFrom3D(h3d,xbin,ybin,zbin):
    Nxbins = h3d.GetNbinsX()
    Nybins = h3d.GetNbinsY()
    return (zbin-1)*Nxbins*Nybins + (ybin-1)*Nxbins + xbin -1

def UnrollHist2D(h2d,inc_y_of=True):
    # inc_y_of = True includes the y over-flow bins
    if inc_y_of: n = 1
    else: n = 0
    Nbins = (h2d.GetNbinsY()+n)*(h2d.GetNbinsX())
    h1d = ROOT.TH1D(h2d.GetName(), '', Nbins, 0, Nbins)
    for i in range(1,h2d.GetNbinsX()+1):
      for j in range(1,h2d.GetNbinsY()+1+n):
        glob_bin = Get1DBinNumFrom2D(h2d,i,j)
        content = h2d.GetBinContent(i,j)
        error = h2d.GetBinError(i,j)
        h1d.SetBinContent(glob_bin+1,content)
        h1d.SetBinError(glob_bin+1,error)
    return h1d

def UnrollHist3D(h3d,inc_y_of=False,inc_z_of=True):
    if inc_y_of: ny = 1
    else: ny = 0
    if inc_z_of: nz = 1
    else: nz = 0

    Nbins = (h3d.GetNbinsZ()+nz)*(h3d.GetNbinsY()+ny)*(h3d.GetNbinsX())
    h1d = ROOT.TH1D(h3d.GetName(), '', Nbins, 0, Nbins)
    for i in range(1,h3d.GetNbinsX()+1):
      for j in range(1,h3d.GetNbinsY()+1+ny):
        for k in range(1,h3d.GetNbinsZ()+1+nz):
          glob_bin = Get1DBinNumFrom3D(h3d,i,j,k)
          content = h3d.GetBinContent(i,j,k)
          error = h3d.GetBinError(i,j,k)
          h1d.SetBinContent(glob_bin+1,content)
          h1d.SetBinError(glob_bin+1,error)
    return h1d
# ----------------------------
is_2d=False
is_3d=False
var_name = args.var.split('[')[0]
var_name = var_name.split('(')[0]
var_name = var_name.replace('/','_over_')
if var_name.count(',') == 1:
    is_2d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]
if var_name.count(',') == 2:
    is_3d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]+'_vs_'+var_name.split(',')[2]

output_name = f'{args.outputfolder}/datacard_{var_name}_{args.channel}_{args.era}.root'
outfile = ROOT.TFile(output_name, 'RECREATE')
# ----------------------------

# set systematics:
# - 1st index sets folder name contaning systematic samples
# - 2nd index sets string to be appended to output histograms
# - 3rd index specifies the weight to be applied
# - 4th lists samples that should be skipped
nodename = args.channel
systematics = OrderedDict()
systematics['nominal'] = ('nominal','','(w_DY_soup*w_WJ_soup*w_Pileup*w_Muon_ID*w_Muon_Reco*w_Muon_Isolation*(w_Muon_Trigger**(1/2)))',[],False)

categories['cat'] = '('+categories['baseline']+')'
categories_unmodified = copy.deepcopy(categories)
sel = args.sel
plot = args.var
plot_unmodified = plot

samples_to_skip_dict = {}
systematic_suffixes = []

max_systematics_per_pass = 10
while len(systematics) > 0:

    analysis = Analysis.Analysis()
    analysis.nodes.AddNode(Analysis.ListNode(args.channel))
    analysis.remaps = {}

    if args.channel in ["mm"]:
        analysis.remaps['Muon'] = 'data_obs'

    previous_systematic_variation = None
    for index, systematic in enumerate(list(systematics.keys())[:max_systematics_per_pass]):
        if previous_systematic_variation is not None and systematics[systematic][0] != previous_systematic_variation: continue
        previous_systematic_variation = systematics[systematic][0]
        print("Processing:", systematic)
        print("")

        plot = args.var
        categories = copy.deepcopy(categories_unmodified)
        wshift=1.0
        systematic_folder_name = systematics[systematic][0]
        systematic_suffix = systematics[systematic][1]
        weight = systematics[systematic][2]
        samples_to_skip = systematics[systematic][3]
        is_FFsyst = systematics[systematic][4]

        systematic_suffixes.append(systematic_suffix)

        for sample_name in data_samples:
            analysis.AddSamples(f'{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root', 'ntuple', None, sample_name)

        for sample_name in ztt_samples + top_samples + vv_samples + wjets_samples:
            analysis.AddSamples(f'{args.input_folder}/{args.era}/{args.channel}/{sample_name}/{systematic_folder_name}/merged.root', 'ntuple', None, sample_name)

        analysis.AddInfo(args.parameter_file, scaleTo='data_obs')

        if systematic == 'nominal': do_data = True
        else: do_data = False

        RunPlotting(analysis, categories['cat'], categories_unmodified['cat'], sel, systematic_suffix, weight, do_data, samples_to_skip, outfile)

        del systematics[systematic]

    analysis.Run()
    analysis.nodes.Output(outfile)

    FixBins(analysis,outfile)
    for n in systematic_suffixes:
      GetTotals(analysis,n,outfile)
    PrintSummary(nodename, ['data_obs'], systematic_suffixes)

outfile.Close()

plot_file = ROOT.TFile(output_name, 'READ')
titles = plotting.SetAxisTitles(args.var,"zmm")
x_title = titles[0]
y_title = titles[1]


plotting.HTTPlot(
  nodename=nodename,
  infile=plot_file,
  channel="zmm",
  ratio_range="0.5,1.5",
  x_title=x_title,
  y_title=y_title,
  plot_name=output_name.replace('.root',''),
)