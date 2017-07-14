#!/usr/bin/env python                                                                                                                                               
import ROOT
from collections import OrderedDict

#cf. http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
PLOTFORMATS=[
    ('#1f78b4',20,1001),
    ('#fb9a99',20,1001),
    ('#e31a1c',20,1001),
    ('#b15928',20,1001),
    ('#33a02c',20,1001),
    ('#b2df8a',20,1001),
    ('#fdbf6f',20,1001),
    ('#ff7f00',20,1001),
    ('#ffff99',20,1001),
    ]

#ffff99 PW+HW++ EE5C
#e31a1c ERD on CR
#b2df8a ISR dn
#fdbf6f FSR up
#fb9a99 QCD based CR
#1f78b4 PW+PY8 CUETP8M2T4
#33a02c ISR up
#b15928 Gluon move CR
#ff7f00 FSR dn

COMPARISONSETS=[
    ('PW+PY8 CUETP8M2T4', [ ('nominal',         ['t#bar{t}']), 
                                ('#deltaCUET8P2MT4',['t#bar{t} UEup',     't#bar{t} UEdn']),
                                ('FSR',             ['t#bar{t} fsr up',   't#bar{t} fsr dn']),
                                ('ISR',             ['t#bar{t} isr up',   't#bar{t} isr dn']),
                                ('hdamp',           ['t#bar{t} hdamp up', 't#bar{t} hdamp dn']),
                                ('CR',              ['t#bar{t} QCDbased', 't#bar{t} ERDon', 't#bar{t} gluon move']) ] 
     ),
    ('FSR up', [ ('nominal',         ['t#bar{t} fsr up'])]),
    ('FSR dn', [ ('nominal',         ['t#bar{t} fsr dn'])]),
    ('ISR up', [ ('nominal',         ['t#bar{t} isr up'])]),
    ('ISR dn', [ ('nominal',         ['t#bar{t} isr dn'])]),
    ('QCD based CR', [ ('nominal',         ['t#bar{t} QCDbased'])]),
    ('ERD on CR',    [ ('nominal',         ['t#bar{t} ERDon'])]),
    ('Gluon move CR',[ ('nominal',         ['t#bar{t} gluon move'])]),
    ('aMC@NLO+PY8 CUETP8M2T4', [ ('nominal', ['t#bar{t} aMC@NLO']) ]),
    ('PW+HW++ EE5C'      , [ ('nominal', ['t#bar{t} Herwig++']) ]),
    ]

PLOTTINGSETS=(
    ['PW+PY8 CUETP8M2T4'],
    ['QCD based CR','ERD on CR','Gluon move CR'],
    ['ISR up','ISR dn'],
    ['FSR up','FSR dn'],
    ['PW+HW++ EE5C']
)




"""
"""
def showProfile(plotColl,sliceAxis,xtitle,ytitle,outName):

    #check ranges and prepare frames
    imax,meanMin,meanMax=0,1.0e9,-1.0e9
    for i in xrange(0,len(plotColl)):
        p=plotColl[i]["Data"]
        meanMin=min(meanMin,p.mean[0])
        if meanMax<p.mean[0]:
            meanMax=p.mean[0]
            imax=i
        
    #check how models should be plotted (order and format)
    modelIndices=OrderedDict()
    for idx in xrange(0,len(PLOTTINGSETS)):
        for key in PLOTTINGSETS[idx]:
            modelIndices[key]=(idx,len(modelIndices))


    #make graphs
    incProfiles={}
    profiles={}
    for i in xrange(0,len(plotColl)):
        for key in plotColl[i]:
            
            #filter only for needed keys
            if key!='Data' and not key in modelIndices : continue

            p=plotColl[i][key]
            mean=p.mean[0]
            statUnc,expUnc,thUnc=p.mean[1]
            totalUnc=ROOT.TMath.Sqrt(statUnc**2+expUnc**2+thUnc**2)

            #start graphs
            if not key in incProfiles:
                incProfiles[key]=[ROOT.TGraphErrors(),ROOT.TGraphErrors()]
                incProfiles[key][0].SetName("prof_%s"%key)

                ci,fill,marker=1,0,20
                if key!='Data': 
                    ci,marker,fill=PLOTFORMATS[ modelIndices[key][1] ]
                    print ci,key
                    ci=ROOT.TColor.GetColor(ci)

                #print key,ci,marker,fill,modelIndices[key][1]
                incProfiles[key][0].SetMarkerStyle(marker)
                incProfiles[key][0].SetFillStyle(fill)
                incProfiles[key][0].SetMarkerColor(ci)
                incProfiles[key][0].SetLineColor(ci)
                incProfiles[key][0].SetFillColor(ci)                
                incProfiles[key][1]=incProfiles[key][0].Clone("profexp_%s"%key)                    
                if key!='Data':
                    incProfiles[key][1].SetMarkerColor(ci+1)
                    incProfiles[key][1].SetLineColor(ci+1)
                    incProfiles[key][1].SetFillColor(ci+1)
                profiles[key]=[x.Clone('exc'+x.GetName()) for x in incProfiles[key]]
            
            #fill plots
            if i==0:          
                xcen,xwid=1,0
                if key!='Data':
                    xwid=1.0/float(len(PLOTTINGSETS))
                    xcen=0.5+(modelIndices[key][0]+0.5)*xwid
                    xwid/=2
                incProfiles[key][0].SetPoint(0,xcen,mean)
                incProfiles[key][0].SetPointError(0,xwid,totalUnc)
                incProfiles[key][1].SetPoint(0,xcen,mean)
                incProfiles[key][1].SetPointError(0,xwid,expUnc)
            else:
                xcen=sliceAxis.GetBinCenter(i)
                xwid=sliceAxis.GetBinWidth(i)*0.5
                if key=='Data' : 
                    xwid=0
                else:
                    xcen-=xwid
                    xwid=2*xwid/float(len(PLOTTINGSETS))
                    xcen+=(modelIndices[key][0]+0.5)*xwid
                    xwid/=2
                np=profiles[key][0].GetN()
                profiles[key][0].SetPoint(np,xcen,mean)
                profiles[key][0].SetPointError(np,xwid,totalUnc)
                profiles[key][1].SetPoint(np,xcen,mean)
                profiles[key][1].SetPointError(np,xwid,expUnc)

    showFormattedUEPlots(incProfiles,profiles,modelIndices,meanMax*1.5,meanMin*0.5,imax,sliceAxis,xtitle,ytitle,outName)
    incProfilesRatio=getRatiosWithRespectTo(incProfiles,'Data')
    profilesRatio=getRatiosWithRespectTo(profiles,'Data')
    showFormattedUEPlots(incProfilesRatio,profilesRatio,modelIndices,1.5,0.7,1,sliceAxis,xtitle,'K^{%s}_{MC/data}'%ytitle,outName+'_ratio')

    for key in incProfiles:
        for i in xrange(0,len(incProfiles[key])):
            incProfiles[key][i].Delete()
            profiles[key][i].Delete()
            incProfilesRatio[key][i].Delete()
            profilesRatio[key][i].Delete()


def getRatiosWithRespectTo(profiles,refKey):
    """
    Compute the ratio using a given key as reference
    """
    profilesRatio={}
    x,y=ROOT.Double(0),ROOT.Double(0)
    xref,yref=ROOT.Double(0),ROOT.Double(0)    
    for key in profiles:
        profilesRatio[key]=[gr.Clone('ratio2%s_%s'%(refKey,gr.GetName())) for gr in profiles[key]]
        for np in xrange(0,profiles[key][0].GetN()):

            profiles[key][0].GetPoint(np,x,y)
            profiles[refKey][0].GetPoint(np,xref,yref)
            ratio=-1 if yref==0 else y/yref

            ex=profiles[key][0].GetErrorX(np)
            for i in xrange(0,len(profiles[key])):
                ey=profiles[key][i].GetErrorY(np)
                eyref=profiles[refKey][i].GetErrorY(np)
                
                ratioUnc=(eyref*y)**2
                if key!=refKey:ratioUnc+=(ey*yref)**2
                ratioUnc=ROOT.TMath.Sqrt(ratioUnc)/(yref**2)

                profilesRatio[key][i].SetPoint(np,x,ratio)
                profilesRatio[key][i].SetPointError(np,ex,ratioUnc)

    return profilesRatio

def showFormattedUEPlots(incProfiles,profiles,modelIndices,meanMax,meanMin,imax,sliceAxis,xtitle,ytitle,outName):
    """
    This method dumps the formatted plots to the canvas
    """
    legAtLeft=True if sliceAxis.GetBinCenter(imax)>(sliceAxis.GetXmax()+sliceAxis.GetXmin())*0.5 else False
    frame=ROOT.TH1F('frame','frame',sliceAxis.GetNbins(),sliceAxis.GetXbins().GetArray())
    frame.GetYaxis().SetRangeUser(meanMin,meanMax)
    incframe=ROOT.TH1F('incframe','incframe',1,0.5,1.5)
    incframe.GetYaxis().SetRangeUser(meanMin,meanMax)

    #start the canvas
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetBottomMargin(0.0)

    c.cd()
    p1=ROOT.TPad('p1','p1',0,0,0.8,0.95)
    p1.SetTopMargin(0.01)
    p1.SetRightMargin(0.01)
    p1.SetLeftMargin(0.12)
    p1.SetBottomMargin(0.1)
    p1.SetGridy()
    p1.Draw()

    c.cd()
    p2=ROOT.TPad('p2','p2',0.8,0,1,0.95)
    p2.SetTopMargin(0.01)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.01)
    p2.SetBottomMargin(0.1)
    p2.SetGridy()
    p2.Draw()

    #inclusive frame
    p2.cd()
    incframe.Draw()
    incframe.GetXaxis().SetNdivisions(0)
    incframe.GetXaxis().SetTitleSize(0.15)
    incframe.GetXaxis().SetTitle('Inclusive')
    incframe.GetXaxis().SetTitleOffset(0.25)
    for key in incProfiles:
        if key =='Data': continue
        incProfiles[key][0].Draw('2')
        #incProfiles[key][1].Draw('2')
    incProfiles['Data'][0].Draw('ep')
    #incProfiles['Data'][1].Draw('ep')

    #legend for uncertainties
    #lines=[]
    #for i in [0]: #,1]:
    #    lines.append(ROOT.TLine())
    #    lines[-1].SetLineWidth(5)
    #    dx,ci=0.25,PLOTFORMATS[0][0]
    #    if ci==1:
    #        dx=0.15
    #        ci+=1
    #    lines[-1].SetLineColor(ci)
    #    lines[-1].DrawLineNDC(0.5-dx,0.93,0.5+dx,0.93)
    #tex=ROOT.TLatex()
    #tex.SetTextFont(42)
    #tex.SetTextSize(0.1)
    #tex.SetNDC()
    #tex.DrawLatex(0.1,0.95,'#it{Total Exp.}')

    p2.RedrawAxis()
                
    #differential frame
    p1.cd()

    frame.Draw()
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetXaxis().SetLabelSize(0.04)
    for key in profiles:
        if key =='Data': continue
        profiles[key][0].Draw('2')
        #profiles[key][1].Draw('2')
    profiles['Data'][0].Draw('ep')
    #profiles['Data'][1].Draw('ep')
        
    #draw the legend
    startLegX=0.15 if legAtLeft else 0.50
    startLegY=0.9  if legAtLeft else 0.95
    leg=ROOT.TLegend(startLegX,startLegY,startLegX+0.4,startLegY-len(PLOTTINGSETS)*0.06) 
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry( incProfiles['Data'][0], 'Data','ep' )
    for key in modelIndices: leg.AddEntry(profiles[key][0],key,'f')
    leg.Draw()
        
    #the header
    tex2=ROOT.TLatex()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.045)
    tex2.SetNDC()
    tex2.DrawLatex(0.16,0.93,'#bf{CMS} #it{preliminary}')

    p1.RedrawAxis()

    #the lumi/sqrts
    c.cd()
    tex3=ROOT.TLatex()
    tex3.SetTextFont(42)
    tex3.SetTextSize(0.04)
    tex3.SetNDC()
    tex3.DrawLatex(0.7,0.96,'#scale[0.8]{35.9 fb^{-1} (#sqrt{s}=13 TeV)}')

    # all done
    c.Modified()
    c.Update()
    for ext in ['pdf','png']: c.SaveAs('%s.%s'%(outName,ext))

    p1.Delete()
    p2.Delete()
    c.Delete()
