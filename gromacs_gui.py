#!/usr/bin/env python
import wx
import os
import re
import subprocess
import numpy
from functools import partial
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as NavigationToolbar

solvent_model=""
base=""
dir=""
version=5.02
#version=4.5
def findThisProcess( process_name ):
    ps     = subprocess.Popen("ps -eaf | grep "+process_name, shell=True, stdout=subprocess.PIPE)
    output = ps.stdout.read()
    ps.stdout.close()
    ps.wait()

    return output

# This is the function you can use  
def isThisRunning( process_name ):
    output = findThisProcess( process_name )
    
    if re.search(process_name, output) is None:
        return False
    else:
        return True

class MyMenuBar(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(900, 775))
        #Menubar
        menubar = wx.MenuBar()

        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_NEW, '&New')
        fileMenu.Append(wx.ID_OPEN, '&Open')
        fileMenu.Append(wx.ID_SAVE, '&Save')
        fileMenu.Append(wx.ID_EXIT, '&Quit')
        fileMenu.AppendSeparator()

        helpMenu=wx.Menu()
        helpMenu.Append(wx.ID_ANY, '&Bevan Lab Tutorial Web Page')
        about_menu_item=helpMenu.Append(wx.ID_ANY, '&About')

#        qmi = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
#        fileMenu.AppendItem(qmi)


        menubar.Append(fileMenu, '&File')
        menubar.Append(helpMenu, '&Help')
        
        #Menubar events
#        self.Bind(wx.EVT_MENU, self.OnOpen, fileMenu)
#        self.Bind(wx.EVT_MENU, self.OnExit, fileMenu)
        self.Bind(wx.EVT_MENU, self.OnAboutDlg, about_menu_item)
        self.SetMenuBar(menubar)

    def OnAboutDlg(self, event):
        info = wx.AboutDialogInfo()
        info.Name = "Quick and Dirty Gromacs Tutorial"
        info.Version = "0.0.1 Beta"
        info.Copyright = "(C) 2014 Mustafa Tekpinar\nEmail: tekpinar@buffalo.edu\nLicence: LGPL"
        wx.AboutBox(info)

class MyMainWindow(wx.Frame):

    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(900, 775))

        #Box parameters
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(vbox)

        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.plottoolbar = NavigationToolbar(self.canvas)

        self.button1 = wx.Button(self, -1, "Plot Total Energy vs. Time")
        self.button1.Bind(wx.EVT_BUTTON, self.changePlot)

        self.button2 = wx.Button(self, -1, "Plot Temperature vs. Time")
#        self.button2.Bind(wx.EVT_BUTTON, self.drawTemperaturevsTime())

        self.button3 = wx.Button(self, -1, "Plot Pressure vs. Time")
#        self.button3.Bind(wx.EVT_BUTTON, self.drawPressurevsTime())

        #Toolbar items
        toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        toolbar.AddSimpleTool(2, wx.Image('./images/stock_open.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Open', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(7, wx.Image('./images/prepare.png', wx.BITMAP_TYPE_PNG).Rescale(75, 40).ConvertToBitmap(), 'Prepare Model', '')
        toolbar.AddSimpleTool(9, wx.Image('./images/solvate.png', wx.BITMAP_TYPE_PNG).Rescale(75, 40).ConvertToBitmap(), 'Solvate', '')
        toolbar.AddSimpleTool(10, wx.Image('./images/ionize.png', wx.BITMAP_TYPE_PNG).Rescale(75, 40).ConvertToBitmap(), 'Ionize', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(11, wx.Image('./images/minimize.png', wx.BITMAP_TYPE_PNG).Rescale(90, 40).ConvertToBitmap(), 'Minimize', '')
        toolbar.AddSimpleTool(12, wx.Image('./images/equilibrate_phase1.png', wx.BITMAP_TYPE_PNG).Rescale(90, 40).ConvertToBitmap(), 'Equilibrate-Phase 1', '')
        toolbar.AddSimpleTool(13, wx.Image('./images/equilibrate_phase2.png', wx.BITMAP_TYPE_PNG).Rescale(90, 40).ConvertToBitmap(), 'Equilibrate-Phase 2', '')
        toolbar.AddSimpleTool(14, wx.Image('./images/run.png', wx.BITMAP_TYPE_PNG).Rescale(90, 40).ConvertToBitmap(), 'Production Run', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(4, wx.Image('./images/stock_exit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Exit', '')
        toolbar.Realize()

        self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer2.Add(self.button1, -1, wx.EXPAND)
        self.sizer2.Add(self.button2, -1, wx.EXPAND)
        self.sizer2.Add(self.button3, -1, wx.EXPAND)

        vbox.Add(toolbar, 0, wx.EXPAND)
        vbox.Add(self.plottoolbar, 0, wx.EXPAND)
        vbox.Add(self.sizer2, 0, wx.EXPAND)
        vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)

#        self.drawLines()
        self.current_draw = 'lines'

        self.statusbar = self.CreateStatusBar()
        self.Centre()

        #Toolbar events
        self.Bind(wx.EVT_TOOL, self.OnNew, id=1)
#        self.Bind(wx.EVT_TOOL, self.OnOpen, id=2)
        self.Bind(wx.EVT_TOOL, self.openfile, id=2)
        self.Bind(wx.EVT_TOOL, self.OnSave, id=3)
        self.Bind(wx.EVT_TOOL, self.OnExit, id=4)
#        self.Bind(wx.EVT_TOOL, lambda event: self.OnPrepare(event, base), id=7) 
#        self.Bind(wx.EVT_TOOL, lambda event: self.OnPrepare(event, base)) 
#        self.Bind(wx.EVT_TOOL, partial( self.OnPrepare, base), id=7)
        self.Bind(wx.EVT_TOOL, self.OnPrepare, id=7)
        self.Bind(wx.EVT_TOOL, self.solvate, id=9)
        self.Bind(wx.EVT_TOOL, self.ionize, id=10)
        self.Bind(wx.EVT_TOOL, self.OnMinimize, id=11)
        self.Bind(wx.EVT_TOOL, self.OnEquilibratePhase1, id=12)
        self.Bind(wx.EVT_TOOL, self.OnEquilibratePhase2, id=13)
#        self.Bind(wx.EVT_TOOL, self.tempentry, id=12)

        #Menubar
        myNewMenuBar=MyMenuBar(self, -1, "test")
        self.SetSize((800, 600))

    def changePlot(self, event):
        if self.current_draw == 'sin' :
            self.drawLines()
            self.current_draw = 'lines'
        else: 
            self.drawPotential(event)
            self.current_draw = 'sin'

        self.Layout()

    def drawLines(self):
        x = numpy.arange(0.0,1000,100.0)
        y = [0,0,0,0,0,0,0,0,0,0]

        self.axes.clear()
        self.axes.plot(x, y)

    def drawPotential(self, event):
        data_file=open(dir+"/potential.xvg", "r")
        lines=data_file.readlines()
        md_columns=[]
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
#                print line
                md_columns=line.split()
                md_pc1.append(float(md_columns[0]))
                md_pc2.append(float(md_columns[1]))
                counter+=1
        self.axes.clear()
        self.axes.plot(md_pc1, md_pc2)
 #       self.x_max=1000.0
        self.axes.set_xlabel("Step") 
        self.axes.set_ylabel("Energy(kJ/mol)") 
        data_file.close()
        self.Layout()

        
    def drawTemperaturevsTime(self, event):
        data_file=open(dir+"/temperature.xvg", "r")
        lines=data_file.readlines()
        md_columns=[]
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
#                print line
                md_columns=line.split()
                md_pc1.append(float(md_columns[0]))
                md_pc2.append(float(md_columns[1]))
                counter+=1
        self.axes.clear()
        self.axes.plot(md_pc1, md_pc2)
#        self.x_max=1000.0
        self.axes.set_xlabel("Time (ps)") 
        self.axes.set_ylabel("Temperature (K)") 
        data_file.close()
        self.Layout()

    def drawPressurevsTime(self, event):
        data_file=open(dir+"/pressure.xvg", "r")
        lines=data_file.readlines()
        md_columns=[]
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
#                print line
                md_columns=line.split()
                md_pc1.append(float(md_columns[0]))
                md_pc2.append(float(md_columns[1]))
                counter+=1
        self.axes.clear()
        self.axes.plot(md_pc1, md_pc2)
#        self.x_max=1000.0
        self.axes.set_xlabel("Time (ps)") 
        self.axes.set_ylabel("Pressure (K)") 
        data_file.close()
        self.Layout()

    def OnNew(self, event):
        self.statusbar.SetStatusText('New Command')

    def OnOpen(self, event):
        wildcard = 'Protein Data Bank Files (*.pdb)|*.pdb|'\
            'All files (*.pdb1)|*.pdb1'
        dlg_fo = wx.FileDialog(None, 'Select input file', os.getcwd(), '', wildcard, wx.OPEN)
        if dlg_fo.ShowModal() == wx.ID_OK:
            print dlg_fo.GetPath()

    def OnSave(self, event):
        self.statusbar.SetStatusText('Save Command')
    
    def OnExit(self, event):
        self.Close()

    def openfile(self, event):
       dlg = wx.FileDialog(self, "Choose an initial pdb file", os.getcwd(), "", "*.pdb", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           path = dlg.GetPath()
           global base
           global dir
           base = os.path.basename(path)
           dir  = os.path.dirname(path)
           myFileDataList=[dir, base]
           #                self.SetStatusText("You selected: %s" % base)
           #                self.SetStatusText("You selected: %s" % dir)
           self.SetStatusText("You selected: %s" % path)
       dlg.Destroy()
       return dlg.GetFilename()

    def waterchoice(self, event):
        wm = ['SPC', 'SPCE', 'TIP3P', 'TIP4P']
        dlg = wx.SingleChoiceDialog(self, 'Select water model', 'Which one?', wm, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You chose: %s\n' % dlg.GetStringSelection())
        dlg.Destroy()
        return dlg.GetStringSelection().lower()

    def ffchoice(self, event):
        ffm = ['amber03            -AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)',\
                   'amber94        -AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)',\
                    'amber96       -AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)',\
                    'amber99       -AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)',\
                    'amber99sb     -AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)',\
                    'amber99sb-ildn-AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)',\
                    'amberGS       -AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)',\
                    'charmm27      -CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)',\
                    'gromos43a1    -GROMOS96 43a1 force field',\
                    'gromos43a2    -GROMOS96 43a2 force field (improved alkane dihedrals)',\
                    'gromos45a3    -GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)',\
                    'gromos53a5    -GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)',\
                    'gromos53a6    -GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)',\
                    'gromos54a7    -GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)',\
                    'oplsaa        -OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)']
        dlg = wx.SingleChoiceDialog(self, 'Select water model', 'Which one?', ffm, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You chose: %s\n' % dlg.GetStringSelection())
        dlg.Destroy()
        return dlg.GetStringSelection()


    def OnPrepare(self, event):
#        initial_pdb=self.openfile(event)
         #This / makes it very unix dependent!!!
        initial_pdb=dir+"/"+base
#        print dir
        global solvent_model 
        if(initial_pdb != ""):
            solvent_model=self.waterchoice(event)
            if(solvent_model!=""):
                force_field=self.ffchoice(event)

                #This / makes it very unix dependent!!!
#                os.system("cd "+dir)
                prepare_command="pdb2gmx -f "+initial_pdb+" -o "+dir+"/"+"processed.pdb -i "+dir+"/"+"posre.itp -p "+dir+"/"+"topol.top -water "\
                    +solvent_model.lower()+" -ff "+force_field[0:14]
                if(version<5.0):
                    os.system(prepare_command)
                elif(version>=5):
                    os.system("gmx "+prepare_command)

    def boxtypechoice(self, event):
        bt = ['Cubic', 'Triclinic', 'Dodecahedron']
        dlg = wx.SingleChoiceDialog(self, 'Select box type', 'Which one?', bt, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You\'ve chosen %s\n' % dlg.GetStringSelection())
        dlg.Destroy()
        return dlg.GetStringSelection().lower()

    def solvate(self, event):
        #Define box properties
        boxtype=self.boxtypechoice(event)
        error_message="ERROR: Hey dude! Something went wrong in solvation precedure! Check log files!"
        presolvate_command="editconf -f "+dir+"/"+"processed.pdb -o "+dir+"/"+"newbox.pdb -bt "+boxtype+" -c -d 1.0"
        if(boxtype != ""):
            if(version<5.0):
                status=os.system(presolvate_command)
                if(status==0):
                    #Produce box
                    if ( solvent_model.lower() == "tip3p"  or  solvent_model == "spc"  or  solvent_model == "spce"):
                        os.system("genbox -cp "+dir+"/"+"newbox.pdb -cs spc216.gro -o "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top")
                    
                    elif (solvent_model == "tip4p"):
                        os.system("genbox -cp "+dir+"/"+"newbox.pdb -cs tip4p.gro -o "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top")
                else:
                    print (error_message)
            elif(version>=5.0):
                status=os.system("gmx "+presolvate_command)
                if(status==0):
                    #Produce box
                    if ( solvent_model.lower() == "tip3p"  or  solvent_model == "spc"  or  solvent_model == "spce"):
                        os.system("gmx solvate -cp "+dir+"/"+"newbox.pdb -cs spc216.gro -o "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top")
                    
                    elif (solvent_model == "tip4p"):
                        os.system("gmx solvate -cp "+dir+"/"+"newbox.pdb -cs tip4p.gro -o "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top")
                else:
                    print (error_message)

    def pos_iontypechoice(self, event):
        lst = ["Na", "K", "Ca", "Mg"]
        dlg = wx.SingleChoiceDialog(self, 'Select positive ion type', 'Which one?',lst, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You\'ve chosen %s\n' % dlg.GetStringSelection())
        dlg.Destroy()
        return dlg.GetStringSelection().upper()

    def neg_iontypechoice(self, event):
        lst = ["Cl"]
        dlg = wx.SingleChoiceDialog(self, 'Select negative ion type', 'Which one?',lst, wx.CHOICEDLG_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You\'ve chosen %s\n' % dlg.GetStringSelection())
        dlg.Destroy()
        return dlg.GetStringSelection().upper()

    def ionize(self, event):
        pos_ion=self.pos_iontypechoice(event)
        neg_ion=self.neg_iontypechoice(event)
        #Prepare for ion addition
        preprocess_command="-c "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top -o "+dir+"/"+"ions.tpr"
        ionize_command="genion -s "+dir+"/"+"ions.tpr -o "+dir+"/"+"solv_ions.pdb -p "+dir+"/"+"topol.top -pname "+\
            pos_ion+" -nname "+neg_ion+" -neutral -conc 0.15"
        error_message="ERROR: Something went wrong in ionization precedure! Check log files!"
        if(version<5.0):
            status=os.system("grompp -f ./scripts/gromacs_less_than_5/ions.mdp "+preprocess_command)
            if(status==0):
                #Produce and place ions
                os.system(ionize_command)
            else:
                print (error_message)

        elif(version>=5.0):
            status=os.system("gmx grompp -f ./scripts/gromacs_5_or_more/ions.mdp "+preprocess_command)
            if(status==0):
                #Produce and place ions
                os.system("gmx "+ionize_command)
            else:
                print (error_message)

    def tempentry(self, event):
        dlg = wx.TextEntryDialog(self, 'Enter temperature in Kelvin','Text Entry')
        dlg.SetValue("300")
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You entered: %s\n' % dlg.GetValue())
        dlg.Destroy()

    def OnMinimize(self, event):
        pre_min_command=" -c "+dir+"/"+"solv_ions.pdb -p "+dir+"/"+"topol.top -o "+dir+"/"+"em.tpr"
        if(version<5.0):
            pre_min_command="grompp -f ./scripts/gromacs_less_than_5/minim.mdp"+pre_min_command
        elif(version>=5.0):
            pre_min_command="gmx grompp -f ./scripts/gromacs_5_or_more/minim.mdp"+pre_min_command

        status1=os.system(pre_min_command)
        if(status1 == 0):
            status2=os.system("mdrun -v -deffnm "+dir+"/"+"em 1>"+dir+"/"+"em.out 2>"+dir+"/"+"em.err")
            if(status2==0):
                status2=os.system("g_energy -f "+dir+"/"+"em.edr -o "+dir+"/"+"potential.xvg")
#                    self.changePlot(event)
                self.drawPotential(event)
            self.SetStatusText('Completed minimization succesfully!\n')
        else:
            print "ERROR: Hey dude! I think something is really wrong here in minimization!"                

    def OnEquilibratePhase1(self, event):
        pre_phase1_command=" -c "+dir+"/"+"em.gro -p "+dir+"/"+"topol.top -o "+dir+"/"+"nvt.tpr"
        if(version<5.0):
            pre_phase1_command="grompp -f ./scripts/gromacs_less_than_5/nvt.mdp"+pre_phase1_command
        elif(version>=5.0):
            pre_phase1_command="gmx grompp -f ./scripts/gromacs_5_or_more/nvt.mdp"+pre_phase1_command

        status1=os.system(pre_phase1_command)
        if(status1 == 0):
            status2=os.system("mdrun -v -deffnm "+dir+"/"+"nvt")
            if(status2 == 0):
                os.system("gmx energy -f "+dir+"/"+"nvt.edr -o "+dir+"/"+"temperature.xvg")
                self.drawTemperaturevsTime(event)
                self.SetStatusText('Completed NVT equilibration succesfully!\n')
        else:
            print "ERROR: Hey dude! I think something went wrong in NVT simulation!"

#        self.SetTitle('Quick and Dirty Gromacs')
#        self.Centre()
#        self.Show(True)

    def OnEquilibratePhase2(self, event):
        pre_phase2_command=" -c "+dir+"/"+"nvt.gro -p "+dir+"/"+"topol.top -o "+dir+"/"+"npt.tpr"
        if(version<5.0):
            pre_phase2_command="grompp -f ./scripts/gromacs_less_than_5/npt.mdp"+pre_phase2_command
        elif(version>=5.0):
            pre_phase2_command="gmx grompp -f ./scripts/gromacs_5_or_more/npt.mdp"+pre_phase2_command

        status1=os.system(pre_phase2_command)
        if(status1 == 0):
            status2=os.system("mdrun -v -deffnm "+dir+"/"+"npt")
            if(status2 == 0):
                os.system("gmx energy -f "+dir+"/"+"npt.edr -o "+dir+"/"+"pressure.xvg")
                self.drawPressurevsTime(event)
                self.SetStatusText('Completed NVT equilibration succesfully!\n')
        else:
            print "ERROR: Hey dude! I think something went wrong in NPT simulation!"

    def OnProduction(self, event):
        pre_production_command="grompp -f md.mdp -c "+dir+"/"+"npt.gro -t "+dir+"/"+"npt.cpt -p "+dir+"/"+"topol.top -o "+dir+"/"+"md_0_1.tpr"
        if(version<5.0):
            pre_production_command=pre_production_command
        elif(version>=5.0):
            pre_production_command="gmx "+pre_production_command

        status1=os.system(pre_production_command)
        if(status1 != 0):
            print "ERROR: Hey dude! I think something went wrong in Pre Production run!"
    
    def OnQuit(self, e):
        self.Close()

class p1(wx.Panel):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="title", size=(50, 50))

        #configure graph
        self.figure = matplotlib.figure.Figure()
        self.axes   = self.figure.add_subplot(111)
        t=numpy.arange(0.0, 10, 1.0)
        s=[0,1,0,1,0,2,1,2,1,0]
        self.y_max =10
        self.axes.plot(t,s)
        self.canvas=FigureCanvas(self, -1, self.figure)

class MyApp(wx.App):
    def OnInit(self):
        frame = MyMainWindow(None, -1, 'Quick and Dirty Gromacs Tutorial')
#        frame = TestFrame(None, 'Hello World!')
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
