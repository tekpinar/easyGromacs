#!/usr/bin/env python
import wx
import os
import webbrowser
import numpy
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as NavigationToolbar

solvent_model=""
base=""
dir=""
version=5.02
#version=4.5

class IonizationDialog(wx.Dialog):
    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Ionization Setup")
        self.SetSize((350, 225))
        pos_lst = ["Na", "K", "Ca", "Mg"]
        self.rbox1=wx.RadioBox(self, wx.ID_ANY, "Select positive ion", (20, 60), wx.DefaultSize, pos_lst, 2, wx.RA_SPECIFY_COLS)
        neg_lst = ["Cl"]
        self.rbox2=wx.RadioBox(self, wx.ID_ANY, "Select negative ion", (200, 60), wx.DefaultSize, neg_lst, 2, wx.RA_SPECIFY_COLS)

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(self.rbox1)
        hbox1.Add(self.rbox2)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)

        self.chk = wx.CheckBox(self, label="Set concentration to ")
        self.chk.SetValue(True)

        hbox2.Add(self.chk, 0, wx.ALL, 5)

        self.concentration = wx.TextCtrl(self)
        self.concentration.SetValue("0.15")
        unit_text=wx.StaticText(self, label='M/L')
        hbox2.Add(self.concentration, 0, wx.ALL, 5)
        hbox2.Add(unit_text, 0, wx.ALL, 5)

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(hbox1, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(hbox2, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
#        vbox.Add(self.btns, 0, wx.ALL | wx.EXPAND, 5)
        vbox.Add(self.btns, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)

class SolvationDialog(wx.Dialog):
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Solvation Setup")
        self.SetSize((405, 200))
        box_lst = ["Triclinic", "Dodecahedron", "Cubic"]
        self.rbox=wx.RadioBox(self, wx.ID_ANY, "Select box type", (0, 60), wx.DefaultSize, box_lst, 1, wx.RA_SPECIFY_ROWS)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        buffer_text=wx.StaticText(self, label='Set buffer distance: ')

        self.buffer_distance = wx.TextCtrl(self)
        self.buffer_distance.SetValue("1.0")
        unit_text=wx.StaticText(self, label='nm')
        hbox.Add(buffer_text, 0, wx.ALL, 5)
        hbox.Add(self.buffer_distance, 0, wx.ALL, 5)
        hbox.Add(unit_text, 0, wx.ALL, 5)

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.rbox, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(hbox, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(self.btns, 0, wx.ALL | wx.EXPAND, 5)

        self.SetSizer(vbox)

class PreparationDialog(wx.Dialog):
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Preparation Setup")
        self.SetSize((270, 250))
        water_lst = ["TIP3P", "SPC", "TIP4P", "SPCE"]
        self.rbox=wx.RadioBox(self, wx.ID_ANY, "Select water model", (0, 60), wx.DefaultSize, water_lst, 2, wx.RA_SPECIFY_ROWS)
        
        ffm = ['amber03', 'amber94', 'amber96', 'amber99', 'amber99sb', \
                   'amber99sb-ildn', 'amberGS', 'charmm27', 'gromos43a1', 'gromos43a2', \
                   'gromos45a3', 'gromos53a5','gromos53a6', 'gromos54a7', 'oplsaa']
    
        self.cb = wx.ComboBox(self, pos=(50, 30), choices=ffm, style=wx.CB_READONLY)
        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

        text=wx.StaticText(self, label='Select force field: ')

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.rbox, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(text,      flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(self.cb,   flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
        vbox.Add(self.btns, flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)


class MyMainWindow(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(800, 600))

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
        toolbar.AddSimpleTool(2, wx.Image('images/stock_open.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Open', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(7, wx.Image('images/prepare.png', wx.BITMAP_TYPE_PNG).Rescale(70, 40).ConvertToBitmap(), 'Prepare Model', '')
        toolbar.AddSimpleTool(9, wx.Image('images/solvate.png', wx.BITMAP_TYPE_PNG).Rescale(70, 40).ConvertToBitmap(), 'Solvate', '')
        toolbar.AddSimpleTool(10, wx.Image('images/ionize.png', wx.BITMAP_TYPE_PNG).Rescale(65, 40).ConvertToBitmap(), 'Ionize', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(11, wx.Image('images/minimize.png', wx.BITMAP_TYPE_PNG).Rescale(85, 40).ConvertToBitmap(), 'Minimize', '')
        toolbar.AddSimpleTool(12, wx.Image('images/equilibrate_phase1.png', wx.BITMAP_TYPE_PNG).Rescale(110, 40).ConvertToBitmap(), 'Equilibrate-Phase 1', '')
        toolbar.AddSimpleTool(13, wx.Image('images/equilibrate_phase2.png', wx.BITMAP_TYPE_PNG).Rescale(110, 40).ConvertToBitmap(), 'Equilibrate-Phase 2', '')
        toolbar.AddSimpleTool(14, wx.Image('images/preduction.png', wx.BITMAP_TYPE_PNG).Rescale(90, 40).ConvertToBitmap(), 'Production Run', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(4, wx.Image('images/stock_exit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Exit', '')
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
        self.Bind(wx.EVT_TOOL, self.openfile, id=2)
#        self.Bind(wx.EVT_TOOL, self.OnSave, id=3)
        self.Bind(wx.EVT_TOOL, self.OnExit, id=4)
        self.Bind(wx.EVT_TOOL, self.OnPrepare, id=7)
        self.Bind(wx.EVT_TOOL, self.OnSolvate, id=9)
        self.Bind(wx.EVT_TOOL, self.OnIonize, id=10)
        self.Bind(wx.EVT_TOOL, self.OnMinimize, id=11)
        self.Bind(wx.EVT_TOOL, self.OnEquilibratePhase1, id=12)
        self.Bind(wx.EVT_TOOL, self.OnEquilibratePhase2, id=13)
        self.Bind(wx.EVT_TOOL, self.OnProduction, id=14)

        #Menubar
        menubar = wx.MenuBar()

        fileMenu = wx.Menu()
        open_menu_item=fileMenu.Append(wx.ID_OPEN, '&Open')
        fileMenu.AppendSeparator()
        quit_menu_item=fileMenu.Append(103, '&Quit', 'Quit Application')

        helpMenu=wx.Menu()
        tutorial_menu_item=helpMenu.Append(wx.ID_ANY, '&Justin Lemkul\'s Tutorial Web Page')
        about_menu_item=helpMenu.Append(wx.ID_ANY, '&About')

        menubar.Append(fileMenu, '&File')
        menubar.Append(helpMenu, '&Help')
        
        #Menubar events
        self.Bind(wx.EVT_MENU, self.openfile, open_menu_item)
        self.Bind(wx.EVT_MENU, self.OnExit, quit_menu_item, id=103)
        self.Bind(wx.EVT_MENU, self.OnAboutDlg, about_menu_item)
        self.Bind(wx.EVT_MENU, self.OnOpenTutorial, tutorial_menu_item)
        self.SetMenuBar(menubar)

    def OnAboutDlg(self, event):
        info = wx.AboutDialogInfo()
        info.Name = "Quick and Dirty Gromacs"
        info.Version = "0.0.3 Beta"
        info.Copyright = "(C) 2014 Mustafa Tekpinar\nEmail: tekpinar@buffalo.edu\nLicence: LGPL"
        wx.AboutBox(info)

    def OnOpenTutorial(self, event):
        url="http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/lysozyme/"
        webbrowser.open_new(url)

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

    def OnExit(self, event):
        self.Close(True)

    def OnNew(self, event):
        self.statusbar.SetStatusText('New Command')

#    def OnSave(self, event):
#        self.statusbar.SetStatusText('Save Command')
    
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

#    def waterchoice(self, event):
#        wm = ['SPC', 'SPCE', 'TIP3P', 'TIP4P']
#        dlg = wx.SingleChoiceDialog(self, 'Select water model', 'Which one?', wm, wx.CHOICEDLG_STYLE)
#        if dlg.ShowModal() == wx.ID_OK:
#            self.SetStatusText('You chose: %s\n' % dlg.GetStringSelection())
#        dlg.Destroy()
#        return dlg.GetStringSelection().lower()

    def ffchoice(self, event):
        ffm = ['amber03            -AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)',\
                   'amber94       -AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)',\
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
        global solvent_model 
        force_field=""
        initial_pdb=dir+"/"+base
###############################
        dlg = PreparationDialog()
        if (dlg.ShowModal()==wx.ID_OK):
            solvent_model=dlg.rbox.GetStringSelection().lower()
            force_field=dlg.cb.GetValue()
            self.SetStatusText('You selected %s force field and %s water model. \n' % (solvent_model, force_field))
            dlg.Destroy()
################################

         #This / makes it very unix dependent!!!
            if(dir == ""):
                dial=wx.MessageDialog(None, 'You have not opened a pdb file to start!\nGo to File-Open menu!', 'ERROR', wx.OK | wx.ICON_ERROR)
                dial.ShowModal()
                dial.Destroy()
            else:
                if(force_field !=""):
                #This / makes it very unix dependent!!!
                    prepare_command="pdb2gmx -f "+initial_pdb+" -o "+dir+"/"+"processed.pdb -i "+dir+"/"+"posre.itp -p "+dir+"/"+"topol.top -water "\
                        +solvent_model.lower()+" -ff "+force_field
                    if(version<5.0):
                        os.system(prepare_command)
                    elif(version>=5):
                        os.system("gmx "+prepare_command)
    def OnSolvate(self, event):
#####################
        boxtype=""
        buffer=""

        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           

        if os.path.isfile(dir+"/solvated.pdb"):
            warning = wx.MessageDialog(None, "It seems like you've already solvated the system!\n Are you sure about solvating the system again?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_NO):
                return (-1)
            warning.Destroy()           

        dlg = SolvationDialog()
        if (dlg.ShowModal()==wx.ID_OK):
            boxtype=dlg.rbox.GetStringSelection().lower()
            buffer=dlg.buffer_distance.GetValue()
            self.SetStatusText('You selected %s box and %s buffer distance. \n' % (boxtype, buffer))
        dlg.Destroy()
#####################

        #Define box properties
        #boxtype=self.boxtypechoice(event)
        #buffer=self.set_buffer_distance(event)
        error_message="ERROR: Hey dude! Something went wrong in solvation precedure! Check log files!"
        presolvate_command="editconf -f "+dir+"/"+"processed.pdb -o "+dir+"/"+"newbox.pdb -bt "+boxtype+" -c -d "+buffer
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

    def OnIonize(self, event):
        posIon=""
        negIon=""
        ion_concentration=""

        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           


        if os.path.isfile(dir+"/solv_ions.pdb"):
            warning = wx.MessageDialog(None, "It seems like you've already ionized the system!\nYou can continue from 'Minimize' or\n reset the system by  redoing 'Solvate'!",'Warning', wx.OK | wx.ICON_INFORMATION)
            warning.ShowModal()
            return (-1)
            warning.Destroy()           

        dlg = IonizationDialog()
        if (dlg.ShowModal()==wx.ID_OK):
            posIon=dlg.rbox1.GetStringSelection().upper()
            negIon=dlg.rbox2.GetStringSelection().upper()
            if(dlg.chk.GetValue()):
                ion_concentration=dlg.concentration.GetValue()
            self.SetStatusText('You selected %s and %s ions. \n' % (posIon, negIon))
        
            #Prepare for ion addition
            preprocess_command="-c "+dir+"/"+"solvated.pdb -p "+dir+"/"+"topol.top -o "+dir+"/"+"ions.tpr"
            ionize_command="genion -s "+dir+"/"+"ions.tpr -o "+dir+"/"+"solv_ions.pdb -p "+dir+"/"+"topol.top -pname "+posIon+" -nname "+negIon+" -neutral"
            
            if(ion_concentration !=""):
                ionize_command=ionize_command+" -conc 0.15"    

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
        dlg.Destroy()

    def set_temperature(self, event):
        dlg = wx.TextEntryDialog(self, 'Enter temperature in Kelvin','')
        dlg.SetValue("300")
        filename=""
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You entered: %s\n' % dlg.GetValue())
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/nvt.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/nvt.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            file.close()
            file=open(filename, "w")
            for line in lines:
                if( (line.find("ref_t") ) !=(-1)):
                    line="ref_t\t\t= "+dlg.GetValue()+"\t"+dlg.GetValue()+" \t; reference temperature, one for each group, in K\n"
                if( (line.find("gen_temp") ) !=(-1)):
                    line="gen_temp\t= "+dlg.GetValue()+"\t\t; temperature for Maxwell distribution\n"
                file.write(line)            
            file.close()
            if(version<5.0):
                file_name="./scripts/gromacs_less_than_5/nvt.mdp"
            elif(version>=5.0):
                file_name="./scripts/gromacs_5_or_more/nvt.mdp"
#            return 0
#        else:
#            return (-1)

#Now, lets set npt.mdp script also to the same temperature. 
######################################################################################################################################
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/npt.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/npt.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            file.close()
            file=open(filename, "w")
            for line in lines:
                if( (line.find("ref_t") ) !=(-1)):
                    line="ref_t\t\t= "+dlg.GetValue()+"\t"+dlg.GetValue()+" \t; reference temperature, one for each group, in K\n"
                file.write(line)            
            file.close()
######################################################################################################################################

#Finally, lets set md.mdp script also to the same temperature. 
######################################################################################################################################
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/md.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/md.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            file.close()
            file=open(filename, "w")
            for line in lines:
                if( (line.find("ref_t") ) !=(-1)):
                    line="ref_t\t\t= "+dlg.GetValue()+"\t"+dlg.GetValue()+" \t; reference temperature, one for each group, in K\n"
                file.write(line)            
            file.close()
######################################################################################################################################
            return dlg.GetValue()
        else:
            return (-1)
        dlg.Destroy()

    def set_pressure(self, event):
        dlg = wx.TextEntryDialog(self, 'Unit: Bar', 'Enter Pressure')
        dlg.SetValue("1.0")
        filename=""
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You entered: %s\n' % dlg.GetValue())
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/npt.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/npt.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            file.close()
            file=open(filename, "w")
            for line in lines:
                if( (line.find("ref_p") ) !=(-1)):
                    line="ref_p\t\t\t= "+dlg.GetValue()+"\t\t\t; reference pressure, in bar\n"
                file.write(line)            
            file.close()
#Finally, lets set md.mdp script also to the same pressure.
######################################################################################################################################
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/md.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/md.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            file.close()
            file=open(filename, "w")

            for line in lines:
                if( (line.find("ref_p") ) !=(-1)):
                    line="ref_p\t\t\t= "+dlg.GetValue()+"\t\t\t; reference pressure, in bar\n"

                file.write(line)            
            file.close()
######################################################################################################################################

            return dlg.GetValue()
        else:
            return (-1)
        dlg.Destroy()

    def set_simulation_time(self, event):
        dlg = wx.TextEntryDialog(self, 'Unit: nanosecond','Simulation Time')
        dlg.SetValue("1.0")
        filename=""
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText('You entered: %s\n' % dlg.GetValue())
            #Finally, lets set total simulation time in md.mdp.
            if(version<5.0):
                filename="./scripts/gromacs_less_than_5/md.mdp"
            elif(version>=5.0):
                filename="./scripts/gromacs_5_or_more/md.mdp"

            file=open(filename, "r")
            lines=file.readlines()
            #Get time step now
            time_step=""
            for line in lines:
                if( (line.find("dt") ) !=(-1)):
                    time_step=line.partition('=')[-1].rpartition(';')[0]
                    print "This is time step"+time_step
            file.close()
        
            number_of_steps=1000*int(float(dlg.GetValue())/float(time_step))
            file=open(filename, "w")

            for line in lines:
                if( (line.find("nsteps") ) !=(-1)):
                    line="nsteps          = "+str(number_of_steps)+"\t\t; 2 * "+str(number_of_steps)+" =  ("+dlg.GetValue()+" ns)\n"
                file.write(line)            
            file.close()
            return dlg.GetValue()
        else:
            return (-1)
        dlg.Destroy()

    def OnMinimize(self, event):
        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           

        if (os.path.isfile(dir+"/em.tpr")) and (os.path.isfile(dir+"/em.log") ) and (os.path.isfile(dir+"/em.gro")):
            warning = wx.MessageDialog(None, "It seems like you've already minimized the system!\n Are you sure about minimizing the system again?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_NO):
                return (-1)
            warning.Destroy()           

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
        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           

        if (os.path.isfile(dir+"/nvt.tpr")) and (os.path.isfile(dir+"/nvt.log") ) and (os.path.isfile(dir+"/nvt.gro")):
            warning = wx.MessageDialog(None, "It seems like you've already performed NVT equilibration!\n Are you sure about doing it again?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_NO):
                return (-1)
            warning.Destroy()           

        status=self.set_temperature(event)
        if(status!=(-1)):
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

    def OnEquilibratePhase2(self, event):
        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           

        if (os.path.isfile(dir+"/npt.tpr")) and (os.path.isfile(dir+"/npt.log") ) and (os.path.isfile(dir+"/npt.gro")):
            warning = wx.MessageDialog(None, "It seems like you've already performed NPT equilibration!\n Are you sure about doing it again?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_NO):
                return (-1)
            warning.Destroy()           

        status=self.set_pressure(event)
        if(status!=(-1)):
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
                    print "ERROR: Hey dude! I think something went wrong in NPT results plotting precedure!"
            else:
                print "ERROR: Hey dude! I think something went wrong in NPT simulation!"
        
    def production_warning(self, event):
        dlg = wx.MessageDialog(None, "Congratulations! \nYou produced a md.tpr file. \nWe suggest you to perform\n production run on a cluster!",'Warning',wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
    
    def OnProduction(self, event):
        global dir
        if ( dir == ""):
            warning = wx.MessageDialog(None, "Hey dude! It seems like you've not opened your project folder. \nWould you like to do that?",\
                                           'Warning', wx.YES_NO | wx.ICON_INFORMATION)
            if(warning.ShowModal()==wx.ID_YES):
                dir_dlg = wx.DirDialog(self, "Choose your working directory:")
                if dir_dlg.ShowModal() == wx.ID_OK:
                    self.SetStatusText('You\'ve chosen : %s\n' % dir_dlg.GetPath())
                    dir = dir_dlg.GetPath()
                    dir_dlg.Destroy()
                else:
                    return (-1)
            else:
                return (-1)
            warning.Destroy()           
        status=self.set_simulation_time(event)
        if(status!=(-1)):
            pre_production_command = dir+"/"+"npt.gro -t "+dir+"/"+"npt.cpt -p "+dir+"/"+"topol.top -o "+dir+"/"+"md_0_1.tpr"
            if(version<5.0):
                pre_production_command="grompp -f ./scripts/gromacs_less_than_5/md.mdp -c "+pre_production_command
            elif(version>=5.0):
                pre_production_command="gmx grompp -f ./scripts/gromacs_5_or_more/md.mdp -c "+pre_production_command

            status1=os.system(pre_production_command)
            if(status1==0):
                self.production_warning(event)
            else:
                print "ERROR: Hey dude! I think something went wrong in Pre Production run!"
    
    def OnQuit(self, e):
        self.Close()

class MyApp(wx.App):
    def OnInit(self):
        frame = MyMainWindow(None, -1, 'Quick and Dirty Gromacs')
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
