#!/usr/bin/env python
import wx
import wx.lib.filebrowsebutton
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

class Energy_Options_Check_Dialog(wx.Dialog):
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Check Available Energy Options")
        self.SetSize((550, 100))

        self.fbb = wx.lib.filebrowsebutton.FileBrowseButton(self, size=(550, -1), labelText="Select an energy file:", fileMask="*.edr")

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.fbb, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)        
        vbox.Add(self.btns, flag=wx.ALIGN_RIGHT|wx.TOP|wx.BOTTOM, border=5)
        self.SetSizer(vbox)

    def openedrfile(self, event):
       dlg = wx.FileDialog(self, "Choose an energy file", os.getcwd(), "", "*.edr", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           return dlg.GetPath()
       dlg.Destroy()

class Energy_Plot_Dialog(wx.Dialog):
    def __init__(self, energy_list):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Plot Energy")
        self.SetSize((550, 150))

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        outfile_text=wx.StaticText(self, label='Output file name   :')
        self.outfile = wx.TextCtrl(self)
        self.outfile.SetValue("energy.xvg")
        hbox1.Add(outfile_text, 0, wx.ALL, 5)
        hbox1.Add(self.outfile, 0, wx.ALL, 5)        

        self.cb = wx.ComboBox(self, choices=energy_list, style=wx.CB_READONLY)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        outfile_text2=wx.StaticText(self, label='Select an option    :')
        hbox2.Add(outfile_text2, 0, wx.ALL, 5)
        hbox2.Add(self.cb, 0, wx.ALL, 5)

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(hbox1, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)
        vbox.Add(hbox2, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)
        vbox.Add(self.btns, flag=wx.ALIGN_RIGHT|wx.TOP|wx.BOTTOM, border=5)

        self.SetSizer(vbox)

class RMS_D_and_F_Dialog(wx.Dialog):
    def __init__(self, which_box_string):
        """Constructor"""
        wx.Dialog.__init__(self, None, title=which_box_string+" Setup")
        self.SetSize((550, 290))
        type_lst = ["CA", "Backbone", "All atoms"]
        self.rbox1=wx.RadioBox(self, wx.ID_ANY, "", (20, 10), wx.DefaultSize, type_lst, 3, wx.RA_SPECIFY_COLS)
        sel_atms_text=wx.StaticText(self, label=' Selected atoms        :')

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(sel_atms_text)
        hbox1.Add(self.rbox1)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        init_frm_text=wx.StaticText(self, label=' Initial frame             : ')
        self.num_entry1 = wx.TextCtrl(self)
        self.num_entry1.SetValue("0")

        fnl_frm_text=wx.StaticText(self, label='Final frame:')
        self.num_entry2 = wx.TextCtrl(self)
        self.num_entry2.SetValue("-1")

        hbox2.Add(init_frm_text, 0, wx.ALL, 0)
        hbox2.Add(self.num_entry1, 0, wx.ALL, 0)

        hbox2.Add(fnl_frm_text, 0, wx.ALL, 0)
        hbox2.Add(self.num_entry2, 0, wx.ALL, 0)

        self.fbb1 = wx.lib.filebrowsebutton.FileBrowseButton(self, size=(550, -1), labelText="Select a reference file:", fileMask="*.pdb")
        self.fbb2 = wx.lib.filebrowsebutton.FileBrowseButton(self, size=(550, -1), labelText="Select a trajectory file:", fileMask="*.xtc")

        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        outfile_text=wx.StaticText(self, label='Output file name     :')
        self.outfile = wx.TextCtrl(self)
        self.outfile.SetValue(which_box_string.lower()+".xvg")
        hbox5.Add(outfile_text, 0, wx.ALL, 5)
        hbox5.Add(self.outfile, 0, wx.ALL, 5)        

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
        vbox = wx.BoxSizer(wx.VERTICAL)

        vbox.Add(self.fbb1, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)        
        vbox.Add(self.fbb2, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)        
        vbox.Add(hbox5, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)

        vbox.Add(hbox1, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)
        vbox.Add(hbox2, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)

        vbox.Add(self.btns, flag=wx.ALIGN_RIGHT|wx.TOP|wx.BOTTOM, border=5)

        self.SetSizer(vbox)

    def openrefpdb(self, event):
       dlg = wx.FileDialog(self, "Choose an initial pdb file", os.getcwd(), "", "*.pdb", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           return dlg.GetPath()
       dlg.Destroy()

    def opentrajectory(self, event):
       dlg = wx.FileDialog(self, "Choose an initial xtc file", os.getcwd(), "", "*.xtc", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           return dlg.GetPath()
       dlg.Destroy()

class PCA_Dialog(wx.Dialog):
    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="PCA Setup")
        self.SetSize((550, 380))
        type_lst = ["CA", "Backbone", "All atoms"]
        self.rbox0=wx.RadioBox(self, wx.ID_ANY, "", (20, 10), wx.DefaultSize, type_lst, 3, wx.RA_SPECIFY_COLS)
        sel_atms_text=wx.StaticText(self, label=' Selected atoms        :')

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(sel_atms_text)
        hbox1.Add(self.rbox0)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        init_frm_text=wx.StaticText(self, label=' Initial frame             : ')
        self.num_entry1 = wx.TextCtrl(self)
        self.num_entry1.SetValue("0")

        fnl_frm_text=wx.StaticText(self, label='Final frame:')
        self.num_entry2 = wx.TextCtrl(self)
        self.num_entry2.SetValue("-1")

        hbox2.Add(init_frm_text, 0, wx.ALL, 0)
        hbox2.Add(self.num_entry1, 0, wx.ALL, 0)

        hbox2.Add(fnl_frm_text, 0, wx.ALL, 0)
        hbox2.Add(self.num_entry2, 0, wx.ALL, 0)

        self.fbb1 = wx.lib.filebrowsebutton.FileBrowseButton(self, size=(550, -1), labelText="Select a reference file:", fileMask="*.pdb")
        self.fbb2 = wx.lib.filebrowsebutton.FileBrowseButton(self, size=(550, -1), labelText="Select a trajectory file:", fileMask="*.xtc")

############################################
        pc_lst = ["1", "2", "3", "4"]
        self.rbox1=wx.RadioBox(self, wx.ID_ANY, "", (20, 10), wx.DefaultSize, pc_lst, 4, wx.RA_SPECIFY_COLS)
        sel_pc1_text=wx.StaticText(self, label=' The first PC              :')

        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        hbox3.Add(sel_pc1_text)
        hbox3.Add(self.rbox1)

        self.rbox2=wx.RadioBox(self, wx.ID_ANY, "", (20, 10), wx.DefaultSize, pc_lst, 4, wx.RA_SPECIFY_COLS)
        self.rbox2.SetSelection(1)  
        sel_pc2_text=wx.StaticText(self, label=' The second PC         :')

        hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        hbox4.Add(sel_pc2_text)
        hbox4.Add(self.rbox2)

############################################
        hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        outfile_text=wx.StaticText(self, label='Output file name     :')
        self.outfile = wx.TextCtrl(self)
        self.outfile.SetValue("pca.xvg")
        hbox5.Add(outfile_text, 0, wx.ALL, 5)
        hbox5.Add(self.outfile, 0, wx.ALL, 5)        

        self.btns = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
        vbox = wx.BoxSizer(wx.VERTICAL)

        vbox.Add(self.fbb1, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)        
        vbox.Add(self.fbb2, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)        
        vbox.Add(hbox5, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)

        vbox.Add(hbox1, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)
        vbox.Add(hbox2, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)

        vbox.Add(hbox3, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)
        vbox.Add(hbox4, flag=wx.ALIGN_LEFT|wx.TOP|wx.BOTTOM, border=5)

        vbox.Add(self.btns, flag=wx.ALIGN_RIGHT|wx.TOP|wx.BOTTOM, border=5)

        self.SetSizer(vbox)

    def openrefpdb(self, event):
       dlg = wx.FileDialog(self, "Choose an initial pdb file", os.getcwd(), "", "*.pdb", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           return dlg.GetPath()
       dlg.Destroy()

    def opentrajectory(self, event):
       dlg = wx.FileDialog(self, "Choose an initial xtc file", os.getcwd(), "", "*.xtc", wx.OPEN)
       if dlg.ShowModal() == wx.ID_OK:
           return dlg.GetPath()
       dlg.Destroy()

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
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(900, 675))

        #Box parameters
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(vbox)

        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.plottoolbar = NavigationToolbar(self.canvas)

        self.button1 = wx.Button(self, -1, "Plot RMSD")
        self.button1.Bind(wx.EVT_BUTTON, self.OnPlotRMSD)

        self.button2 = wx.Button(self, -1, "Plot RMSF")
        self.button2.Bind(wx.EVT_BUTTON, self.OnPlotRMSF)

        self.button3 = wx.Button(self, -1, "Plot PCA")
        self.button3.Bind(wx.EVT_BUTTON, self.OnPlotPCA)

        self.button4 = wx.Button(self, -1, "Plot Energy")
        self.button4.Bind(wx.EVT_BUTTON, self.OnPlotEnergy)

        self.button5 = wx.Button(self, -1, "Plot Rg")
#        self.button5.Bind(wx.EVT_BUTTON, self.OnPlotEnergy)

        self.button6 = wx.Button(self, -1, "Plot Clusters")
#        self.button6.Bind(wx.EVT_BUTTON, self.OnPlotEnergy)

        self.button7 = wx.Button(self, -1, "Plot SS")
#        self.button7.Bind(wx.EVT_BUTTON, self.OnPlotEnergy)


        #Toolbar items
        toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
#        toolbar.AddSimpleTool(2, wx.Image('images/stock_open.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Open', '')
        toolbar.AddSimpleTool(2, wx.Image('icons/open/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Open', '')
        toolbar.AddSeparator()
#        toolbar.AddSimpleTool(7, wx.Image('new_images/prepare.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Prepare Model', '')
#        toolbar.AddSimpleTool(7, wx.Image('new_images/prepare.png', wx.BITMAP_TYPE_ANY).Rescale(85, 50).ConvertToBitmap(), 'Prepare Model', '')
        toolbar.AddSimpleTool(7, wx.Image('icons/prepare/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Prepare Model', '')
#        toolbar.AddSimpleTool(7, wx.Image('images/Prepare_New2.png', wx.BITMAP_TYPE_PNG).Rescale(56, 50).ConvertToBitmap(), 'Prepare Model', '')
        toolbar.AddSimpleTool(9, wx.Image('icons/solvate/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Solvate', '')
#        toolbar.AddSimpleTool(9, wx.Image('png_2/solvate.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Solvate', '')
#        toolbar.AddSimpleTool(9, wx.Image('images/Solvate.png', wx.BITMAP_TYPE_PNG).Rescale(64, 80).ConvertToBitmap(), 'Solvate', '')
        toolbar.AddSimpleTool(10, wx.Image('icons/ionize/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Ionize', '')
#        toolbar.AddSimpleTool(10, wx.Image('png_2/ionize.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Ionize', '')
#        toolbar.AddSimpleTool(10, wx.Image('images/Ionize_new.png', wx.BITMAP_TYPE_PNG).Rescale(56, 80).ConvertToBitmap(), 'Ionize', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(11, wx.Image('icons/minimize/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Minimize', '')
#        toolbar.AddSimpleTool(11, wx.Image('png_2/minimize.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Minimize', '')
        toolbar.AddSimpleTool(12, wx.Image('icons/equilibrate_nvt/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Equilibrate-Phase 1', '')
#        toolbar.AddSimpleTool(12, wx.Image('png_2/equilibrate_nvt.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Equilibrate-Phase 1', '')
        toolbar.AddSimpleTool(13, wx.Image('icons/equilibrate_npt/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Equilibrate-Phase 2', '')
 
#        toolbar.AddSimpleTool(12, wx.Image('png_2/equilibrate_npt.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Equilibrate-Phase 1', '')
        toolbar.AddSimpleTool(14, wx.Image('icons/preduction/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Preduction Run', '')
#        toolbar.AddSimpleTool(14, wx.Image('png_2/preduction.png', wx.BITMAP_TYPE_PNG).Rescale(85, 50).ConvertToBitmap(), 'Preduction Run', '')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(4, wx.Image('icons/close/96x96.png', wx.BITMAP_TYPE_PNG).Rescale(96, 96).ConvertToBitmap(), 'Exit', '')
#        toolbar.AddSimpleTool(4, wx.Image('png_2/close.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Exit', '')
        toolbar.Realize()

        self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer2.Add(self.button1, -1, wx.EXPAND)
        self.sizer2.Add(self.button2, -1, wx.EXPAND)
        self.sizer2.Add(self.button3, -1, wx.EXPAND)
        self.sizer2.Add(self.button4, -1, wx.EXPAND)
        self.sizer2.Add(self.button5, -1, wx.EXPAND)
        self.sizer2.Add(self.button6, -1, wx.EXPAND)
        self.sizer2.Add(self.button7, -1, wx.EXPAND)


        vbox.Add(toolbar, 0, wx.EXPAND)
        vbox.Add(self.sizer2, 0, wx.EXPAND)
        vbox.Add(self.plottoolbar, 0, wx.EXPAND)
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
        self.Bind(wx.EVT_ENTER_WINDOW, self.OnEnter, id=7)
        self.Bind(wx.EVT_LEAVE_WINDOW, self.OnLeave, id=7)



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

    def OnEnter(self, event):
        btn = event.GetEventObject()   
#self.panel.SetBackgroundColour('Green')
#        self.panel.Refresh()     
        btn.SetBackgroundColour('Green')
        btn.Refresh()
        
    def OnLeave(self, event):
        #            btn = e.GetEventObject()
        self.SetBackgroundColour('DARKGREY', id=7)
#            btn.Refresh()

    def OnAboutDlg(self, event):
        info = wx.AboutDialogInfo()
        info.Name = "Easy GROMACS"
        info.Version = "0.0.5 Beta"
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
    
    def OnPlotRMSD(self, event):
        dlg = RMS_D_and_F_Dialog("RMSD")
        if (dlg.ShowModal()==wx.ID_OK):
            atom_selection=dlg.rbox1.GetStringSelection().upper()
            beg_frame=dlg.num_entry1.GetValue()
            end_frame=dlg.num_entry2.GetValue()
            if(atom_selection == "CA"):
                echo_string="echo 3 3|"
            elif (atom_selection == "BACKBONE"):
                echo_string="echo 4 4|"
            else:
                echo_string=""

            ref_file=dlg.fbb1.GetValue()
            trj_file=dlg.fbb2.GetValue()
            out_file=dlg.outfile.GetValue()
        
            #Prepare for rmsd calculation
            if(end_frame!="-1"):
                rmsd_command="rms -s "+ref_file+" -f "+trj_file+" -o "+out_file+" -tu ns "+"-b "+beg_frame+" -e "+end_frame
            else:
                rmsd_command="rms -s "+ref_file+" -f "+trj_file+" -o "+out_file+" -tu ns "+"-b "+beg_frame

            print rmsd_command
            error_message="ERROR: Something went wrong in rmsd procedure! Check log files!"
            if(version<5.0):
                status=os.system(echo_string+"g_"+rmsd_command)
            elif(version>=5.0):
                status=os.system(echo_string+"gmx "+rmsd_command)

            if(status==0):
           #Plot RMSD on canvas
                self.drawBlahvsBlah(out_file, "Time (ns)", "RMSD (nm)")
            else:
                print (error_message)

        dlg.Destroy()

####################
    def OnPlotRMSF(self, event):
        dlg = RMS_D_and_F_Dialog("RMSF")
        if (dlg.ShowModal()==wx.ID_OK):
            atom_selection=dlg.rbox1.GetStringSelection().upper()
            beg_frame=dlg.num_entry1.GetValue()
            end_frame=dlg.num_entry2.GetValue()
            if(atom_selection == "CA"):
                echo_string="echo 3 3|"
            elif (atom_selection == "BACKBONE"):
                echo_string="echo 4 4|"
            else:
                echo_string=""

            ref_file=dlg.fbb1.GetValue()
            trj_file=dlg.fbb2.GetValue()
            out_file=dlg.outfile.GetValue()
        
            #Prepare for rmsf calculation
            if((end_frame!="-1") and (beg_frame!="0")):
                rmsf_command="rmsf -s "+ref_file+" -f "+trj_file+" -res -o "+out_file+"-b "+beg_frame+" -e "+end_frame
            else:
                rmsf_command="rmsf -s "+ref_file+" -f "+trj_file+" -res -o "+out_file

            print rmsf_command
            error_message="ERROR: Something went wrong in rmsf procedure! Check log files!"
            if(version<5.0):
                status=os.system(echo_string+"g_"+rmsf_command)
            elif(version>=5.0):
                status=os.system(echo_string+"gmx "+rmsf_command)

            if(status==0):
            #Plot RMSF on canvas
                self.drawBlahvsBlah(out_file, "Residues", "RMSF (nm)")
            else:
                print (error_message)

        dlg.Destroy()
####################
    def OnPlotPCA(self, event):
        dlg = PCA_Dialog()
        if (dlg.ShowModal()==wx.ID_OK):
            atom_selection=dlg.rbox0.GetStringSelection().upper()
            pcA=dlg.rbox1.GetStringSelection()
            pcB=dlg.rbox2.GetStringSelection()
            beg_frame=dlg.num_entry1.GetValue()
            end_frame=dlg.num_entry2.GetValue()
            if(atom_selection == "CA"):
                echo_string="echo 3 3|"
            elif (atom_selection == "BACKBONE"):
                echo_string="echo 4 4|"
            else:
                echo_string=""

            ref_file=dlg.fbb1.GetValue()
            trj_file=dlg.fbb2.GetValue()
            out_file=dlg.outfile.GetValue()
            pca_command1="covar -f "+trj_file+" -s "+ref_file+" -o eigenval.xvg -v eigenvec.trr"
            pca_command2="anaeig -f "+trj_file+" -s "+ref_file+" -v eigenvec.trr -2d "+out_file+" -first "+pcA+" -last "+pcB        

#            #Prepare for rmsf calculation
#            if((end_frame!="-1") and (beg_frame!="0")):
#                rmsf_command="rmsf -s "+ref_file+" -f "+trj_file+" -res -o "+out_file+"-b "+beg_frame+" -e "+end_frame
#            else:
#                rmsf_command="rmsf -s "+ref_file+" -f "+trj_file+" -res -o "+out_file

            error_message="ERROR: Something went wrong in PCA procedure! Check log files!"
            if(version<5.0):
                status1=os.system(echo_string+"g_"+pca_command1)
                status2=os.system(echo_string+"g_"+pca_command2)
            elif(version>=5.0):
                status1=os.system(echo_string+"gmx "+pca_command1)
                status2=os.system(echo_string+"gmx "+pca_command2)

            if((status1==0) and (status2==0)):
            #Plot RMSF on canvas
                self.drawPC_A_vs_PC_B(out_file, "PC"+pcA+" (nm)", "PC"+pcB+" (nm)")
            else:
                print (error_message)

        dlg.Destroy()

####################
    def OnPlotEnergy(self, event):
        edr_file=""
        out_file=""
        dlg = Energy_Options_Check_Dialog()
        if (dlg.ShowModal()==wx.ID_OK):

            edr_file=dlg.fbb.GetValue()

            #First, get available parameters!
            echo_string="echo 0 0|"
            check_options_command="energy -f "+edr_file+" 2>options.dat"

            if(version<5.0):
                status1=os.system(echo_string+"g_"+check_options_command)
            elif(version>=5.0):
                status1=os.system(echo_string+"gmx "+check_options_command)

        dlg.Destroy()
        
        options_file=open("options.dat", "r")
        lines=options_file.readlines()
            
        #Find lines between two dashed lines. 
        counter=0
        beg_and_end=[]
        
        for line in lines:
            if line[0]=='-':
                beg_and_end.append(counter)
            counter+=1
        energy_options=[]
        for i in range((beg_and_end[0]+1), beg_and_end[1]):
            energy_options.extend(lines[i].split())

        options_file.close()
        energy_options_pairs=[]
        for i in range (0, len(energy_options)/2):
            temp_string=str(energy_options[2*i])+" "+str(energy_options[2*i+1])
            energy_options_pairs.append(temp_string)

        dlg=Energy_Plot_Dialog(energy_options_pairs)
        if (dlg.ShowModal()==wx.ID_OK):

            out_file=dlg.outfile.GetValue()
            myindex=dlg.cb.GetValue().split()

            #Now, plot selected index
            echo_string="echo "+myindex[0]+" 0|"
            print echo_string
            energy_command="energy -f "+edr_file+" -o "+out_file

            if(version<5.0):
                status2=os.system(echo_string+"g_"+energy_command)
            elif(version>=5.0):
                status2=os.system(echo_string+"gmx "+energy_command)

            error_message="ERROR: Something went wrong in energy calculation procedure! Check log files!"

            if(status2==0):
            #Plot Energy on canvas
                self.drawBlahvsBlah(out_file, "Time", "Energy (kJ/mol)")
            else:
                print (error_message)

        dlg.Destroy()

    def drawBlahvsBlah(self, out_file, x_label_string, y_label_string):
        data_file=open(out_file, "r")
        lines=data_file.readlines()
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
                md_pc1.append(float(line.split()[0]))
                md_pc2.append(float(line.split()[1]))
                counter+=1

        self.axes.clear()
        self.axes.set_xlim([min(md_pc1),max(md_pc1)])
        self.axes.plot(md_pc1, md_pc2)
        self.axes.grid()
        self.axes.set_xlabel(x_label_string) 
        self.axes.set_ylabel(y_label_string) 
        self.axes.fill_between(md_pc1, 0, md_pc2, facecolor='blue', alpha=0.3)
        data_file.close()
        self.Layout()

    def drawPC_A_vs_PC_B(self, out_file, x_label_string, y_label_string):
        data_file=open(out_file, "r")
        lines=data_file.readlines()
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
                md_pc1.append(float(line.split()[0]))
                md_pc2.append(float(line.split()[1]))
                counter+=1

#        print "This file has %d lines"%counter
        md_time_percentage=[]
        for i in range(counter):
            time_percentage=0.0
            time_percentage=(float(i)/float(counter-1))
#            time_percentage=(float(i)/)
            md_time_percentage.append(time_percentage)

        self.axes.clear()
#        self.axes.set_xlim([min(md_pc1),max(md_pc1)])
#        self.axes.scatter(md_pc1, md_pc2)
        self.axes.scatter(md_pc1, md_pc2, marker='o', c=[md_time_percentage], s=10)    
#        cb=self.axes.colorbar()
#        cb.set_label('Time Percentage')

        self.axes.grid()

        self.axes.set_xlabel(x_label_string) 
        self.axes.set_ylabel(y_label_string) 
        data_file.close()
        self.Layout()

    def changePlot(self, event):
        if self.current_draw == 'sin' :
            self.drawLines()
            self.current_draw = 'lines'
        else: 
            self.drawPotential(out_file)
            self.current_draw = 'sin'

        self.Layout()

    def drawLines(self):
        x = numpy.arange(0.0,1000,100.0)
        y = [0,0,0,0,0,0,0,0,0,0]

        self.axes.clear()
        self.axes.plot(x, y)

    def drawPotential(self, out_file):
        data_file=open(out_file, "r")
        lines=data_file.readlines()
        md_pc1=[]
        md_pc2=[]
        counter=0
        for line in lines:
            if((line[0]!='#') and (line[0]!='@')):
                md_pc1.append(float(line.split()[0]))
                md_pc2.append(float(line.split()[1]))
                counter+=1

        self.axes.clear()
        self.axes.set_xlim([min(md_pc1),max(md_pc1)])
        self.axes.plot(md_pc1, md_pc2)
        self.axes.grid()
        self.axes.set_xlabel("Step") 
        self.axes.set_ylabel("Potential Energy (kJ/mol)") 
        data_file.close()
        self.Layout()

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
        error_message="ERROR: Hey dude! Something went wrong in solvation procedure! Check log files!"
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

                error_message="ERROR: Something went wrong in ionization procedure! Check log files!"
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
#                self.drawPotential(event)
                self.drawBlahvsBlah("potential.xvg", "Step", "Potential Energy (kJ/mol)")
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
                    print "ERROR: Hey dude! I think something went wrong in NPT results plotting procedure!"
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
        frame = MyMainWindow(None, -1, 'Easy GROMACS')
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
