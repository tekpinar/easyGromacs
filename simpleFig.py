import wx
import numpy 
import matplotlib

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as NavigationToolbar
class TestFrame(wx.Frame):
    def __init__(self,parent,title):
        wx.Frame.__init__(self,parent,title=title,size=(500,500))
        self.p2 = MatplotPanel(self)
#        self.statusbar = self.CreateStatusBar()
#        self.statusbar.SetStatusText('Oi')
      

class MatplotPanel(wx.Panel):

    def __init__(self, parent):     
        wx.Panel.__init__(self, parent,-1,size=(50,50))

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.vbox)

        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.plottoolbar = NavigationToolbar(self.canvas)

        self.toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        self.toolbar.AddSimpleTool(1, wx.Image('stock_exit.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap(), 'Exit', '')
        self.toolbar.Realize()


        self.button = wx.Button(self, -1, "Change plot")
        self.button.Bind(wx.EVT_BUTTON, self.changePlot)


        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.Add(self.plottoolbar, 0, wx.EXPAND)
        self.vbox.Add(self.button, 0, wx.EXPAND)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)

        self.drawSin()
        self.current_draw = 'sin'

#       self.drawLines()

    def changePlot(self, event):

        if self.current_draw == 'sin' :
            self.drawLines()
            self.current_draw = 'lines'
        else: 
            self.drawSin()
            self.current_draw = 'sin'

        self.Layout()

    def drawLines(self):

        x = numpy.arange(0.0,10,1.0)
        y = [0,1,0,1,0,2,1,2,1,0]

        self.axes.clear()
        self.axes.plot(x, y)

    def drawSin(self):

        x = numpy.arange(0.0,10,0.1)
        y = numpy.sin(x)

        self.axes.clear()
        self.axes.plot(x, y)

app = wx.App(redirect=False)
frame = TestFrame(None, 'Hello World!')
frame.Show()
app.MainLoop()
