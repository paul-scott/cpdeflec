import wx

class ScrollWinEvent(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title)
        panel = wx.Panel(self, -1)
        #self.st = wx.StaticText(panel, -1, '0', (30,0))
        panel.Bind(wx.EVT_SCROLLWIN, self.OnScroll)
        panel.SetScrollbar(wx.VERTICAL, 1, 6, 50)
        self.Centre()
        self.Show(True)

    def OnScroll(self, evt):
        y = evt.GetPosition()
        #self.st.SetLabel(str(y)) 
app = wx.App()
ScrollWinEvent(None, -1, 'scrollwinevent.py')
app.MainLoop()

