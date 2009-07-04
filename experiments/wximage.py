import wx
import os

ID_OPEN=111
ID_EXIT=110

class MainFrame(wx.Frame):
	""" We simply derive a new class of Frame. """
	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title, size=(800,700))
		#self.Maximize()

		self.pan = ListPanel(self, -1)
		self.imwin = ImageWindow(self, -1)
		self.sizer = wx.BoxSizer(wx.HORIZONTAL)
		self.sizer.Add(self.pan, 0, 10)
		self.sizer.Add(self.imwin, 0, wx.EXPAND | wx.ALL)

		self.sb = self.CreateStatusBar()
		# Create menus.
		filemenu = wx.Menu()
		filemenu.Append(ID_OPEN, '&Open', 'Open file')
		filemenu.AppendSeparator()
		filemenu.Append(ID_EXIT, 'E&xit', 'Terminate the program')
		#exit = wx.MenuItem(filemenu, ID_EXIT, 'E&xit', 'Terminate the program')
		#exit.SetBitmap(wx.Bitmap('icons/exit.png'))
		#filemenu.AppendItem(exit)
		# Create menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		self.SetMenuBar(menuBar)
		wx.EVT_MENU(self, ID_OPEN, self.OnOpen)
		wx.EVT_MENU(self, ID_EXIT, self.OnExit)

		self.Bind(wx.EVT_PAINT, self.imwin.OnPaint)
		self.SetSizer(self.sizer)
		self.Show(True)
		
	def OnOpen(self, e):
		"""Open a file"""
		self.dirname = ''
		dlg = wx.FileDialog(self, "Choose a file", self.dirname,
				"", "*.*", wx.MULTIPLE)
		if dlg.ShowModal() == wx.ID_OK:
			self.filenames = dlg.GetFilenames()
			# Filenames as "panel_ID-h.tiff" or "panel_ID-v.tiff".
			self.dirname = dlg.GetDirectory()
			#f = open(os.path.join(self.dirname, self.filename), 'r')
			#f.close()
		dlg.Destroy()
	def OnExit(self, e):
		self.Close(True)

class ListPanel(wx.Panel):
	""" Window container for bitmap. """
	def __init__(self, parent, id):
		wx.Panel.__init__(self, parent, id)
		text = wx.StaticText(self, -1, 'Start of list')
	
class ImageWindow(wx.ScrolledWindow):
	""" Window container for bitmap. """
	def __init__(self, parent, id):
		wx.ScrolledWindow.__init__(self, parent, id)
		self.parent = parent
		self.list = []
		self.bitmap = wx.Bitmap('1d-tun-edit.tiff', 
				type=wx.BITMAP_TYPE_TIF)
		sz = self.bitmap.GetSize()
		incr = 40
		self.SetScrollbars(incr, incr, sz[0]/incr+1, sz[1]/incr+1)

		self.statbit = wx.StaticBitmap(self, -1, self.bitmap)
		#refresh calling this first?
		#might want to create transparent canvas on top of sw.

		#self.statbit.Bind(wx.EVT_MOTION, self.OnMove)
		#self.statbit.Bind(wx.EVT_LEFT_UP, self.OnClick)

		self.statbit.Bind(wx.EVT_MOTION, self.OnMove)
		self.statbit.Bind(wx.EVT_LEFT_UP, self.OnLClick)
		self.statbit.Bind(wx.EVT_RIGHT_UP, self.OnRClick)

	def OnPaint(self, e):
		dc = wx.PaintDC(self)
		dc.SetBrush(wx.Brush(wx.Colour(255,255,255), style=wx.TRANSPARENT))
		for pos in self.list:
			x, y = self.CalcScrolledPosition(pos)
			dc.SetPen(wx.Pen(wx.Colour(150,150,150), 6))
			dc.DrawCircle(x, y, 9)
			dc.SetPen(wx.Pen('black', 2))
			dc.DrawCircle(x, y, 9)
	def OnMove(self, e):
		pos = self.CalcUnscrolledPosition(e.GetPosition())
		self.parent.sb.SetStatusText(str(pos))
	def OnLClick(self, e):
		pos = self.CalcUnscrolledPosition(e.GetPosition())
		self.list.append(pos)
		self.OnPaint(-1)
	def OnRClick(self, e):
		if self.list != []:
			self.list.pop()


# Refresh is drawing static bitmap on top of circles.
		#self.parent.Refresh()

		#x, y = pos
		#dc = wx.MemoryDC()
		#dc.SelectObject(self.bitmap)
		#dc.SetBrush(wx.Brush(wx.Colour(255,255,255), style=wx.TRANSPARENT))
		#dc.SetPen(wx.Pen(wx.Colour(150,150,150), 6))
		#dc.DrawCircle(x, y, 9)
		#dc.SetPen(wx.Pen('black', 2))
		#dc.DrawCircle(x, y, 9)
		#self.statbit = wx.StaticBitmap(self.sw, -1, self.bitmap)


app = wx.PySimpleApp()
frame = MainFrame(None, wx.ID_ANY, 'Deflectometry')
app.MainLoop()
