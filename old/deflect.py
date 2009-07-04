import wx
import os
import profile

ID_OPEN=111
ID_SOLVE=101
ID_EXIT=110

# "http://www.zetcode.com/wxpython/tips/" for how to write config file.

class MainFrame(wx.Frame):
	""" We simply derive a new class of Frame. """
	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title, size=(800,700))
		#self.Maximize()

		self.files = {}
		self.outpath = ''
		self.dots = {}
		self.calibfile = ''
		self.surfaces = {}

		self.pan = wx.Panel(self, -1)
		self.listb = wx.ListBox(self.pan, -1, size=(150, 600))
		self.imwin = ImageWindow(self.pan, -1)
		self.radh = wx.RadioButton(self.pan, -1, 'Horizontal',
				style=wx.RB_GROUP)
		self.radv = wx.RadioButton(self.pan, -1, 'Vertical')

		self.hsiz = wx.BoxSizer(wx.HORIZONTAL)
		self.vsiz = wx.BoxSizer(wx.VERTICAL)
		self.vsiz.Add(self.radh, 0, wx.ALL)
		self.vsiz.Add(self.radv, 0, wx.ALL)
		self.vsiz.Add(self.listb, 0, wx.EXPAND | wx.TOP, 5)
		self.hsiz.Add(self.vsiz, 0, wx.EXPAND | wx.ALL, 5)
		self.hsiz.Add(self.imwin, 1, wx.EXPAND | wx.ALL)

		self.sb = self.CreateStatusBar()
		# Create menus.
		filemenu = wx.Menu()
		filemenu.Append(ID_OPEN, '&Open', 'Open file')
		filemenu.Append(ID_SOLVE, '&Solve', 'Solve profiles')
		filemenu.AppendSeparator()
		filemenu.Append(ID_EXIT, 'E&xit', 'Terminate the program')
		#exit = wx.MenuItem(filemenu, ID_EXIT, 'E&xit', 'Terminate the program')
		#exit.SetBitmap(wx.Bitmap('icons/exit.png'))
		#filemenu.AppendItem(exit)
		# Create menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(filemenu, "&File")
		self.SetMenuBar(menuBar)
		self.Bind(wx.EVT_MENU, self.OnOpen, id=ID_OPEN)
		self.Bind(wx.EVT_MENU, self.OnExit, id=ID_EXIT)
		self.Bind(wx.EVT_MENU, self.OnSolve, id=ID_SOLVE)
		self.Bind(wx.EVT_RADIOBUTTON, self.OnSelect)

		self.Bind(wx.EVT_LISTBOX, self.OnSelect)
		self.pan.Bind(wx.EVT_PAINT, self.imwin.OnPaint)
		self.pan.SetSizer(self.hsiz)
		self.Show(True)
		
	def OnOpen(self, e):
		"""Open a file"""
		dlg = wx.FileDialog(self, "Choose a file", self.outpath,
				"", "*.*", wx.MULTIPLE)
		if dlg.ShowModal() == wx.ID_OK:
			fns = dlg.GetFilenames()
			# Filenames as "panel_ID-h.tiff" or "panel_ID-v.tiff".
			path = dlg.GetDirectory()
			self.outpath = path
			for fn in fns:
				id = fn[0:fn.rfind('.')-2]
				ori = fn[fn.rfind('.')-1]
				ind = -1
				if (ori == 'h') | (ori == 'H'):
					ind = 0
				elif (ori == 'v') | (ori == 'V'):
					ind = 1
				else:
					print('Invalid File Name')
				if self.files.has_key(id):
					self.files[id][ind] = os.path.join(path, fn)
				else:
					pair = [-1, -1]
					pair[ind] = os.path.join(path, fn)
					self.files[id] = pair
				if not(self.dots.has_key(id)):
					self.dots[id] = []
		dlg.Destroy()
		#self.ids = self.files.keys()
		self.listb.SetItems(self.files.keys())
		if self.files != {}:
			self.listb.SetSelection(0)
			self.OnSelect(self.listb)
	def OnSolve(self, e):
		dotsok = True
		for id in self.dots.keys():
			if len(self.dots[id]) != 3:
				dotsok = False
		if dotsok:
			for id in self.files.keys():
				self.surfaces[id] = profile.findprofile(self.files[id],
					self.dots[id], self.outpath, self.calibfile)
		else:
			print('Not all dots have been located')
	def OnExit(self, e):
		self.Close(True)
	def OnSelect(self, e):
		ind = 0
		if self.radv.GetValue():
			ind = 1
# For some reason this is getting called on opening a new file where the
# selected item indice is -1. OnSelect gets called again and works second time.
# To replicate bug open one image, and then a second.
		self.imwin.setimage(self.files[self.listb.GetStringSelection()][ind])
	
class ImageWindow(wx.ScrolledWindow):
	""" Window container for bitmap. """
	def __init__(self, parent, id):
		wx.ScrolledWindow.__init__(self, parent, id)
		self.parent = parent

		self.Bind(wx.EVT_MOTION, self.OnMove)
		self.Bind(wx.EVT_LEFT_UP, self.OnLClick)
		self.Bind(wx.EVT_RIGHT_UP, self.OnRClick)
		#self.Bind(wx.EVT_SCROLLWIN, self.OnScroll)
	
	def setimage(self, filename):
		self.bitmap = wx.Bitmap(filename, type=wx.BITMAP_TYPE_TIF)
		sz = self.bitmap.GetSize()
		incr = 40
		self.SetScrollbars(incr, incr, sz[0]/incr+1, sz[1]/incr+1)

	def OnPaint(self, e):
		dc = wx.PaintDC(self)
		dc.Clear()
		id = self.parent.GetParent().listb.GetStringSelection()
		if id != '':
			x, y = self.CalcScrolledPosition((0, 0))
			dc.DrawBitmap(self.bitmap, x, y)
			#dc.SetBrush(wx.Brush(wx.Colour(255,255,255), style=wx.TRANSPARENT))
			for pos in self.parent.GetParent().dots[id]:
				x, y = self.CalcScrolledPosition(pos)
				dc.SetPen(wx.Pen('black', 2))
				length = 10
				dc.DrawLine(x-5, y, x-5-length, y)
				dc.DrawLine(x+5, y, x+5+length, y)
				dc.DrawLine(x, y-5, x, y-5-length)
				dc.DrawLine(x, y+5, x, y+5+length)
	def OnMove(self, e):
		pos = self.CalcUnscrolledPosition(e.GetPosition())
		self.parent.GetParent().sb.SetStatusText(str(pos))
	def OnLClick(self, e):
# Must select dots starting at top left and going anti-clockwise.
		id = self.parent.GetParent().listb.GetStringSelection()
		if id != '':
			pos = self.CalcUnscrolledPosition(e.GetPosition())
			if len(self.parent.GetParent().dots[id]) < 3:
				self.parent.GetParent().dots[id].append(pos)
			self.OnPaint(self)
	def OnRClick(self, e):
		id = self.parent.GetParent().listb.GetStringSelection()
		if id != '':
			if self.parent.GetParent().dots[id] != []:
				self.parent.GetParent().dots[id].pop()
			self.OnPaint(self)
	def OnScroll(self, e):
# Takes over parent event handler.
		print('scrolling')
		self.OnPaint(self)


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
