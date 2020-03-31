from __future__ import print_function
from builtins import input
import os,sys,re,glob,datetime,pst,configparser
import pylab as pl
import healpy as hp
from pst.configure import *
from matplotlib.figure import Figure
import numpy as np
if sys.version_info>(3,0,0): 
   import tkinter as Tk
   from tkinter import filedialog, font
else: 
   import Tkinter as Tk
   import tkFileDialog as filedialog
   import tkFont as font

filetypes = ['general config',
             'tel config',
             'healpix fits',
             'voevent xml',
             'npz data']
actions = ['[1] Set configuration',
           '[2] Online alert',
           '[3] Offline alert',
           '[4] Documentation']
_nlmax = 20

###############################     MAIN WINDOW ##########################
class MainMenu(Tk.Frame):

   def __init__(self, parent):
      Tk.Frame.__init__(self, parent)
      self.parent = parent 
      self.datadir = 'local'
      # popup window
      self.w,self.h,self.x,self.y = 1200, 500, 50, 50
      # file menu
      # meg menu
      self.msgw, self.msgh = 30, 15
      self.initUI()

   def initUI(self):    ##########################################

      """ INIT Frames: divide Parent Menu """

      """ ROW 1: 2 spans: """
      # current dir
      self.wdirectory = Tk.Frame(self.parent)
      self.wdirectory.grid(row=0,column=0,columnspan=2)

      # file list keys: name, info
      self.tb = Tk.Label(self.wdirectory,\
            text='data directory: %s (%s)'%\
            (self.datadir,dirlist[self.datadir]))
      self.tb.grid(row=0,column=0,columnspan=2)

      """ ROW 2: file list """
      self.wlist = Tk.Frame(self.parent)
      self.wlist.grid(row=1,column=0,columnspan=2)

      """ ROW 3: Main menu"""
      self.mainmenu = Tk.Frame(self.parent)
      self.mainmenu.grid(row=2,column=0)

      ''' Left '''
      # action list
      self.actionFrame = Tk.LabelFrame(self.mainmenu, text="Actions")   
      self.actionFrame.grid(row=1,column=0,pady=3)
      self.actionMenu = Tk.Menubutton(self.actionFrame,text='options',width=10,
                                      height=2,activebackground='slate grey') 
      self.actionMenu.grid(row=0,column=0,pady=3)
      picks   = Tk.Menu(self.actionMenu)               
      self.actionMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in actions:
         picks.add_command(label=t,\
                  command=lambda k=t: self.defaction(k))

      # tel name...
      self.telFrame = Tk.LabelFrame(self.mainmenu, text="telescope name")
      self.telFrame.grid(row=1,column=1,pady=3,sticky=Tk.NW)
      self.telFrameTxt = Tk.Entry(self.telFrame, width=12, justify=Tk.CENTER)
      self.telFrameTxt.grid(row=0, column=0, columnspan=2, sticky=Tk.SW, pady=3)
      self.telFrameBut = Tk.Button(self.telFrame,text="submit",width=3,\
                                   command=self.ActionListFrame)
      self.telFrameBut.grid(row=0, column=1, columnspan=2, sticky=Tk.SE, pady=3)

      self.parparams = Tk.IntVar()
      self.generalButton = Tk.Radiobutton(self.telFrame, text="gen", \
                                          variable=self.parparams, value=1)
      self.telButton = Tk.Radiobutton(self.telFrame, text="tel", \
                                      variable=self.parparams, value=2)
      self.checkButton = Tk.Radiobutton(self.telFrame, text="check", \
                                      variable=self.parparams, value=3)
      self.generalButton.grid(row=1,column=0,columnspan=3,sticky=Tk.W)
      self.telButton.grid(row=1,column=1,columnspan=3,sticky=Tk.W)
      self.checkButton.grid(row=1,column=2,columnspan=3,sticky=Tk.E)

      # file type list
      self.fileFrame = Tk.LabelFrame(self.mainmenu, text="Check Files")   
      self.fileFrame.grid(row=2,column=0,pady=3)
      self.fileMenu = Tk.Menubutton(self.fileFrame, text='types',width=10,
                                    height=2,activebackground='slate grey') 
      self.fileMenu.grid(row=0,column=0,pady=3)
      picks   = Tk.Menu(self.fileMenu)               
      self.fileMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in filetypes+['NONE']:
         picks.add_command(label=t,\
                     command=lambda k=t: self.FileListFrame(k))

      # display, check, delete, switch buttons
      self.fmenu = Tk.LabelFrame(self.mainmenu, text="")  
      self.fmenu.grid(row=2,column=1,padx=3,columnspan=4,sticky=Tk.W)
      self.dB = Tk.Button(self.fmenu,text="Display",bg='pale green',\
         command= self.display,width=4,height=2,activebackground='medium sea green')
      self.dB.grid(row=0,column=0,padx=0,sticky=Tk.NW)

      self.cB = Tk.Button(self.fmenu,text="Check",bg='light blue',\
         command= self.check,width=4,height=2,activebackground='medium sea green')
      self.cB.grid(row=1,column=0,padx=0,sticky=Tk.SW)

      self.eB = Tk.Button(self.fmenu,text="Remove",bg='light yellow',\
         command= self.remove,width=4,height=2,activebackground='medium sea green')
      self.eB.grid(row=0,column=1,padx=0,sticky=Tk.NE)

      self.eB = Tk.Button(self.fmenu,text="Switch",bg='misty rose',\
         command= self.switch,width=4,height=2,activebackground='medium sea green')
      self.eB.grid(row=1,column=1,padx=0,sticky=Tk.SE)

      # directory list
      self.dirmenu = Tk.LabelFrame(self.mainmenu, text="Change dir")  
      self.dirmenu.grid(row=3,column=0,padx=3,columnspan=2,sticky=Tk.W)
      self.dirMenuBut = Tk.Menubutton(self.dirmenu, text='Typ',width=2,
                                   height=1,activebackground='slate grey') 
      picks   = Tk.Menu(self.dirMenuBut)
      self.dirMenuBut.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in dirlist.keys():
         picks.add_command(label=t,\
                     command=lambda k=t: self.DirListFrame(k))
      self.dirMenuBut.grid(row=0,column=0,pady=3)
      # change dir
      self.folder_path = Tk.StringVar()
      button2 = Tk.Button(self.dirmenu, text="Browse", width=5,
                        height=1, command=self.browse_button)
      button2.grid(row=0, column=1)

      # quit button
      self.quitBut = Tk.Button(self.mainmenu,text="quit",fg="red",\
                     command=self.quit,width=6,height=2,\
                     activebackground='slate grey')
      self.quitBut.grid(row=3,column=1,sticky=Tk.W)

      ''' Right '''
      # message menu
      self.wmess = Tk.Frame(self.parent)
      self.wmess.grid(row=2,column=1)

      """ update values """
      # select files
      self.nn = 0

      ###  list windows
      self.sblist = Tk.Scrollbar(self.wlist)
      self.sblist.pack(side='right', fill='y')
      self.lbox = Tk.Listbox(self.wlist,width=70,height=15,\
               yscrollcommand=self.sblist.set,selectmode='extended')

      ###  message windows
      self.sbmess = Tk.Scrollbar(self.wmess)
      self.sbmess.pack(side='right', fill='y')
      self.textbox = Tk.Text(self.wmess,width=self.msgw,height=self.msgh,\
                  yscrollcommand=self.sbmess.set)
      self.textbox.config(foreground = "red")
      self.textbox.insert(Tk.END,6*"\n"+"!!! WATCH FOR MESSAGES HERE !!!")
      self.sbmess.config(command=self.textbox.yview)
      self.textbox.config(state='disabled')
      self.textbox.pack()

   def browse_button(self):
      self.filename = filedialog.askdirectory()
      self.folder_path.set(self.filename)
      print(self.filename)

   ##########################################
   # var: actparams, tel name, gen/tel, submit
   def defaction(self,_type): 
      self.actparams = _type
      self.Message(self.actparams)

   def ActionListFrame(self):
      if not hasattr(self,'actparams'): 
         self.Message('choose actions first...')
         return
      if self.actparams == actions[0]:
         self.setconfig()
      if self.actparams == actions[1]:
         self.setalert()
      if self.actparams == actions[2]:
         self.offalert()
      if self.actparams == actions[3]:
         self.tutorial()

   ''' 1 config part '''
   def setconfig(self): # action: 0

      _tel = self.telFrameTxt.get()
      if _tel == '': 
         self.Message('set telescope name first')
         return
      if self.parparams.get()==0:
         self.Message('set par type first')
         return
      if self.parparams.get() == 1:_class='general'
      elif self.parparams.get() == 2: _class = 'telescope'
      else:
         _msg=''
         for _tel1 in _tel.split(','):
            _msg+='for %s:\n'%_tel1
            _config = pst.load_config(_tel1)
            for _cl in ['general','telescope']:
               _msg += '- %s: %s \n'%(_cl,_config[_cl]['config'])
         self.Message(_msg)
         return

      _flist = {'general':[],'telescope':[]}
      for _tel1 in _tel.split(','):
         # for one telescope
         _config = pst.load_config(_tel1,_indir=dirlist[self.datadir])
         _flist[_class].append(os.path.basename(_config[_class]['config']))
         if hasattr(self,'filelist'):
            self.filelist[_config[_class]['config']] = _config[_class]['check']

      # update list
      if hasattr(self, 'typefile'):
         self.filelist = group_data(dirlist[self.datadir],self.typefile)
      self.lbox.delete(0,Tk.END)
      if hasattr(self, 'filelist'):
         for _name in self.filelist:
            _info = self.filelist[_name]
            imgtext = '%-45s       %-5s' %\
                      (os.path.basename(_name)[:45],_info)
            self.lbox.insert(Tk.END,imgtext)

      # pop up
      _msg = ''
      for _class in _flist: 
         for _ff in np.unique(_flist[_class]):
            _msg+='%s\n'%_ff
            self.popup(_ff,_class)
      self.Message(_msg)

   def popup(self,_configfile,_class): 

      # Init frame
      self.wpop = Tk.Toplevel()
      self.wpop.title('%s dir: %s (%s)'%(self.datadir,_configfile,_class))
      self.wpop.geometry("%sx%s+%s+%s"%(self.w,self.h,self.x,self.y))
      Tk.Label(self.wpop, text="Revise config: \n")

      # read docs
      _docs = configparser.ConfigParser()
      _docs.read(filelist[_class]['doc'])

      # read config
      config = configparser.ConfigParser()
      config.read('%s/%s'%(dirlist[self.datadir],_configfile))

      # go to element
      self.varlist = {}
      for s in config.sections():
         _varFrame = Tk.LabelFrame(self.wpop, text="%s"%s)
         self.varlist[s] = {}
         for o in config.options(s):
            _help,_fmt,_opt = _docs.get(s,o).replace('\n','').split('//')
            _varFrame1 = Tk.LabelFrame(_varFrame, text="%s"%o)
            # entry
            if config.get(s,o) != '': v = Tk.StringVar(value=config.get(s,o))
            else: v = Tk.StringVar(value=_opt.replace(' ',''))
            _ent = Tk.Entry(_varFrame1, textvariable=v, width=5, justify=Tk.CENTER)
            _ent.pack(side=Tk.LEFT,expand=True)
            self.varlist[s][o] = _ent

            # button
            _but = Tk.Button(_varFrame1,text="Help",\
                             command=lambda k=[s,o,_class]: self.helpVar(k),width=4,height=1)
            _but.pack(side=Tk.LEFT,expand=True)
            _varFrame1.pack(side=Tk.LEFT,expand=True,anchor=Tk.W,padx=10)
         _varFrame.pack(side=Tk.TOP,expand=True,anchor=Tk.W)
      Tk.Button(self.wpop, text='submit', \
                command=lambda k=[_configfile,_class]: self.updateVar(k)).pack(side=Tk.LEFT,expand=True)
      Tk.Button(self.wpop, text='close', command=self.wpop.destroy).pack(side=Tk.LEFT,expand=True)

   def updateVar(self,_ll):

      _configfile, _class = _ll
      if hasattr(self, 'varlist'):
         self.Message('update: %s/%s\n'%\
                      (dirlist[self.datadir],_configfile))
         _mdir = []
         for s in self.varlist:
            for o in self.varlist[s]:
               if len(self.varlist[s][o].get())>0:
                  _mdir.append(([s,o],self.varlist[s][o].get()))
         # modify
         pst.modify_config('%s/%s'%(dirlist[self.datadir],_configfile),_mdir)

         # remain pos and size
         self.w = self.wpop.winfo_width()
         self.h = self.wpop.winfo_height()
         self.x = self.wpop.winfo_x()
         self.y = self.wpop.winfo_y()
         self.wpop.destroy()
         self.popup(_configfile,_class)

   def helpVar(self,_ll):

      if hasattr(self,'varlist'):
         # read docs
         _docs = configparser.ConfigParser()
         _docs.read(filelist[_ll[2]]['doc'])
         _help,_fmt,_opt = _docs.get(_ll[0],_ll[1]).replace('\n','').split('//')
         self.Message('Help:%s\n\nFormat:%s\n\ndefault:%s'%(_help,_fmt,_opt))

   ''' serve alert '''
   def setalert(self):

      # serve alert
      self.saFrame = Tk.LabelFrame(self.mainmenu, text="Alert actions")   
      self.saFrame.grid(row=1,column=2,pady=3)
      self.saMenu = Tk.Menubutton(self.saFrame,text='options',width=10,
                                      height=2,activebackground='slate grey') 
      self.saMenu.grid(row=0,column=0,pady=3)
      picks   = Tk.Menu(self.saMenu)               
      self.saMenu.config(menu=picks,relief=Tk.RAISED,bd=5)  
      for t in ['serve local','monitor local','monitor LVC']:
         picks.add_command(label=t,\
                  command=lambda k=t: self.AlertListFrame(k))

   def AlertListFrame(self, _action):

      from gcn.cmdline import serve_main
      import gcn

      if 'local' in _action:
         _s1,_s2 = pst.lvc_server('local', 'PUBLIC', True)
      elif 'LVC' in _action:
         _s1,_s2 = pst.lvc_server('eApps', 'PUBLIC', True)
      else:
         self.Message('wrong server/port')
         return

      if _action == 'serve local':
         if hasattr(self,'nextfile'):
            if len(self.nextfile)>0:
               for _xml in self.nextfile:
                  self.Message('Serve Mode:\thost:%s:port:%s'%(_s1,_s2))
                  if '.xml' in _xml:
                     serve_main(args=['%s/%s'%(dirlist[self.datadir],_xml),\
                                      '--host', '%s:%s'%(_s1,_s2)])
         self.Message('select voevent-xml file first...')
      else:
         self.Message('Monitor Mode:\thost:%s:port:%s'%(_s1,_s2))
         gcn.listen(handler = pst.process_gcn, host=_s1, port=_s2)

   ''' offline serch '''
   def offalert(self):

      _tel = self.telFrameTxt.get()
      if _tel == '': 
         self.Message('set telescope name first')
         return
      if hasattr(self,'nextfile'):
         if len(self.nextfile)>0:
            for _fits in self.nextfile:
               if '.fits' in _fits or '.xml' in _fits:
                  _fits = '%s/%s'%(dirlist[self.datadir],_fits)
                  self.Message('Offline search for %s'%_fits)
                  pst.man_search(_fits, _tel, True)
      pst.man_search(None, _tel, True)

   ''' tutorial '''
   def tutorial(self):
      _msg = 'github: https://github.com/saberyoung/pstools\n\n'
      _msg += 'tutorial: https://pstools-documentation.readthedocs.io/en/latest/\n\n'
      _msg += 'GraceDB: https://gracedb.ligo.org/superevents/'
      self.Message(_msg)

   ##########################################
   def FileListFrame(self,_type):

      # group files
      self.typefile = _type
      self.filelist = group_data(dirlist[self.datadir],self.typefile)

      # show keys
      imgtext = '%-50s        status OK?'%(_type[:50])
      textboxb =  Tk.Label(self.wdirectory,text=imgtext)
      textboxb.grid(row=1,column=0,sticky='W')

      # show values
      self.lbox.delete(0,Tk.END)
      for _name in self.filelist:
         _info = self.filelist[_name]
         imgtext = '%-55s       %-5s' % (os.path.basename(_name)[:55],_info)
         self.lbox.insert(Tk.END,imgtext)

      self.lbox.pack(side=Tk.LEFT,expand=True, fill=Tk.BOTH)
      self.sblist.config(command=self.lbox.yview)

      self.nextfile = []
      def productselect(*event):
         nn = self.lbox.curselection()
         self.nextfile = [self.lbox.get(n).split()[0] for n in nn]  
         try:self.nn = nn[0]
         except:self.nn = 0
      self.lbox.bind("<<ListboxSelect>>", productselect)      


   def DirListFrame(self,_dir):   ##########################################

      self.datadir = _dir
      # update list
      if hasattr(self, 'typefile'):
         self.filelist = group_data(dirlist[self.datadir],self.typefile)
         self.lbox.delete(0,Tk.END)
         for _name in self.filelist:
            _info = self.filelist[_name]
            imgtext = '%-45s       %-5s' % (os.path.basename(_name)[:45],_info)
            self.lbox.insert(Tk.END,imgtext)

      # update tb
      self.tb.destroy()
      self.tb = Tk.Label(self.wdirectory,\
            text='data directory: %s (%s)'%\
            (self.datadir,dirlist[self.datadir]))
      self.tb.grid(row=0,column=0,columnspan=2)

   def Message(self,message):    #########################################

      self.textbox.config(state='normal',foreground='black')
      self.textbox.delete(1.0,Tk.END)
      self.textbox.insert(Tk.END,message)
      self.sbmess.config(command=self.textbox.yview)
      self.textbox.config(state='disabled')

   def display(self):   #########################################
      # 'general config','tel config','healpix fits','voevent xml','npz data'
      
      if not hasattr(self,'nextfile') or len(self.nextfile)==0:
         self.Message('select files first...')
         return

      _msg = ''
      for _ff in self.nextfile:
         _ff = '%s/%s'%(dirlist[self.datadir],_ff)
         if 'config' in self.typefile:
            _msg += '-'*5
            _msg += '%s:\n'%_ff
            # read config
            config = configparser.ConfigParser()
            config.read(_ff)
            for s in config.sections():
               _msg += '[%s]\n'%s
               for o in config.options(s):
                  _msg += '- %s:\t%s\n'%(o,config.get(s,o))
            _msg += '\n'

         if 'fits' in self.typefile:         
            pl.figure()
            _map = hp.read_map(_ff)
            hp.mollview(map=_map, fig=1)
            pl.show()

         if 'npz' in self.typefile:
            _msg += '-'*5
            _msg += '%s (show %i lines):\n'%(_ff,_nlmax)
            # load
            _dict = np.load(_ff)
            _length,_sl,_kl = [],'\n',[]
            for _key in _dict.keys():
               try:
                  _length.append(len(_dict[_key]))
                  _msg += '%s\t'%_key
                  _kl.append(_key)
               except:
                  _sl += '%s:%s\t'%(_key,_dict[_key])
            _msg += _sl
            _msg += '\n'            
            for ii in range(max(_length)):
               if ii>_nlmax:continue
               for _key in _kl:
                  try:
                     _msg += '%s\t'%_dict[_key][ii]
                  except:pass
               _msg += '\n'
            _msg += '\n'
      self.Message(_msg)

   def check(self):
      # 'general config','tel config','healpix fits','voevent xml','npz data'

      if not hasattr(self,'nextfile') or len(self.nextfile)==0:
         self.Message('select files first...')
         return
      _msg,_c = '',None
      for _ff in self.nextfile:
         _dir = dirlist[self.datadir]
         if 'config' in self.typefile: 
            if 'pstools.par' in _ff:
               _class = 'general'
               _config = pst.load_config(_inclass=_class, \
                                 _indir=_dir, _verbose=True)
               _c = _config['general']['check']
            else:
               _class = 'telescope'
               _tel = _ff.replace('pst_','').replace('.par','')
               _config = pst.load_config(_tel=_tel, _inclass=_class, \
                                         _indir=_dir, _verbose=True)
               _c = _config['telescope']['check']
         if 'fits' in self.typefile:
            try:              
               hp.read_map('%s/%s'%(_dir,_ff))
               _c = True
            except:
               _c = False
         if 'npz' in self.typefile:
            try:              
               np.load('%s/%s'%(_dir,_ff))
               _c = True
            except:
               _c = False

         self.filelist[_ff] = _c
         _msg += '%s/%s: %s\n'%(_dir,_ff,_c)
      self.Message(_msg)               

      # update list
      self.lbox.delete(0,Tk.END)
      for _name in self.filelist:
         _info = self.filelist[_name]
         imgtext = '%-45s       %-5s' % (_name[:45],_info)
         self.lbox.insert(Tk.END,imgtext)

   def remove(self):
      # 'general config','tel config','healpix fits','voevent xml','npz data'
      
      if not hasattr(self,'nextfile') or len(self.nextfile)==0:
         self.Message('select files first...')
         return
      _msg=''
      for _ff in self.nextfile: 
         os.remove('%s/%s'%(dirlist[self.datadir],_ff))
         del self.filelist[_ff]
         _msg += 'remove %s'%('%s/%s'%(dirlist[self.datadir],_ff))
      self.Message(_msg)

      # update list
      self.lbox.delete(0,Tk.END)
      for _name in self.filelist:
         _info = self.filelist[_name]
         imgtext = '%-45s       %-5s' % (_name[:45],_info)
         self.lbox.insert(Tk.END,imgtext)

   def switch(self):
      # switch dir
      if not hasattr(self,'nextfile') or len(self.nextfile)==0:
         self.Message('select files first...')
         return
      _msg = ''
      for _ff in self.nextfile: 
         for _dirn in dirlist:
            if _dirn == self.datadir: 
               _currentdirn = _dirn
               _currentdir = dirlist[_dirn]
            else: 
               _targetdirn = _dirn
               _targetdir = dirlist[_dirn]
         os.rename('%s/%s'%(_currentdir,_ff),'%s/%s'%(_targetdir,_ff))
         _msg += 'move %s to %s'%('%s/%s'%(_currentdir,_ff),'%s/%s'%(_targetdir,_ff))
      self.Message(_msg)

      # update
      self.DirListFrame(_targetdirn)

   def quit(self):    #########################################

      self.parent.quit()
      self.parent.destroy()

def group_data(data_dir, _type):
   _filelist = {}
   for _infile in glob.glob(data_dir+'/*'):
      _file = os.path.basename(_infile) 
      _dir = '%s/'%os.path.dirname(_infile)
      if 'pstools.par' in _file:
         if _type == 'general config':
            _config = pst.load_config(_inclass='general', _indir=_dir)
            _filelist[_file] = _config['general']['check']
      elif re.match("pst_.*.par", _file):
         if _type == 'tel config':
            _tel = _file.replace('pst_','').replace('.par','')
            _config = pst.load_config(_tel=_tel, _indir=_dir)
            _filelist[_file] = _config['telescope']['check']
      elif '.fits' in _file:
         if _type == 'healpix fits':                 
            _filelist[_file] = None
      elif '.xml' in _file:
         if _type == 'voevent xml':                  
            _filelist[_file] = None
      elif '.npz' in _file:
         if _type == 'npz data':
            try:
               import numpy as np
               np.load(_file)
               _c = True
            except:
               _c = False
            _filelist[_file] = _c
      else:
         if _type == 'NONE':
            _filelist[_file] = None
   return _filelist

################################   call to main #########################
def main():
   root = Tk.Tk()
   root.geometry("+50+50")
   root.title('PSTools GUI (version '+pst.__version__+')')
   default_font = font.nametofont("TkFixedFont")
   default_font.configure(weight="bold")
   root.option_add( "*font", default_font)
   root.protocol('WM_DELETE_WINDOW', root.quit)
   app = MainMenu(root)
   app.mainloop()
