#!/usr/bin/env python
from __future__ import print_function
from builtins import input

import os, sys, time, pst
if sys.version_info>(3,0,0): 
    from slack import WebClient as SlackClient
else: 
    from slackclient import SlackClient

_pstpath = pst.__path__[0]
_docfile = '%s/default/fmt/doc_slack.txt'%_pstpath

# define databse
cc = {}
cc['db'] = 'gw'
cc['host'] = 'localhost'
cc['user'] = 'syang'
cc['passwd'] = 'augusto90'

class handle_command(object):
    """
        Receives commands directed at the bot and determines if they
        are valid commands. If so, then acts on the commands. If not,
        returns back what it needs for clarification.
    """
    def __init__(self,token,channel,command):

        self.command = command

        # def output
        self.msg, self.files, self.link = '', [], []

        # parser command
        self.parse_cmd()

        # do tasks
        self.def_cmd()
        
        # output
        pst.slack(token, channel, self.msg, _files=self.files)

    def parse_cmd(self):
        import configparser

        # read healper from doc
        self.docp = {}
        _docs = configparser.ConfigParser()
        _docs.read(_docfile)
        for s in _docs.sections():
            self.docp[s] = {}
            for o in _docs.options(s):
                # check format
                _help,_cmd = _docs.get(s,o).split('//')
                self.docp[s][o] = [_help,_cmd]

    def def_cmd(self):

        self.command = self.command.strip()
        _cmd = self.command.split(':')[0].strip()
        if len(self.command.split(':'))==1:_dict = None
        else: _dict = self.command.replace('%s:'%_cmd,'').strip()
        _def = False
        for s in self.docp:
            for o in self.docp[s]:
                if _cmd == o:
                    eval(self.docp[s][o][1])
                    _def = True
        if not _def:
            self.msg += 'command %s not defined!\n\n'%_cmd
            self.helper_general()

    def helper(self,_dict):
        if _dict is None: # show all helper
            self.helper_general()
        else: # show specific task
            for _cmd in _dict.split(','):
                _in = False
                for s in self.docp:
                    for o in self.docp[s]:
                        if _cmd == o:
                            self.msg += '%s\n'%self.docp[s][o][0]
                            _in = True
                if not _in:
                    self.msg+='*%s* not defined\n\n'%_cmd.strip()

    def helper_general(self):
        for s in self.docp:
            self.msg += ' - *%s*: \t'%s
            for o in self.docp[s]:
                self.msg += '`%s`\t'%o
            self.msg += '\n\n'
        self.msg += "\n\n\nuse `help:command,command,...` "+\
                    "to see more detailed tutorials for specific command"

##################
def parse_dict(_dict):
    try:_dict = dict(e.strip().split('=') for e in _dict.split(';'))
    except:_dict=_dict
    return _dict

def url(_dict):
    _dict = parse_dict(_dict)
    print(_dict)

def cpar(_dict,_keys):
    _ok=False
    if not _ok:
        _okp=0
        _dict = parse_dict(_dict)
        for _kk in _keys:
            if not _kk in _dict:
                self.msg += "\n%s not defined, input:"%_kk
                command, channel, user = \
                        parse_slack_output(slack_client.rtm_read(),AT_BOT)
                if channel:
                    command = escape(command)
                else:
                    print ('..')
            else:
                _okp += 1
        if _okp==len(_keys):_ok=True

def conf(_dict):

    cpar(_dict, ['tel','class'])
    for _tel in _dict['tel'].split(','):
        if len(_tel) == 0:continue
        for _class in ['general','telescope']:  
            _config = pst.load_config(_tel)
            _d = False
            while not _d:
                answ = input('Modify %s-%s file? (Y/N/Q)'%(_tel,_class))
                if answ in ['y','Y']:
                    pst.choose(_tel,_class,_config)
                    _d = True
                elif answ in ['n','N']:
                    _d = True
                elif answ in ['q','Q']:
                    sys.exit()
                else:
                    print ('wrong option')

def parse_slack_output(slack_rtm_output,AT_BOT):
    """
        The Slack Real Time Messaging API is an events firehose.
        this parsing function returns None unless a message is
        directed at the Bot, based on its ID.
    """
    output_list = slack_rtm_output
    command,channel,user = None, None,None   
    if output_list and len(output_list) > 0:                
        for output in output_list:
            if output and 'text' in output and AT_BOT in output['text']: 
                command,channel,user = output['text'].split(AT_BOT)[1].strip(), \
                                       output['channel'],output['user']
    return command,channel,user

def escape(word):   
    replace_with = {
        '&gt;'  : '<',
        '&lt;'  : '>',
        '&amp;' : '&',
        '&quot;': '"', 
        '&#39'  : "'"}   
    for oo in replace_with: word = word.replace(oo,replace_with[oo])
    return word

################################   call to main #########################
def main():

    _config = pst.load_config()
    if not _config: sys.exit()

    # bot token and name
    token = _config['general']['params']['slack']['token']
    bname = _config['general']['params']['slack']['botname']

    if token=='' or bname=='':
        sys.exit('set slack infos in pstools.par')

    # read all user name and id
    _idlist,_idrlist = pst.read_id(token)
    input(_idlist)

    # get bot id
    if bname in _idlist:BOT_ID = _idlist[bname]
    else:sys.exit('bot name not in %s'%_idlist.keys())
    AT_BOT = "<@" + BOT_ID + ">"
    slack_client = SlackClient(token)

    # 1 second delay between reading from firehose
    READ_WEBSOCKET_DELAY = 0

    if slack_client.rtm_connect():
        print("Bot connected and running!")
        while True:
            try: 
                command, channel, user = \
                    parse_slack_output(slack_client.rtm_read(),AT_BOT)
            except Exception as e:
                pst.slack(token, channel, str(e))
            if channel:
                command = escape(command)
                print('usr:%s\tcmd:%s'%(_idrlist[user], command))

                # main code
                try:
                    handle_command(token, channel, command)
                except Exception as e:
                    pst.slack(token, channel, str(e))
            time.sleep(READ_WEBSOCKET_DELAY)
    else:
        print("Connection failed. Invalid Slack token or bot ID")
