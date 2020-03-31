#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/circulate/circulate.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
from kobe import utils

__all__ = ['circulate']

class circulate(utils):
    """provide various functions to circulate KOBE products.
    inherient from kobe.utils;
    will be inheriented by kobe.trigger/galaxies/tilings/candidates.  
    """
    
    defkwargs = {'emailpass'    : '',
                 'emailsmtp'    : '',
                 'subject'      : 'GW Alerts',
                 'fromaddr'     : '',
                 'toaddrs'      : '',                 
                 'slacktoken'   : None,
                 'slackto'      : None,                
    }        
    '''
    default parameter settings for: email, slack, esop2 and so on.
    '''
    defkwargs = {**defkwargs,
                 **utils.defkwargs}

    texts = ''
    '''
    texts to be circulated
    '''
    
    attachments = []
    '''
    attachments to be circulated
    '''
    
    images = []
    '''
    iamges to be circulated
    '''
    
    def send_email(self, **kwargs):
        """send out KOBE products via email
        """
        kwargs = self.setkeys(kwargs)
        
        import smtplib
        from email.mime.text import MIMEText
        from email.mime.image import MIMEImage
        from email.mime.multipart import MIMEMultipart, MIMEBase
        from email import encoders
        from email.header import Header
        
        sender = kwargs['fromaddr']
        if type(kwargs['toaddrs']) is list:
            receivers = kwargs['toaddrs']
        elif type(kwargs['toaddrs']) is str:
            receivers = kwargs['toaddrs'].split(',')
        else:
            return
        
        # settings
        msg = MIMEMultipart()
        for ii,jj in zip(['Subject', 'From', 'To'],
                         ['subject','fromaddr', 'toaddrs']):
            assert type(kwargs[jj]) is str, '%s should be a string'%jj
            msg[ii] = Header(kwargs[jj], 'utf-8')            

        # text
        text = MIMEText(self.texts, 'plain', 'utf-8')
        msg.attach(text)

        # attachments
        for FileName in self.attachments:
            part = MIMEBase('application', "octet-stream")
            part.set_payload(open(FileName, "rb").read())
            encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"'%FileName)
            msg.attach(part)

        # image
        for ImgFileName in self.images:
            img_data = open(ImgFileName, 'rb').read()
            image = MIMEImage(img_data, name=os.path.basename(ImgFileName))
            msg.attach(image)

        # send
        try:
            smtpObj = smtplib.SMTP(kwargs['emailsmtp'])
            smtpObj.starttls()
            smtpObj.login(kwargs['fromaddr'],kwargs['emailpass'])
            smtpObj.sendmail(sender, receivers, msg.as_string())
            smtpObj.quit()
        except smtplib.SMTPException as e:
            self.logger.info (e)                
        
    def _read_slackid(self, **kwargs):        
        if self.ckpython() == 3:
            from slack import WebClient as SlackClient
        else: 
            from slackclient import SlackClient
            
        #check bot id
        kwargs = self.setkeys(kwargs)
        _token = kwargs['token']
        if _token is None:
            self.logger.info ('Warning: set token first')
            return
        else:
            slack_client = SlackClient(_token)
            
        _idl = {}        
        for _list,_att in zip(["users.list", "channels.list"], \
                              ['members', 'channels']):
            api_call = slack_client.api_call(_list)            
            if api_call.get('ok'):
                # retrieve all users so we can find our bot            
                for user in api_call.get(_att):
                    if 'name' in user:
                        _idl[user.get('name')] = user.get('id')
        return _idl

    def send_slack(**kwargs):
        """send out KOBE products via slack
        """
        self.ckpython()
        if self.pythonversion == 3:
            from slack import WebClient as SlackClient
        else: 
            from slackclient import SlackClient
            
        # self.texts, self.attachments = [], []
        kwargs = self.setkeys(kwargs)
        _token = kwargs['token']
        if _token is None:
            self.logger.info ('Warning: set token first')
            return
        else:
            slack_client = SlackClient(_token)
            
        _idlist = self._read_slackid()
        if channel in _idlist: 
            # channel is name, need to find id
            channel = _idlist[channel]
        if _msg:
            if sys.version_info>(3,0,0):
                slack_client.api_call("chat.postMessage", \
                                      json={'channel':channel,'text':content})
            else:
                slack_client.api_call("chat.postMessage", channel=channel,
                                      text=content, as_user=True)
        for _file in _files:
            if sys.version_info>(3,0,0):
                slack_client.files_upload(channels=channel,file=_file)
            else:
                slack_client.api_call("files.upload",channels=channel,
                                      file=open(_file, 'rb'),filename=_file)
                
    def send_sms(_account,_token,_from,_to,_txt):
        """send out KOBE products via sms (twilio package needed)
        """
        try:from twilio.rest import Client
        except: return False
        
        account_sid = _account
        auth_token  = _token
        client = Client(account_sid, auth_token)
        message = client.messages.create(
            to=_to,
            from_=_from,
            body=_txt)
        print(message.sid)
        return True
