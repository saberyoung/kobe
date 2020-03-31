"""############################################################################ 
2019/1/30 Start
define some functions for linking soft/email/etc
""" ############################################################################
from __future__ import print_function
from builtins import input
import os,sys

def sendemail(email,emailpass,emailsmtp,subj,\
              fromaddr,toaddrs,text,Imglist=[],Filelist=[]):
    """
    send email
    """
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.image import MIMEImage
    from email.mime.multipart import MIMEMultipart, MIMEBase
    from email import encoders

    try:
        # settings
        msg = MIMEMultipart()
        msg['Subject'] = subj
        msg['From'] = fromaddr
        msg['To'] = toaddrs

        # text
        text = MIMEText(text)
        msg.attach(text)

        # attachments
        for FileName in Filelist:
            part = MIMEBase('application', "octet-stream")
            part.set_payload(open(FileName, "rb").read())
            encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="%s"'%FileName)
            msg.attach(part)

        # image
        for ImgFileName in Imglist:
            img_data = open(ImgFileName, 'rb').read()
            image = MIMEImage(img_data, name=os.path.basename(ImgFileName))
            msg.attach(image)

        server = smtplib.SMTP(emailsmtp)
        server.starttls()
        server.login(email,emailpass)
        server.sendmail(fromaddr,toaddrs,msg.as_string())
        server.quit()
        return True
    except:
        return False

def createSSHClient(server, port, user, password):

    import paramiko
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client

def sftp_put_dir(server, port, user, password, filelist, remote_dir):

    import paramiko    
    t = paramiko.Transport(sock=(server, port))
    t.connect(username=user, password=password)
    sftp = paramiko.SFTPClient.from_transport(t)   

    try:sftp.chdir(remote_dir)  # Test if remote_path exists
    except IOError:sftp.mkdir(remote_dir)  # Create remote_path       

    for x in filelist:
        filename = os.path.split(x)[-1]
        remote_filename = remote_dir + '/' + filename
        print('Put %s...' % filename)
        sftp.put(x, remote_filename)
    sftp.close()

def wechat(_name,_msg,_img):
    # pip install -U wxpy
    # pip install -U wxpy -i "https://pypi.doubanio.com/simple/"
    from wxpy import Bot

    # initialize
    bot = Bot()

    # find friend
    if len(bot.friends().search(_name)) == 1:answ = 0
    else:
        for _nn,_friend in enumerate(bot.friends().search(_name)):print(_nn+1,_friend)
        answ = int(input('which one?'))
    my_friend = bot.friends().search(_name)[answ-1]
    
    # send msg
    if _msg: my_friend.send(_msg)
    # send img
    if _img: my_friend.send_image(_img)

def read_id(_token):
    #check bot id
    if sys.version_info>(3,0,0): 
        from slack import WebClient as SlackClient
    else: 
        from slackclient import SlackClient
    _idl,_idl1 = {}, {}
    slack_client = SlackClient(_token)       
    # user and channel list
    for _cc,_dd in zip(["users.list", "channels.list"], \
                       ['members', 'channels']):
        api_call = slack_client.api_call(_cc)        
        if api_call.get('ok'):
            # retrieve all users so we can find our bot
            users = api_call.get(_dd)
            for user in users:
                if 'name' in user:
                    _idl[user.get('name')] = user.get('id')
                    _idl1[user.get('id')] = user.get('name')
    return _idl,_idl1

def slack(token, channel, content, _files=[], _msg=True):
    # pip install slackclient
    if sys.version_info>(3,0,0): 
        from slack import WebClient as SlackClient
    else: 
        from slackclient import SlackClient
    slack_client = SlackClient(token)

    _idlist,_idrlist = read_id(token)
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
    return True

def phone(_account,_token,_from,_to,_txt):
    # pip install twilio
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
