#!/usr/bin/python
import sys
import smtplib
import os
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
import time
import getpass
from email import Encoders


def send_html(
    user, password, recipients, subject, body,server='smtp.gmail.com',port=465):
    """
    Send an email through smtp server
    Arguments:
    - `user`: sender email address
    - `password`: sender password
    - `recipients`: A list of email addresses
    - `subject`: Subject line
    - `body`: HTML body
    - `server`: smtp server
    - `port`: smtp port
    """
    bulk=2
    wait=1*60
    msg = MIMEMultipart('alternative')
    msg['Subject'] = subject
    msg['From'] = user
    part2 = MIMEText(body, 'html')
    msg.attach(part2)

    smtp = smtplib.SMTP_SSL(server,port)
    smtp.login(user.split('@')[0],password)
    count=0
    for addr in recipients:
        msg['To'] = addr
        smtp.sendmail(user,addr,msg.as_string())
        count+=1
        if count==bulk:
            count=0
            time.sleep(wait)

