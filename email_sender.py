from sys import argv
from email.mime.text import MIMEText
import smtplib
import logging
logger = logging.getLogger('main')

def send_email(sender, receiver, subject, content):
    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = receiver
    s = smtplib.SMTP("mxout.tau.ac.il")
    s.send_message(msg)
    s.quit()

if __name__ == '__main__':
    if len(argv) < 5:
        logger.error('Usage: python ' + argv[0] + ' <sender> <receiver> <subject> <content>')
        exit()
    else:
        send_email(*argv[1:])
