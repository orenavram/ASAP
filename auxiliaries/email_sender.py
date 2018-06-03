from sys import argv
import logging
logger = logging.getLogger('main')


def send_email(smtp_server, sender, receiver, subject='', content=''):
    from email.mime.text import MIMEText
    from smtplib import SMTP
    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = receiver
    s = SMTP(smtp_server)
    s.send_message(msg)
    s.quit()

if __name__ == '__main__':
    if len(argv) < 4:
        logger.error('Usage: python ' + argv[0] + ' <smtp_server> <sender> <receiver> <?subject> <?content>')
        exit()
    else:
        send_email(*argv[1:])
