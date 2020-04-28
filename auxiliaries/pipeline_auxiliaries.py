import os

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


def measure_time(total):
    hours = total // 3600
    minutes = (total % 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes}:{seconds} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds} minutes'
    else:
        return f'{seconds} seconds'


def create_dir(path):
    if not os.path.exists(path):
        import logging
        logger = logging.getLogger('main')
        logger.info("Creating directory:\n" + path)
        os.makedirs(path)