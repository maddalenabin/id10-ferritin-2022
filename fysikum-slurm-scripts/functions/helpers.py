import yaml
import numpy as np
import smtplib, ssl
    
def dump_filelist(filelist, filename='filelist.yml', nprocs=16):
    files_per_proc = int(np.ceil(len(filelist)/nprocs))
    iproc = 0
    to_file = {}
    for i in range(len(filelist)):
        if i % files_per_proc == 0:
            if i > 0:
                iproc += 1
            to_file[iproc] = []

        to_file[iproc].append(filelist[i])


    with open(filename, 'w') as f:
        yaml.dump(to_file, f)


def email(body="Job done now", subject="slurm job done", filename="./confidential.txt"):
    """Send an email. You can edit the subject, body and the filename where the sender and receiver emails are with the password."""
    port = 587  # For starttls
    smtp_server = "smtp.office365.com"
    # load txt file where you write sender, receiver and password
    sender_email, receiver_email, password = np.loadtxt(filename, dtype='str')

    context = ssl.create_default_context()

    # -- this works
    message = f"""Subject: {subject}
    From: Maddalena <{sender_email}>
    To: <{receiver_email}>

    {body}"""

    with smtplib.SMTP(smtp_server, port) as server:
        server.ehlo()  # Can be omitted
        server.starttls(context=context)
        server.ehlo()  # Can be omitted
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email.split(','), message)