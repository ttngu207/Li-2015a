FROM datajoint/pydev:python3.6

ENV PATH "$PATH:/root/.local/bin"

RUN pip install --upgrade pip setuptools

ADD . /home/Li2015
RUN pip install --user -e /home/Li2015
RUN pip install --user -r /home/Li2015/requirements.txt

COPY ./jupyter_notebook_config.py /etc/jupyter/jupyter_notebook_config.py
