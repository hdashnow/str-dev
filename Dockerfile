FROM python:3.8

RUN git clone https://github.com/quinlan-lab/STRling.git

RUN cd STRling && mkdir bin && cd bin && \
    wget https://github.com/quinlan-lab/STRling/releases/download/v0.3.0/strling && \
    chmod +x strling

ENV PATH="/STRling/bin:${PATH}"

RUN strling -h
