FROM ubuntu:16.04

# Make Kovah and main script available
COPY include/* /tmp/Kavosh
COPY MC_script_gt_v2.py /tmp/MC_script_gt_v2.py

# Install Kavosh
RUN apt-get update \
&& apt install -y aptitude \
&& apt-key adv --keyserver pgp.skewed.de --recv-key 98507F25 \
&& echo 'deb http://downloads.skewed.de/apt/xenial xenial universe' | tee -a  /etc/apt/sources.list \
&& echo 'deb-src http://downloads.skewed.de/apt/xenial xenial universe' | tee -a  /etc/apt/sources.list \
&& apt-get update \
&& aptitude -y install \
g++ \
make \
python3 \
python3-pip \
python3-graph-tool \
&& rm -rf /var/lib/apt/lists/* \
&& pip3 install networkx \
&& cd /tmp/Kavosh \
&& make \
&& cp Kavosh /usr/local/bin/Kavosh \
&& cp /tmp/MC_script_gt_v2.py /usr/local/bin/MC_script_gt_v2 \
&& chmod +x /usr/local/bin/MC_script_gt_v2
ENTRYPOINT ["python3", "/usr/local/bin/MC_script_gt_v2"]
CMD ["-h"]
