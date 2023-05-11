FROM matthewfeickert/docker-python3-ubuntu
WORKDIR /daniotalk
COPY Data/fetch_data.sh /daniotalk/Data/
COPY requirements.txt /daniotalk/requirements.txt
COPY Assets/* /daniotalk/Assets/
COPY create_pairs.py /daniotalk/

WORKDIR /daniotalk/Data/
RUN sudo bash fetch_data.sh
WORKDIR /daniotalk
RUN sudo python -m pip install -r requirements.txt
RUN mkdir -p /daniotalk/Database

