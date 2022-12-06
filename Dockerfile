FROM matthewfeickert/docker-python3-ubuntu
WORKDIR /daniotalk
COPY Data/fetch_data.sh /daniotalk/Data/
COPY Assets/* /daniotalk/Assets/
COPY create_pairs.py /daniotalk/

WORKDIR /daniotalk/Data/
RUN sudo ./fetch_data.sh
RUN python -m pip install pandas sqlalchemy openpyxl
RUN mkdir -p /daniotalk/Database
WORKDIR /daniotalk

