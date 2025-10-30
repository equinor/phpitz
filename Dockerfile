FROM intel/fortran-essentials
COPY . work
WORKDIR /work

RUN cd src && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_EXE_LINKER_FLAGS="-qmkl" /work && \
    make

RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt update && \
    apt install -y python3.12 python3.12-venv

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN python3.12 -m venv /opt/server_env
ENV PATH="/opt/server_env/bin:$PATH"

RUN cd api && \
    pip install .

ENV RUN_FOLDER=/work/src/build
USER 1001
EXPOSE 8008

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8008"]
