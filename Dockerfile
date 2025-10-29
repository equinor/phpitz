FROM intel/fortran-essentials as build-env
COPY src src
WORKDIR /src

RUN mkdir build
RUN cd build
RUN cmake -DCMAKE_EXE_LINKER_FLAGS="-llapack -lblas" -DCMAKE_Fortran_FLAGS="-std=legacy -ffree-line-length-none" ..
RUN make

FROM python:3.11-slim

COPY --from=build-env /src/build/CO2Transport /usr/local/bin/CO2Transport
COPY api api
WORKDIR /api

RUN pip install .

USER 1001

EXPOSE 8001

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

CMD ["uvicorn", "phPhitz_api.app:app", "--host", "0.0.0.0", "--port", "8008"]
