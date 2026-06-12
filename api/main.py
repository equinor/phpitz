import os
from fastapi import FastAPI
from pydantic import BaseModel, Field

import subprocess
from contextlib import chdir

app = FastAPI()


RUNFOLDER = os.getenv("RUN_FOLDER", ".")
TEMPLATE = """ 5
 0
 0
 1
 1
 0
 1
 1
 12
 1
 1
 1
 1
 1
 1
 1
 1
 0
 0
 0
 0
 1
 0
 1
 {pressure}           {temperature}
 {co2}
 {h2o}
 {n2}
 {o2}
 {h2s}
 {no}
 {no2}
 {so2}
 {h2so4}
 {hno3}
 {s8}
 {nh3}
 {n2o}
 {n2o4}
 {nh4hso4}
 {hcho}
 {ch3cho}
 {ch3coch3}
 {hcooh}
 {ch3cooh}
"""


@app.get("/")
async def read_root():
    return {"message": "running server"}


class Concentrations(BaseModel):
    h2o: float = Field(default=0)
    o2: float = Field(default=0)
    so2: float = Field(default=0)
    no2: float = Field(default=0)
    h2s: float = Field(default=0)
    no: float = Field(default=0)
    h2so4: float = Field(default=0)
    hno3: float = Field(default=0)
    s8: float = Field(default=0)
    n2: float = Field(default=0)
    nh3: float = Field(default=0)
    n2o: float = Field(default=0)
    n2o4: float = Field(default=0)
    nh4hso4: float = Field(default=0)
    hcho: float = Field(default=0)
    ch3cho: float = Field(default=0)
    ch3coch3: float = Field(default=0)
    hcooh: float = Field(default=0)
    ch3cooh: float = Field(default=0)


class SimulationRequest(BaseModel):
    temperature: int
    pressure: int
    concentrations: Concentrations


@app.post("/api/run")
async def calculate(simulationRequest: SimulationRequest):
    result = simulationRequest.model_dump()
    concentrations = result["concentrations"]
    concentrations["co2"] = 1000000 - sum(concentrations.values())

    with chdir(RUNFOLDER):
        with open("../datafiles/CalcInput.dat", "w") as fh:
            fh.write(
                TEMPLATE.format(
                    **concentrations,
                    pressure=result["pressure"],
                    temperature=result["temperature"],
                )
            )

        result = subprocess.run("CO2Transport", text=True, stdout=subprocess.PIPE)

    result = result.stdout.strip()
    return result
