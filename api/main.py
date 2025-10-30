import os
from fastapi import FastAPI
from pydantic import BaseModel, Field

import subprocess
from contextlib import chdir

app = FastAPI()


RUNFOLDER = os.getenv("RUN_FOLDER", ".")
TEMPLATE = """
 {pressure}           {temperature} 
CO2
 {co2}        0 
H2O
 {h2o}             0 
N2
 {n2}            0 
O2
 {o2}             0 
H2S
 {h2s}            0 
NO
 {no}             0 
NO2
 {no2}             0 
SO2
 {so2}             0 
H2SO4
 {h2so4}             0 
HNO3
 {hno3}             0 
S8
 {s8}             0 
 5 
 1 
 1 
 1 
 1 
 1 
 5 
 1 
 1 
 1 
 1 
 0 
 1 
 {pressure}           {temperature} 
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
    hno2: float = Field(default=0)
    n2: float = Field(default=0)


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
        with open("../../datafiles/CalcInput.dat", "w") as fh:
            fh.write(
                TEMPLATE.format(
                    **concentrations,
                    pressure=result["pressure"],
                    temperature=result["temperature"],
                )
            )

        result = subprocess.run(
            "./CO2phaseReactions", text=True, stdout=subprocess.PIPE
        )

    result = result.stdout.strip()
    return result
