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
 {reactive}
 {solubility}
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
 1
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
    temperature: float
    pressure: float
    concentrations: Concentrations
    reactive: bool = Field(default=True)
    solubility: bool = Field(default=True)


def parse_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        if "-" in value[1:] and "e" not in value.lower():
            idx = value.index("-", 1)
            return float(value[:idx] + "e" + value[idx:])
        raise


def parse_output(text: str) -> dict:
    result = {
        "raw": text,
        "phase_reactions": {},
        "aqueous_phase": None,
        "solubility": None,
    }

    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if "MOLE RATES AFTER" in line:
            result["solubility"] = {}
            dash_count = 0
            i += 1
            while i < len(lines) and dash_count < 2:
                if "---" in lines[i]:
                    dash_count += 1
                i += 1
            while i < len(lines) and "---" not in lines[i]:
                parts = lines[i].split()
                if len(parts) >= 4:
                    name = parts[0]
                    total_mol = parse_float(parts[1])
                    gas_mol = parse_float(parts[2])
                    water_mol = parse_float(parts[3])
                    pct_in_aqueous = parse_float(parts[4]) if len(parts) >= 5 else 0.0
                    result["solubility"][name] = {
                        "total_mol": total_mol,
                        "gas_mol": gas_mol,
                        "water_mol": water_mol,
                        "pct_in_aqueous": pct_in_aqueous,
                    }
                i += 1

        elif "No aquous phase present" in line:
            result["solubility"] = {}

        elif "CO2 PHASE REACTIONS" in line:
            i += 1
            while i < len(lines) and "---" not in lines[i]:
                i += 1
            i += 1
            while i < len(lines) and "===" not in lines[i]:
                parts = lines[i].split()
                if len(parts) >= 4:
                    result["phase_reactions"][parts[0]] = {
                        "feed_mol": parse_float(parts[1]),
                        "equil_mol": parse_float(parts[2]),
                        "equil_ppm": parse_float(parts[3]),
                    }
                i += 1

        elif "AQUEOUS PHASE RESULTS" in line or "ACID PHASE RESULTS" in line:
            result["aqueous_phase"] = {}
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if "Water pH" in line:
                    result["aqueous_phase"]["ph"] = parse_float(
                        line.split(":")[1].strip()
                    )
                elif "Water phase" in line:
                    result["aqueous_phase"]["water_mass_kg"] = parse_float(
                        line.split(":")[1].replace("kg", "").strip()
                    )
                elif "Water activity" in line:
                    result["aqueous_phase"]["water_activity"] = parse_float(
                        line.split(":")[1].strip()
                    )
                elif "H2SO4 in water" in line:
                    result["aqueous_phase"]["h2so4_wt_pct"] = parse_float(
                        line.split(":")[1].replace("wt%", "").strip()
                    )
                elif "HNO3  in water" in line:
                    result["aqueous_phase"]["hno3_wt_pct"] = parse_float(
                        line.split(":")[1].replace("wt%", "").strip()
                    )
                elif "NH3   in water" in line:
                    result["aqueous_phase"]["nh3_wt_pct"] = parse_float(
                        line.split(":")[1].replace("wt%", "").strip()
                    )
                elif "CO2   in water" in line:
                    result["aqueous_phase"]["co2_wt_pct"] = parse_float(
                        line.split(":")[1].replace("wt%", "").strip()
                    )
                elif "---" in line and result["aqueous_phase"]:
                    break
                i += 1

        i += 1

    return result


@app.post("/api/run")
async def calculate(simulationRequest: SimulationRequest):
    data = simulationRequest.model_dump()
    concentrations = data["concentrations"]
    concentrations["co2"] = 1000000 - sum(concentrations.values())

    with chdir(RUNFOLDER):
        with open("../datafiles/CalcInput.dat", "w") as fh:
            fh.write(
                TEMPLATE.format(
                    **concentrations,
                    pressure=data["pressure"],
                    temperature=data["temperature"],
                    reactive=int(data["reactive"]),
                    solubility=int(data["solubility"]),
                )
            )

        proc = subprocess.run("CO2Transport", text=True, stdout=subprocess.PIPE)

    return parse_output(proc.stdout)
