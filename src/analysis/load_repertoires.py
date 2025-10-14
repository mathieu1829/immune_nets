
from src.models import Repertoire,Dataset
from src.db import engine
from sqlalchemy.orm import Session
from pathlib import Path
from src.factories import RepertoireFactory


root_dir = Path(__file__).parent.parent.parent

leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

with Session(engine) as session: 
    healthy = Dataset(name="healthy", description=" ")
    leukemia = Dataset(name="leukemia", description=" ")
    covid = Dataset(name="covid", description=" ")

    healthy.repertoires = [
            RepertoireFactory.fromCSV(name="healthy",desc=" ",path=healthy_path)
            ]
    covid.repertoires = [
            RepertoireFactory.fromCSV(name="covid",desc=" ",path=covid_path)
            ]
    leukemia.repertoires = [
            RepertoireFactory.fromCSV(name="leukemia",desc=" ",path=leukemia_path)
            ]
    session.add_all([healthy,leukemia,covid])
    session.commit()
