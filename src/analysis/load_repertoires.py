
from src.orm.models import Repertoire,Dataset
from src.orm.db import engine
from sqlalchemy.orm import Session
from pathlib import Path


root_dir = Path(__file__).parent.parent.parent

leukemia_path = root_dir  / "tests/test_data/leukemia_test_clonotypes.csv" # leukemia
covid_path = root_dir / "tests/test_data/covid_test_clonotypes.csv" # covid
healthy_path = root_dir / "tests/test_data/healthy_test_clonotypes_1.csv" #healthy

with Session(engine) as session: 
    healthy = Dataset(name="healthy", description=" ")
    leukemia = Dataset(name="leukemia", description=" ")
    covid = Dataset(name="covid", description=" ")

    healthy.repertoires = [
            Repertoire.from_csv(name="healthy",desc=" ",path=healthy_path)
            ]
    covid.repertoires = [
            Repertoire.from_csv(name="covid",desc=" ",path=covid_path)
            ]
    leukemia.repertoires = [
            Repertoire.from_csv(name="leukemia",desc=" ",path=leukemia_path)
            ]
    session.add_all([healthy,leukemia,covid])
    session.commit()
