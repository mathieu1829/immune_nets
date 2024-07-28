from pathlib import Path

class pathManager:
    def __init__(self):
        self.testDataPath = Path(__file__).parent.parent.parent.parent / "tests/test_data"

if __name__ == "__main__":
    a = pathManager()
    print(a.testDataPath)
