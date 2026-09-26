
from io import TextIOWrapper

class TextIO(TextIOWrapper):
    def __init__(self, fileName, mode='r', encoding='utf-8', delete=False):
        super().__init__()
        self.__delete = delete

    def __enter__(self):
        return self.fh

    def __exit__(self, exc_type, exc_value, traceback):
        self.__close__()

    def __del__(self):
        if self.__delete:
            os.remove(self.name)

    def close(self):
        self.fh.close()
        self.__delete__()


