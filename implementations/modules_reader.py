import sys
sys.path.insert(0, "../")

def read_modules(file_path):
    modules = []

    with open(file_path, 'r') as f:
        line = f.readline()
        while line != "":
            line = f.readline()
            if line.startswith("cc"):
                modules.append(f.readline().strip()[1:-1].split(', '))

        return modules
