import yaml

docs = yaml.load_all(open('first.yaml'))

for doc in docs:
    for k,v in doc.items():
                print (k, "->", v)
