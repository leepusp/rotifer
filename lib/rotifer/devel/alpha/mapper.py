# coding: utf-8
from collections import OrderedDict

# Default behaviour
class configuration(dict):
    def __init__(self, **kwargs):
        self.include = OrderedDict({
            "annotation":[],
            "dbxref":[],
            "feature":[],
            "qualifier":[],
            "sequence":False,
        })
        self.exclude = OrderedDict({
            "annotation":["comment", "structured_comment", "sequence", "references"],
            "dbxref":[],
            "feature":[],
            "qualifier":[],
            "sequence":True,
        })

    def __getitem__(self,path):
        a = self.__dict__
        for x in path.split("."):
            if x in a:
                a = a[x]
        return a


    def to_yaml(handle):
        pass

# Methods
def seqrecords_to_dataframe(seqrecords, mapping=SeqRecordToDataFrameMap()):
    """
    Generic mapping of Bio.SeqRecord objects to a Pandas DataFrame.
    """
    import pandas as pd
    annotation_prefix='Annotation'
    feature_prefix='Qualifier'
    dbxref_prefix='Dbxref'

    # Loop over seqrecords
    d = []
    fd = []
    dicts = []
    for s in seqrecords:
        # Prepare base dictionary
        data = {
            "id": s.id,
            "sequence": None,
            "type": "sequence",
            "length": len(s),
            "description": s.description,
        }
        if (mapping.include["sequence"] or not mapping.exclude["sequence"]):
            data["sequence"] = str(s.seq)

        # Annotations
        for k, v in s.annotations.items():
            k = k.replace("-", "_")
            if (not include[k] or exclude[k]):
                continue
            key = f'{annotation_prefix}__{k}'
            if (include and (key not in include)) or (key in exclude):
                continue
            if k == "comment":
                data[k] = v.replace("\n", '\n').strip()
            elif k == "structured_comment":
                for k1, v1 in v.items():
                    k1 = k1.replace("-", "_")
                    k1 = f'{annotation_prefix}__{k}__{k1}'
                    data[k1] = v1.strip()
                continue
            elif k == "taxonomy":
                data[k] = "; ".join(v).strip()
            elif k == "db_source":
                data[k] = v.replace("accession ", "")
            elif isinstance(v, list):
                if len(v) == 1:
                    data[k] = v[0].strip()
                else:
                    data[k] = [ x.strip() for x in v ]
            elif isinstance(v, dict):
                data[k] = v
            else:
                data[k] = v.strip()
            data[key] = data[k]
            del data[k]

        # Database references
        for y in s.dbxrefs:
            y = y.split(":")
            data[f'{dbxref_prefix}__{y[0]}'] = y[1]

        # Features
        for f in s.features:
            # Standard feature data
            fdata = {
                "id": s.id,
                "type": f.type,
                "start": f.location.start + 1,
                "end": f.location.end,
                "strand": f.location.strand,
                "length": len(f),
            }

            # Process feature qualifiers
            for k, v in f.qualifiers.items():
                if isinstance(v, list):
                    if len(v) == 1:
                        v = v[0]
                k = f'{feature_prefix}__{f.type}__{k}'
                fdata[k] = v

            # Store qualifiers of features spanning the entire
            # sequence length as annotations
            if (fdata["start"] == 1
                and fdata["end"] == len(s)
                and (fdata["strand"] == None or fdata["strand"] == 1)
            ):
                for k in fdata.keys():
                    if k in ['id','type','start','end','strand','length']:
                        continue
                    data[k] = fdata[k]
            else:
                fd.append(fdata)

        # Append to stack
        d.append(data)

    # Build dataframes
    fd = pd.DataFrame(fd)
    d = pd.DataFrame(d)

    # Return dataframe
    return (d, fd)
