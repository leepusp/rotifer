#!/usr/bin/env python3

from sqlalchemy import Column, String, Integer, Date, Float, MetaData

from sqlalchemy.dialects import registry
registry.register('clickhouse', 'base', 'dialect')
# import rotifer.alchemy.sqlbase
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

metadata = MetaData(Base)
class pfam():
    __tablename__ = 'Pfam_32_0'


    sequence = Column(String)
    domain = Column(String)

    start = Column(Integer)
    end = Column(Integer)

    evalue = Column(Float)
    qstart = Column(Float)
    qend = Column(Integer)
    qcov = Column(Integer)
    iteration = Column(Integer)
    score = Column(Float)
