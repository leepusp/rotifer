( hmmer2table -c model=version inputs/hmmscan.out; phobius2table -e 0.0101 inputs/phobius.out ) | domain2architecture -e 0.0101 | architecture2table | padtable
