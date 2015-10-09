# RESP tools for multiconformer fitting using GAMESS and RESP

```
John Chodera 
2007-03-12
Department of Chemistry, Stanford University

DIRECTORIES
python/     Python code
resp/       Source code for RESP, redimensioned to handle many more conformers
examples/   Example of use on Guthrie test set

TOOLS
python/setupBatchResp.py     Set up batch queue GAMESS calculations for RESP given a multi-molecule mol2 file
python/analyzeBatchResp.py   Analyze results of batch GAMESS calculations

PYTHON CLASSES
python/BatchResp.py    RESP and GAMESS interaction methods used by the tools above
python/Resp.py         Incomplete Resp class -- ignore this for now

TESTS
python/testResp.py     Test harness for python/Resp.py -- ignore this for now
```
