SSRL pyMCA viewer
=================


Running
-----
You can run the tool with the following: `pydm -m DEVICE=XPRESS_DEMO SSRL_MCA.py`
Be sure to change the name of the device to the proper device name. 

Alternatively, you can run with various files in the scripts folder such as:
```
./scripts/run_me.sh
```

Different `.sh` files run different versions of the UI, which could involve different fitting protocols, etc.  
You may have to `chmod +x *.sh` the files in the scripts folder to modify permissions. 

Note: you may run the tool in a debugging mode using the DEBUG flag via the command line. For example:
`pydm --log_level DEBUG SSRL_MCA.py`
or
`pydm -m DEVICE=XPRESS_DEMO --log_level DEBUG SSRL_MCA.py`
These debug commands could be added to the running scripts.

Usage
-----
The tool has two operational modes: live and static data processing. These different modes can be selected within the tool on the tab-based 
panel on the right side of the UI. The two tabs can be swapped to change the type of data being processed. 

On the live data panel, there are options to start and stop displaying data, and the exposure level and count.

If the user attempts to switch to a static file while live data is being received, a confirmation window will pop up to ensure that the user wishes to close the connection. If the user confirms the switch, but later
chooses to switch back to live data, a connection will attempt to be recreated with the previous data channel. 

Static Data Files
----------
Included in the repository are two sample data files to use within the tool. 
They can be interacted with via the file tab on the right side of the UI. 


Demo IOC
--------

This demo IOC uses pcaspy to simulate an EPICS IOC.

Open a terminal and start the IOC with the following command:
`python demo/testing-ioc`

In another terminal, run the following command to start the tool itself:
`pydm -m DEVICE=XPRESS_DEMO SSRL_MCA.py`

Debugging 
--------
When running the tool, there are a few common issues that may occur. The following points may be useful in debugging the environment.

1. Modify firewall access permissions on the local device to allow connections.

2. Set the max array bytes to accomodate larger values. 
    For example, on mac run: `export EPICS_CA_MAX_ARRAY_BYTES=100000000`