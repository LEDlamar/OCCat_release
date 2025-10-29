# README

## Prepartation

The structure of the repository at https://github.com/LEDlamar/OCCat_release is:

```markdown
# Root directory
.
├── bin                        
│   └── occat                  # Binary executable
├── desc_nanoalloy             # Python package of catalytic performance analyzer
│   ├── calcRate.py            # Introduce more descriptors
│   ├── orr.py                 # Demo of analyze a specific reaction
│   └── ...
└── example                    
    ├── Ag40O4_Lasp            
    │   ├── Ag40.xyz           # XYZ coordinate file of base cluster
    │   ├── AgO.pot            # NN potential file for Ag-O system
    │   ├── config.toml        # Configuration file of OCCat
    │   ├── lasp.in            # LASP input file
    │   └── submit.sh          # Shell script for job submission
    ├── Ag8_Vasp               
    │   ├── INCAR_1            # INCAR for VASP
    │   ├── config.toml        # Configuration file of OCCat
    │   └── submit.sh          # Shell script for job submission
    ├── Co62_Gupta             
    │   └── config.toml        # Configuration file of OCCat
    └── LJ21                   
        └── config.toml        # Configuration file of OCCat
```

## Quick start

### 1. Requirements

+ Golang >= 1.18
+ Python 3 with Numpy >= 1.20 (When using catalytic performance analyzer)
+ Unix-like environment (Windows is not supported)
+ If using VASP/LASP, the corrsponding environment is requested. See [VASP - Vienna Ab initio Simulation Package](https://www.vasp.at/) and [LASP Software](http://www.lasphub.com/#/lasp/laspHome)

### 2. How to run

+ **Prepare the configuration file**:

  - The configuration file can be named arbitrarily (e.g., `config.toml`, `my_config.toml`, etc.).
  - Ensure the file is correctly formatted and parsed.

+ **Run the global optimization**:
   Execute the following command to start the optimization process:

  ```sh
  occat --config <config.toml>  
  ```



