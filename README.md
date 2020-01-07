# XLM-Tools

An integrated tool to score model PDB structures based on crosslink and monolink data.

Scoring for Crosslinks is performed using MNXL and Monolinks are scored by MoDS. If both crosslinks and monolinks are being scored the combination score XLMO is calculated.

Information on these scores can be found in [publication citation here].

## Getting Started

Run the following command to install the XLMTools and Jwalk modules:

```
python setup.py install
```

To run please copy the [xlmtools.v1.0.py](xlmtools.v1.0.py) script into the folder containing either pdb structures to generate Depth and/or Jwalk result files or the Depth/Jwalk files themselves. 

Flags can be listed by running python xlmtools.v1.0.py -h, which will show the below help message:

XLM Tools tools: A tool to score model protein structures according to
crosslink and monolink data.

optional arguments:
  -h, --help            show this help message and exit

XLM Arguments:
  -xl_list XL_LIST      input list of crosslinks and monolinks

  -jwalk_files [JWALK_FILES [JWALK_FILES ...]]
                        jwalk files for MNXL scoring

  -depth_files DEPTH_FILES [DEPTH_FILES ...]
                        depth files for MoDS scoring

  -outfile_name OUTFILE_NAME
                        specify output file name

  -sep SEP              separator in output file, default is tab

  -pdb PDB [PDB ...]    specify pdb files for Jwalk/Depth run

Jwalk Arguments:

  -jwalk                flag to use if starting from .pdb files and running
                        Jwalk

  -vox VOX              specify voxel size of grid

  -surface              use higher accuracy method to calculate solvent
                        accessibility - requires Freesasa installation

  -ncpus NCPUS          specify number of cpus to use

Depth Arguments:

  -depth                flag to use if starting from .pdb files and running
                        Jwalk

  -depth_source DEPTH_SOURCE
                        specify depth source


The only required argument is -xl_list, this contains the crosslinks and monolinks which the model structures are scored against. The file is structured as below:

aa1|c1|aa2|c2 <- Crosslink

aa1|c1|aa2|c2

aa1|c1|aa2|c2

aa1|c1|aa2|c2

aa1|c1        <- Monolink

aa1|c1

aa1|c1

### Prerequisites

Python 3.6

Scipy

Numpy

Biopython

DEPTH - Install from http://cospi.iiserpune.ac.in/depth/htdocs/download.html. Either set up an alias named "DEPTH" to the program or specify the path to your source using the -depth_source flag. Note: may not work on all OS, see website for details.

## Authors

* **Matthew James Sinnott**
* **Dr Joshua Matthew Allen Bullock**
* **Prof. Konstantinos Thalassinos**
* **Prof. Maya Topf**

## License

This project is licensed under the GNU GPL v3.0 License - see the [LICENSE](LICENSE) file for details

## Citing

If this software is used please cite:
```
citation
```

## Acknowledgments

* **Dr David Holdershaw**
* **Dr Sony Malhotra**
* **Dr Mark Williams**