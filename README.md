# CombiChem

<!--- These are examples. See https://shields.io for others or to customize this set of shields. You might want to include dependencies, project status and licence info here --->
![GitHub repo size](https://img.shields.io/github/repo-size/karanicolaslab/combichem)
![GitHub contributors](https://img.shields.io/github/contributors/karanicolaslab/combichem)
![GitHub stars](https://img.shields.io/github/stars/karanicolaslab/combichem?style=social)
![GitHub forks](https://img.shields.io/github/forks/karanicolaslab/combichem?style=social)

## Abstract

Our approach uses fragment screening techniques to identify the best fragments that are then concatenated into lead-like compounds. This allows us to significantly reduce the computational resources while performing a more comprehensive screen due to the smaller conformational space of fragments as compared to whole compounds

## Prerequisites
[RDKit](https://www.rdkit.org)
```
pip install rdkit 
OR
conda install -c rdkit rdkit 
```
[OEChem](https://www.eyesopen.com/oechem-tk) or [OpenBabel](http://openbabel.org/)
```
conda install -c openeye openeye-toolkits 
OR
conda install -c conda-forge openbabel 
```
[PyRosetta](http://www.pyrosetta.org). Please, download package from `<Latest PyRosetta Versions>` section
```
pip install <PyRosettaDistro>.tar.bz2
```
[Mol2Params](http://www.pyrosetta.org/scripts#TOC-DNA-Docking), `Small Molecule Docking` section
```
wget http://graylab.jhu.edu/pyrosetta/downloads/scripts/toolbox/molfile2params.tar.gz
tar xf molfile2params.tar.gz -C <combichem directory>
```
[Omega2](https://www.eyesopen.com/omega), if you prefer OpenEye tools to generate conformers

After installation you need to add `environment variables` for tools and scripts into your shell initializing file (e.g. `~/.bashrc`):
```
export OECHARGE=<path-to-assigncharges.py> # you can select only OECHARGE (paid commercial tool) or
export OBABEL=<path-to-obabel>             # OBABEL(free open source)

export OEOMEGA=<path-to-omega2>            # (not required if you prefer to use RDKit 
                                           # conformers generator)
```

## Navigating the Repository

- `Example` contains example input file to present how works each step of the approach
- `Utils` service script to support IO, search and structural manipulation with molecular structures
- `Alignment.py`, `Conformers.py`, `Merging.py`, `Params.py`, `Reactor.py`, `Screening.py` are separate steps of protocol
- `FragmentsMerging.py` and `FragmentsScreening.py` - multistep-protocols for fragments prioritization and merging to lead-like compounds

## Fragments library generation

If you don't have a fragments library yet, you can use our script for generating. It assepts buliding blocks and SMARTS reaction template as input and then returns fragments (i.e. functional group of building blocks merged with hinge-binding core) 

```
python Reactor.py -reactants <R1> ... <RN> -reaction <SMARTS_template> -out <SMILESOutput>
```

Alternatively, you can use [ChemAxon Reactor](https://chemaxon.com/products/reactor). This tools analagous to out script, but provides GUI to draw the transformation SMARTS template. Also, we can recommend to use [Dimorphite-DL](https://durrantlab.pitt.edu/dimorphite-dl/) to assign correct ionization state for fragments.

## Single steps of protocol

```
python Conformers.py -inp <SMILESInput> -out <SDFOutput>
python Params.py -ref_pdb <Reference_PDB> -ref_params <Reference_Param> -target_pdb <Target_PDB> -target_params <Target_Param> -out <Output_Param>
python Screening.py -pdb <Protein/LigandComplex> -params <LigandParamsFile> -output <OutputPDBFile> -residues <KLIFSseq> -beta True > <LogFile>
python Alignment.py -ref <Reference_SDF> -target <Target_SDF> -out <OUTPUT>
python Merging.py -ref_pose <Reference_SDF> -target_pose <Target_SDF#1> <Target_SDF#N> -substructure <HingeCore> -output <OUTPUT>
```

## Fragments screening

For simplicity of fragments virtual screening, we have implemented protocol containing separate manipulation with fragments structure. The script takes fragments in SMILES format, substructure (hinge-binding core in the paper), its `pdb` and `params` files, protein structure in `PDB` format, and 85-residue sequence from [KLIFS](https://klifs.net/) for given kinase and returns minimized conformers of fragments in complex with target protein.

```
python FragmentsScreening.py -smi examples/FragmentsScreening/R1.smi \
                             -substructure examples/FragmentsScreening/ref.sdf \
                             -ref_pdb examples/FragmentsScreening/ref.pdb \
                             -ref_param examples/FragmentsScreening/ref.params \
                             -protein examples/FragmentsScreening/CDK9.pdb \
                             -klifs_seq 23 24 25 26 27 28 29 30 31 32 33 34 35 45 46 47 48 49 50 62 63 64 65 66 67 68 69 70 71 72 73 74 76 77 78 79 80 81 82 83 84 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 165 166 167 168 169 170 171
```

## Fragments merging
Also we created protocol for merging of fragments and screening new lead-like compounds. The script accepts minimized protein-fragment complexes for each R-group (e.g. R1,R2,R3), their params files, substructure (to define overlap between fragments) and 85-residue sequence from [KLIFS](https://klifs.net/) for given kinase structure. Please, note that we selected one of the fragments (R1 in this case) as reference to realign other fragments and to use protein structure for minimization.

```
python ../../FragmentsMerging.py -ref_pdb examples/FragmentsMerging/R1.pdb \
                                 -ref_params examples/FragmentsMerging/R1.params \
                                 -target_pdb examples/FragmentsMerging/R2.pdb \
                                 -target_params examples/FragmentsMerging/R2.params \
                                 -substructure examples/FragmentsMerging/Core.sdf \
                                 -klifs_seq 23 24 25 26 27 28 29 30 31 32 33 34 35 45 46 47 48 49 50 62 63 64 65 66 67 68 69 70 71 72 73 74 76 77 78 79 80 81 82 83 84 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 165 166 167 168 169 170 171
```

## Contact
The corresponding author for this work is John Karanicolas who can be reached at john.karanicolas@fccc.edu.
## License
This project uses the following license: [MIT License](https://github.com/karanicolaslab/combichem/blob/main/LICENSE).

