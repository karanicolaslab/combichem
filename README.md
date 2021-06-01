# CombiChem

<!--- These are examples. See https://shields.io for others or to customize this set of shields. You might want to include dependencies, project status and licence info here --->
![GitHub repo size](https://img.shields.io/github/repo-size/karanicolaslab/combichem)
![GitHub contributors](https://img.shields.io/github/contributors/karanicolaslab/combichem)
![GitHub stars](https://img.shields.io/github/stars/karanicolaslab/combichem?style=social)
![GitHub forks](https://img.shields.io/github/forks/karanicolaslab/combichem?style=social)

## Abstract

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

## Navigating the Repository

## Input files for fragments screening

## Input files for fragments merging

## Example commands



```
python Reactor.py -reactants <R1> <RN> -reaction <SMARTS_template> -out <SMILESOutput>
python Conformers.py -inp <SMILESInput> -out <SMILESOutput>
python Params.py -ref_pdb <Reference_PDB> -ref_params <Reference_Param> -target_pdb <Target_PDB> -target_params <Target_Param> -out <Output_Param>
python Screening.py -pdb <Protein/LigandComplex> -params <LigandParamsFile> -output <OutputPDBFile> -residues <KLIFSseq> -beta True > <LogFile>
python Alignment.py -ref <Reference_SDF> -target <Target_SDF> -out <OUTPUT>
python Merging.py -ref_pose <Reference_SDF> -target_pose <Target_SDF#1> <Target_SDF#N> -substructure <HingeCore> -output <OUTPUT>

```

## Contact
The corresponding author for this work is John Karanicolas who can be reached at john.karanicolas@fccc.edu.
## License
This project uses the following license: [MIT License](https://github.com/karanicolaslab/combichem/blob/main/LICENSE).

