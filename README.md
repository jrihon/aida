# AIDA 

This repository hosts the software employed in the following [publication](https://academic.oup.com/nar/article/35/19/6611/2402188) : 
 
*Backbone-base inclination as a fundamental determinant of nucleic acid self- and cross-pairing* , by Pallan et al.

## DISCLAIMER !

A problem was discovered where several output parameters hold an incorrect sign.
This problem will be looked into in the foreseeable future.

## Dependencies

```shell
# `csh` old shell dialect, similar to bash
# `fpc` free pascal compiler
#
$ sudo apt install fpc csh
```
### Versions tested 
- `fpc 3.2.2`
- `fpc 3.0.4`
- `csh 20110502-7`

## Software
This repository has been successfully installed on `Ubuntu20`, `Pop! OS 22` (which used Ubuntu22 at its core), `Ubuntu22` and `Ubuntu24`. Only tested in `Linux`.

### Shell scripts

- run/run_inclination.sh (run file for calculations)
- bin/scherung.sh (called indirectly through run_inclination.sh)
- bin/analyse.sh (called indirectly through run_inclination.sh)

### Awk script
- bin/prepPDB.awk
- bin/getC.awk
- bin/getBP.awk
- bin/getB.awk
- bin/getPP.awk

### Pascal source files
Modified pascal files to allow compilation using the `fpc`
- src/BSPabstand.p
- src/best_plane.p
- src/dehedral.p
- src/BSScherung.p
- src/PPScherung.p
- src/evaluate.p
		
## Compilation
```shell
$ cd src/
$ ./compile.sh
```

## Run calculation

### Preparation

Before you execute the program, add two files to `run/`.

- file with the name of the pdb-file, without its extension ( ! `.pdb`).
- file with the basepairing information, with the extension (`.bp`), ignoring the terminal ends.

### e.g.
- pdb-file : `adh027s`
- file with basepairing information : `adh027s.bp`

Content of with `adh027s.bp` :
```
G 2 C 15
G 3 C 14
C 4 G 13
G 5 C 12
C 6 G 11
C 7 G 10
```

###	Execute 
 
- $1 : pdbfile (no file-extension)
- $2 : Amount of atoms in the leading strand you want to include (number)
- $3 : Atom you want to target (C1', P, etc.)
```shell
$ cd run/
# script                $1      $2 $3
$ ./run_inclination.sh  adh027s 7  P
```

## Disclaimer !

I am not the original author of the code. The code belongs to the authors of the following publication : 
*Backbone-base inclination as a fundamental determinant of nucleic acid self- and cross-pairing*

If you use this software, or any of its parts, please cite the following publication : 
```bib
@article{Pallan_2007,
    title={Backbone-base inclination as a fundamental determinant of nucleic acid self- and cross-pairing},
    volume={35},
    ISSN={0305-1048},
    url={http://dx.doi.org/10.1093/nar/gkm612},
    DOI={10.1093/nar/gkm612},
    number={19},
    journal={Nucleic Acids Research},
    publisher={Oxford University Press (OUP)},
    author={Pallan, Pradeep S.
        and Lubini, Paolo
        and Bolli, Martin
        and Egli, Martin},
    year={2007},
    month=sep,
    pages={6611–6624}
}
```

### Authorisation
The author of this repository has received formal approval from `prof. dr. M. Egli`, the corresponding author, to host the updated version of the code on GitHub!

## Changes to the original code
1. Installation
- Repackaged directory
- Improve original README
- Use the free pascal compiler, available in common repositories `fpc`
- Figure out to use the Objective Pascal `-Mobjfpc` flag to compile correctly

2. Software
- Introduce `string` primitive type instead of self-declared `string` type (deprecated)
- Introduce `FileHandler` functionality with the `AssignFile` procedure
- Close any open-ended comments (`pascal` comments are of the type `{...}`)
- Solve type problems (`string`, `TextFile`)
- Correctly format pascal code through `ptop`
- Removed unnecessary `escape sequence` characters in `awk` files
- Remove functions being passed as arguments and explicitly call the function inside a procedure


### Future updates (maybe)
- Translated German to English
- Convert `csh` to the more contemporary `bash`



## Known issues
If you encounter problems of the following type : `runtime error 106` while running the `run_inclination.sh` script, this can be resolved through changing `locale` settings on your linux distribution. More information on runtime errors can be found [here](https://www.freepascal.org/docs-html/user/userap4.html).
```shell
# compile the pascal program with adjusted locale settings
$ cd src/
$ LC_NUMERIC="en_US.UTF-8"
$ ./compile.sh 

# run the command 
$ cd ../run/
$ ./run_inclination.sh adh027s 7 P
```
This change to the `LC_NUMERIC` variable is local and resets every time your close the environment.


## License
This repository is made public under the GPL v3 License. Additionally, this software is shared solely for research purposes and and cannot be used to develop a product or service that is meant to be sold commercially
