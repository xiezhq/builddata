The package include interf and some other utilities (asa2pep, asa2prot and searchpdb). The interf calculates the interface between subunits in a protein complex structure. The codes were basically written in 2005-2006 but a few of them were written almost 20 years ago (1997). I currently don't have a plan to improve these obsolete codes.

I didn't write the detailed help function for these programs as they were initially developed for my own use. If you have some experiences in C, you can have a look at the C codes and can easily figure out what the programs do and how to run them.

The parameter files of atomic VDW (Van Der Waals) radii are placed in par directory. The files in peptide directory are atomic coordinate file (PDB format) of the standard tripeptide structures I built using InsightII 2005 package, and the dipeptides for N- and C-terminal residues are also included. 

Zhiqun Xie - Mar 9, 2017

## How to install and run

Under Linux/Unix, you can simply compile all programs in the package.
```sh
make
```
If you just like to compile the specific program, e.g. interf, you should:
```sh
make interf
```
Once the compiler sucha as gcc successfully compiled the codes, you can run the program now.
For example, you want to run interf:
```sh
./interf
```
