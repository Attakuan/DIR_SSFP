
# Introduction
Unofficial Double Inversion Recovery (DIR) for b-SSFP and nb-SSFP. The code is from the original repository: [*Pulseq*](https://github.com/pulseq/pulseq). 


*Pulseq* is an open source framework for the development, 
representation and execution of magnetic resonance (MR) sequences. A central contribution 
of this project is an **open file format** to compactly describe MR sequences 
suitable for execution on an MRI scanner or NMR spectrometer. 
MATLAB and C++ source code is provided for reading and writing sequence files. The main homepage for Pulseq is
http://pulseq.github.io/. 

The directories are organized as follows:

* `doc/` - Contains the file specification and HTML source code documentation
* `examples/` - Contains example sequence files (`*.seq`)
* `src/` - C++ class for reading sequence files
* `matlab/` - MATLAB code for reading, writing, modifying and visualizing sequence files

 ## How to Use

The main file of interest for applying DIR SSFP is in MATLAB directory with the name `DI_bSSFP`. You basically have to run the code with the desired parameters. The code will generate a sequence file with the name `xx.seq` in the same directory. Then, you can use this sequence file to run the sequence on your scanner via USB or Ethernet cable. 

Note that the code is written for 3T Siemens scanners. Thus, you have to modify the code for other scanners. Furthermore, DIR is implemented using adiabatic inversion pulses. In the original source code, there is a function for generating adiabatic pulses. However, it is not compatible with Windows. Thus, I have modified the source code to be able to generate adiabatic pulses using cmd and anaconda. You have to check the main code for the path of the adiabatic pulse which is in `matlab/+mr/makeAdiabaticPulse`. Modify this code according to your path and it should work fine. 


## System requirements

I have modified *Pulseq* source code to be able to generate adiabatic pulse using cmd and anaconda. Thus, this version of *Pulseq* is compatible only with Windows for DIR SSFP but can be used for other sequences on Linux and MacOS. Otherwise, main requirements is the same as [*Pulseq*](https://github.com/pulseq/pulseq).

## Example Results
![Transformation Preview](DIR_SSFP_Example.jpg)

