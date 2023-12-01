# BsHcalling
Explorating methods for best hit calling generated from homology-based scanning tools.

* [Introduction](#introduction)
* [Version](#version)
* [Requirement](#requirement)
* [Installation](installation)
* [Usage](usage)

## Introduction

A collection of utils for calling best hit from results generated from homology-based scanning tools. 
In this version, we only included the AIMH (Adjusted Identity of Merged Hits), a statics designed for calling best hit from blast result (blastn, blastx, blastp)

## Version
+ 0.1 

## Requirement
+ Python 3 ( tested in 3.8.3)
+ pandas >= 1.3.4

## Installation
Just git clone this project and pip install
```
git clone https://github.com/ZehanDai/BsHcalling.git
pip install .
```

## Usage
To print help information:
```
python3 -m AIMH -h
```

To use it on input of format 6 blast output:
```
python3 -m AIMM \
     -i path_to_input_table_file\
     -o path_to_output_file1\
     -O path_to_output_file2.fmt6\
     -s input_sep\
     -S output_sep\
     -l 'original' 
```

