# BsHcalling
Explorating methods for best hit calling generated from homology-based scanning tools.


* [Introduction](#introduction)
* [Version](#version)
* [Requirement](#requirement)


## Introduction

Calling best hits for blast output is a diffcult task due the fact that abundant segemental alignments would be generated using a large scale of database. Though the stats bitscores and evalues serve as good indexes for filtering low quality hits, manual inspection is still needed to determine final best hits. To call the best hit for each query contig from the raw blast result generated from a large data processing, we introduced a statistic named 'adjusted identity of merge hits' (AIMH). Given a query contig Q and a subject sequence S (reference matched by query), equal to or more than one hit (H) was observed. The AIMH was calculated by merging HSPs and calculated a combined identiy weighted by the segmental HSPs.


## Version
+ AIMH 0.1 

## Requirement
------------
+ Python 3 (tested on 3.8.3 / 3.6.9 )
+ pandas 1.3.4 


