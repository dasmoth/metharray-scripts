Scripts for normalization and analysis of Illumina methylation array data.

Requires recent Java and ANT.  To compile, type:

       ant

Then add metharray.jar, bjv2-core-0.1.jar and colt.jar to your CLASSPATH

These scripts expect data as tab-delimited matrices with a probe ID in the first
column and numerical values in all other columns.  You can add some header lines
(starting with '#') if you like.


Normalization
-------------

       java data.t1d.MatrixQuantileApplication matrix.txt >matrix.qn.txt

Calling differences
-------------------

Samples are divided into "foreground" (treatment, mutant, affected, whatever)
and "background" (controls, normals).  They are identified by *zero-based*
column indices.

       java data.t1d.PairwiseDMRCallerWApplication -fg 1,3,5 -bg 2,4,6 matrix.qn.txt >dmrcalls.txt

P-values are in column 3 of the output file.