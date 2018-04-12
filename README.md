![CG2](https://raw.githubusercontent.com/labordynamicsinstitute/cg2/master/doc/images/cg2_logo.png)

# CG2 Fixed Effects Estimation Software

For more details, see [https://labordynamicsinstitute.github.io/cg2](http://labordynamicsinstitute.github.io/cg2) (latest).

## Overview
CG2 is a package of Fortran 90 and SAS programs that estimate non-nested 2 component fixed effect models.  The programs have been successfully compiled and tested on multiple Unix platforms (IA32 (SuSe and Redhat Linux), Itanium (SuSe Linux), and Sparc (Solaris)), but the code should work fine on other platforms as well (ie. IA32 (Windows)). A port to [Stata](http://www.stata.com/) (a2reg) is maintained by [Amine Ouzad](http://www.ouazad.com/) (an early version is available in the [Stata branch](https://github.com/labordynamicsinstitute/cg2/tree/stata), the latest version can be found at https://ideas.repec.org/c/boc/bocode/s456942.html).

The estimation algorithms were developed to solve large-scale fixed person and firm effect wage models.  In the typical scenario, a person's earnings and place of employment are observed over time, with mobility of persons across firms.  This mobility or "connectedness" of persons enables the estimation of the model, although for problems with millions of persons and firms, obtaining results require a substantial amount of computing power.  For users with large problems to solve, the main constraint will be acquiring a computing platform with sufficient  physical memory or RAM.  For example, the largest problem we have solved required 137GB of physical memory and took zzzz hours of CPU time on a Sun Fire 12k server.  However, this does not imply that CG2 is inefficient, since no known software can solve the same size problem using fewer system resources.  To  further reduce single system image memory requirements, a cluster aware version of CG2 is under development.

> The description here was first created several years ago. References to obsolete hardware (IA32, [Sun Fire 12k server](http://en.wikipedia.org/wiki/Sun_Fire_15K) ) are to be excused.

If you need additional information, please see the [References][#References] section below.

## Requirements
In order to use the CG2 package, certain basic software must be available on a current Unix platform (The software will likely port to Windows just fine, but we have not devoted time to testing this assertion.  If anyone succeeds in using this software on IA32 Windows, please let us know and we will post the information here).

* suitable Fortran 90/95 compiler.  The software is known to work with Intel Compilers on IA32 and Itanium and with the Sun compilers on Sparc/Solaris.  The [Intel compilers](http://www.intel.com/cd/software/products/asmo-na/eng/compilers/index.htm) are available for Linux and Windows and are free of charge for non-commercial users.  The [Sun compilers](http://www.sun.com/software/products/studio/index.xml) (obsolete?) are available for Linux and Sparc/Solaris, but are NOT available free of charge.  Other Fortran compilers will likely work as well, but are not supported.
* SAS statistical/data management software package from the [SAS Institute](http://www.sas.com).  SAS is available for a [wide range of hardware/software combinations](http://support.sas.com/documentation/installcenter/). It is used to pre- and post-process data.  If SAS is not available for your platform, other software packages MIGHT be used to pre and post-process the data, but this is not supported.
* A suitable shell environment such as Bash or Ksh, while not strictly required is very desirable.  Either shell should be available on virtually all Unix platforms.

Even if the above requirements can be met, the user must insure that their computing platform has sufficient physical memory. The [calculator](doc/html/calc.html) will allow such an estimation.

## Download
Download the distribution using the buttons on this page, or clone the repository.

## Installation

After download, there should be three directories:

* _runcg_
* _runcg_out_, and
* _cgpost_

## Using CG2

Before actually estimating model parameters, you must
properly structure the input file and calculate the size of the
problem. With the input data ready, the next step is to run CG2
and calculate parameter estimates. The final step is
post-processing the estimates, although this may not be required
depending on the set of parameters you are interested in.

### Test data
The SAS programs as packaged are setup to use a purely fictional,
synthetic dataset. It is a good idea to run through each HowTo on
the synthetic data first before attempting a problem of your own.
This will give you an opportunity to become familiar with each
stage of the process as well as test your installation.

### Typographical convention
We have tried to be consistent throughout in the use of typographical convention:

* _italics_ designate a _file_ or _directory_
* `type-writer` characters are used to designate
    * either `variables` within programs
    * or `commands` you type at the command prompt

## Howtos

* [Pre-process input data](doc/html/HowTo_preproc.html)
* [Run CG2](HowTo_runcg2.html)
* [Post-process the CG2 output](HowTo_postproc.html)


## Acknowledgements
John Abowd @johnmabowd, Robert Creecy, Kevin McKinney, and Lars Vilhuber @larsvilhuber.  Census Bureau and the rest of the LEHD staff.

## Contact
If you have any questions or comments please contact Kevin McKinney at kevinm@ccrdc.ucla.edu.

## References

* John M. Abowd & Robert H. Creecy & Francis Kramarz, 2002.
"Computing Person and Firm Effects Using Linked Longitudinal Employer-Employee Data,"
[Longitudinal Employer-Household Dynamics Technical Papers 2002-06](https://ideas.repec.org/p/cen/tpaper/2002-06.html), Center for Economic Studies, U.S. Census Bureau.
