# Boron_red
## LA-MC-ICPMS boron isotope reduction script

**Disclaimer** : It only work with *Neptune plus* files and *Geostar* laser log files. This script will look for very specific columns and details of the files produced by these programs.

If needs be, I can make the code less rigid to fit other machines. Just contact me.

Methodology originally outlined in:

*Evans, Gerdes, Coenen, Marschall, & Mueller (2021) Accurate correction for the matrix 
interference on laser ablation MC-ICPMS boron isotope measurements in CaCO3 and s
ilicate matrices. JAAS 36:1607.*

(doi: 10.1039/d1ja00073j)

This code stemmed from David Evans original version in Matlab: (https://github.com/dbjevans/LA_d11B). Go check it out!
Although these two script are using the same methodology, there lie a slight difference of less than 0.1‰. 

The Python version uses a OO approach, and is thus best used in an interactive environment (e.g. Jupyter Notebook, IPython console, ...)

Here is a short outline of the steps taken by the script:

1. Go through the .exp files and laser log files to make them into a usable format.
2. Calculate the time offset between the laser and MS to detect where the samples are.
3. Remove outliers and substract background from analysis.
4. Perform a primary standard bracketing mass bias correction using NIST612.
5. Calculate 'raw' $\delta ^{11}\mathrm{B}$, [B] and B/Ca.
6. Calculate the inacurracy ($\Delta \delta ^{11} \mathrm{B}$) of three well studied carbonate standards, namely MACS-3, JCt and JCp. 
7. Derive a power relationship between observed inacurracy as a function of ~B/Ca (11/10.035) for this session.
8. Use this relationship to correct the samples.
9. Check accuracy of three secondary standards
10. Export the reduced data.

The script is still undergoing constant improvements and feature additions.

A full working example with selected files can be found in `./example/` as a Jupyter Notebook.

If you want an interactive GUI for selecting the files, you need the `tkinter` package.