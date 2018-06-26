Instructions for the execution of the stand-alone version:

This program is cross-platform. In order to run it you should have installed:
1. Python 3.6 (or higher)
2. MiXCR

In addition, the following Python modules are needed for the program to run: (try to run this line in your Python she
import Bio, difflib, logging, os, matplotlib, numpy, re, regex, seaborn, scipy, shutil, subprocess, sys, time, traceback, umi_tools, weblogolib

Once everything is set, you need to write a parameters file for the pipeline (see parameters.txt file for example)

Finally, run the following command in the terminal:
cd <pipeline folder>
python NGS_analyzer.py <path to a parameters file>

P.s., if you have more than one sample and you would like to run this script several times, you can write a script that iterates over the parameters files (you might also want to generate these automatically) and fetches a system call of the same format above (NGS_analyzer.py execution).

If you have any question/thought/suggestion(s), please do contact us at bioSequence@tauex.tau.ac.il
