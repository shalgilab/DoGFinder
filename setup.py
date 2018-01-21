import os
from setuptools import setup
from distutils.core import setup
from setuptools.command.install import install
import commands
import sys

class CustomInstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):
        bed_strng=commands.getstatusoutput('bedtools--version')
        if bed_strng[1].find("command not found")!=-1:
            print "Warning: Couldn't find bedtools on path"
        sam_strng=commands.getstatusoutput('samtools')
        
        if sam_strng[1].find("command not found")!=-1:
            print "Warning: Couldn't find samtools on path"
        install.run(self)

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
        
setup(name='DoGFinder',
      version='1.0.0',
      description='DoGs pipeline',
      install_requires=['pybedtools','numpy','pandas','pysam','RSeQC'],
      scripts=['Get_loci_annotation','Get_DoGs','Get_DoGs_rpkm','Pre_Process','Union_DoGs_annotation','Common_DoGs_annotation'],
      packages=['DoGs_functions','tests','annotations'],
      package_data={'tests': ['data/*'],'annotations': ['*']},
      cmdclass={'install': CustomInstallCommand},
      long_description=read('README.md'),
      test_suite="tests", 
      )


