import os
from setuptools import setup
from distutils.core import setup
from setuptools.command.install import install
import sys

class CustomInstallCommand(install):
    """Customized setuptools install command - prints a friendly greeting."""
    def run(self):

        python_ver=sys.version
        if python_ver[0]=='3':
            print ("Error: you are using python 3, should be using python 2")
            sys.exit()
        bed_strng=os.popen('bedtools --version').read()
        if bed_strng.find("command not found")!=-1:
            print ("Warning: Couldn't find bedtools on path")
        sam_strng=os.popen('samtools').read()
        
        if sam_strng.find("not found")!=-1:
            print ("Warning: Couldn't find samtools on path")
        install.run(self)

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
        
setup(name='DoGFinder',
      version='1.0.1',
      description='DoGs pipeline',
      install_requires=['pybedtools','numpy','pandas','pysam','RSeQC==2.6.6'],
      scripts=['Get_loci_annotation','Get_DoGs','Get_DoGs_rpkm','Pre_Process','Union_DoGs_annotation','Common_DoGs_annotation'],
      packages=['DoGs_functions','tests','annotations'],
      package_data={'tests': ['data/*'],'annotations': ['*']},
      cmdclass={'install': CustomInstallCommand},
      long_description=read('README.md'),
      test_suite="tests", 
      )


