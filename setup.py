import distutils
from distutils.core import setup
import glob
import os

bin_files = glob.glob("bin/*")

# Build the structure of bpz/etc
etc_files = {}
for root, dirs, files in os.walk("etc", topdown=False):
    for name in files:
        try:
            etc_files[root].append(os.path.join(root, name))
        except:
            etc_files[root] = [os.path.join(root, name)]

data_files = [(k,v) for k, v in etc_files.iteritems()]

# The main call
setup(name='bpz',
      version ='1.0',
      license = "GPL",
      description = "DESDM implementation of BPZ",
      author = "Ben Hoyle, Felipe Menanteau",
      author_email = "benhoyle1212@gmail.com, felipe@illinois.edu",
      packages = ['bpz'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=data_files,
      )
