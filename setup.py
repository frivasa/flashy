from setuptools import setup, find_packages

setup(name='flashy',
      version='0.4',
      description='FLASH output analysis tools',
      url='http://github.com/frivasa/flashy',
      author='Fernando Rivas',
      author_email='rivas.aguilera@gmail.com',
      license='MIT',
      packages=find_packages(exclude=['docs']),
#       package_dir={'':'flashy'},
      zip_safe=False)
