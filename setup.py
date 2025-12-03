from setuptools import setup, find_packages

setup(
      name='pyviblum',
      version='1.0',
      description='Code to calculate vibronic structure of luminescence in lanthanide complexes',
      author='Vsevolod D. Dergachev, Liviu F. Chibotaru, Sergey A. Varganov',
      url='https://github.com/svarganov/pyVibLum',
      license='MIT',
      package_dir={"": "src"},
      packages=find_packages(where="src"),
      install_requires=['scipy', 'numpy'],
      )
