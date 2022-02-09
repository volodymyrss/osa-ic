from setuptools import setup

setup(name='osa-ic',
      version='0.1.0',
      packages=find_packages(),
      entry_points={
                    "console_scripts": [
                        "osa-ic = osaic.integralicindex:cli"
                    ]
                   },
      )
