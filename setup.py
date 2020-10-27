from setuptools import setup
setup(name='osa-ic',
      version='0.1.0',
      py_modules=['integralicindex','icdeposit'],
      entry_points={
                    "console_scripts": [
                        "osa-ic = integralicindex:main"
                    ]
                   },
      )
