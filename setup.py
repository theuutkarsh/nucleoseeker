from setuptools import setup, find_packages

# Read the contents of requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='nucleoseeker',
    version='0.1.1',
    author='Utkarsh Upadhyay',
    author_email='u.upadhyay@fz-juelich.de',
    description='Precision filtering of RNA databases to curate high-quality datasets',
    long_description=open('docs/README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/theuutkarsh/nucleoseeker',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'nucleoseeker = nucleoseeker.dataset_creator:main'
        ],
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',  # For Unix-based systems
        'Operating System :: POSIX :: Linux',  # Specifically for Linux
        'Operating System :: MacOS :: MacOS X',  # Specifically for macOS
        'Operating System :: Unix',  # Generic Unix OS
        'Programming Language :: Python :: 3.12',
    ],
    python_requires='>=3.6',
    install_requires=requirements,  # Automatically load dependencies from requirements.txt
)
