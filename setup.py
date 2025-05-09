from setuptools import setup, find_packages

setup(
    name="metapilot",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.23.5",
        "pandas>=2.0.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
        "tqdm>=4.65.0",
        "torch>=2.0.0",
        "torchvision>=0.15.0",
        "torchaudio>=2.0.0",
        "tensorflow>=2.12.0",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A next-generation proteomics data analysis platform",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/metapilot",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.9",
) 