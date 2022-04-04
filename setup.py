from setuptools import find_packages, setup

setup(
    name="src",
    version="0.0.1",
    author="Hamish Campbell",
    author_email="author@example.com",
    description="Can artificial intelligence methods be used to improve the assessment of how hunting pressure effects the biodiversity co-benefits of carbon credit projects? I aim to use a CNN model to develop a species distribution model and look for anomolies that can be atributed to hunting pressure.",
    url="url-to-github-page",
    packages=find_packages(),
    test_suite="src.tests.test_all.suite",
)
