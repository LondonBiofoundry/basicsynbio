import setuptools
import basicsynbio as app


def main():
    with open("README.md", "r") as fh:
        long_description = fh.read()

    setuptools.setup(
        name=app.__project__,
        version=app.__version__,
        url=app.__url__,
        license=next(
            (
                classifier.rsplit("::", 1)[1].strip()
                for classifier in app.__classifiers__
                if classifier.startswith("License ::")
            )
        ),
        author=app.__author__,
        author_email=app.__author_email__,
        description=app.__description__,
        long_description=long_description,
        long_description_content_type="text/markdown",
        packages=setuptools.find_packages(exclude=("tests")),
        include_package_data=True,
        classifiers=app.__classifiers__,
        python_requires=app.__python_version__,
        install_requires=app.__requires__,
        project_urls=app.__project_urls__,
    )


if __name__ == "__main__":
    main()
