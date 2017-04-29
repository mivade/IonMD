from skbuild import setup


setup(
    name="ionmd",
    version="2.0.dev",
    author="Michael V. DePalatis",
    author_email="mike@depalatis.net",
    url="https://github.com/mivade/IonMD",
    description="Ion trap molecular dynamics simulation toolkit",
    # long_description='',
    cmake_args=["-DBUILD_PY=ON"],
    zip_safe=False,
)
