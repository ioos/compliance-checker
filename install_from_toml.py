import subprocess
from itertools import chain

import toml

f = toml.load("pyproject.toml")

deps = list(chain.from_iterable(f["project"]["optional-dependencies"].values())) + f["build-system"]["requires"] + f["project"]["dependencies"]
deps = [dep.split(";")[0] for dep in deps if "=='emscripten'" not in dep]

subprocess.call(["micromamba", "install"] + deps)
