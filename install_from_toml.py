import subprocess
from itertools import chain

import toml

f = toml.load("pyproject.toml")

deps = list(chain.from_iterable(f["project"]["optional-dependencies"].values())) + f["build-system"]["requires"] + f["project"]["dependencies"]
deps = [v.split(";")[0] for v in deps]

subprocess.call(["micromamba", "install"] + [v.split(";")[0] for v in deps])
