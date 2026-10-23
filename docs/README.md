<!--
  ~ Copyright 2026 QHARM Team
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~ http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
-->

Building the Documentation
==========================

The documentation is written with Sphinx. The Sphinx build requires several
Python packages (see the `requirements.txt` file), all of which are available
from PIP. The workflow below shows how to set up a Python virtual environment,
install the dependencies, and build the documentation:

```.bash
# In the docs/ directory
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
make html
```

You can then point your web browser at `docs/build/html/index.html` to view
the documentation.

To check the documentation the same way CI does (warnings are treated as
errors), run:

```.bash
make html SPHINXOPTS="-W --keep-going -n"
```
