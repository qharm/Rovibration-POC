# Copyright 2026 QHARM Team
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Physical constants and unit conversions.

All calculations are done in atomic units. The values here are the CODATA
2018 recommended values and are only needed to convert inputs to, and
results from, atomic units.
"""

#: Conversion factor from hartree to wavenumbers (cm^-1).
HARTREE_TO_CM = 219474.6313632

#: Conversion factor from unified atomic mass units to electron masses.
AMU_TO_ME = 1822.888486209
