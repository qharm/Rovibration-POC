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

import unittest

from rovib import Hamiltonian, HamiltonianTerm


class TestHamiltonian(unittest.TestCase):
    def setUp(self):
        self.t1 = HamiltonianTerm()
        self.t2 = HamiltonianTerm()
        self.t3 = HamiltonianTerm()

    def test_default(self):
        h = Hamiltonian()
        self.assertEqual(h.terms, ())
        self.assertEqual(len(h), 0)

    def test_terms(self):
        h = Hamiltonian([self.t1, self.t2])
        self.assertEqual(h.terms, (self.t1, self.t2))
        self.assertEqual(len(h), 2)
        self.assertEqual(list(h), [self.t1, self.t2])

    def test_non_term(self):
        with self.assertRaises(TypeError):
            Hamiltonian((self.t1, "not a term"))

    def test_add_term(self):
        h = Hamiltonian((self.t1,)) + self.t2
        self.assertEqual(h.terms, (self.t1, self.t2))

    def test_add_hamiltonian(self):
        h = Hamiltonian((self.t1,)) + Hamiltonian((self.t2, self.t3))
        self.assertEqual(h.terms, (self.t1, self.t2, self.t3))

    def test_chained_addition(self):
        h = self.t1 + self.t2 + self.t3
        self.assertEqual(h.terms, (self.t1, self.t2, self.t3))

    def test_add_does_not_modify_operands(self):
        h = Hamiltonian((self.t1,))
        _ = h + self.t2
        self.assertEqual(h.terms, (self.t1,))

    def test_add_unsupported(self):
        with self.assertRaises(TypeError):
            Hamiltonian() + 1
