GetFEM GWFL matrix notation in Python
=====================================

.. _gwfl-matrix-convention:

.. warning::

   **A GetFEM expression string is not a Python or NumPy matrix.**
   In the Generic Weak-Form Language (GWFL), the nested literal
   ``[[a,b],[c,d]]`` contains **two column vectors**. In Python,
   the apparently identical list ``[[a,b],[c,d]]`` is normally
   interpreted as **two rows**.

   This difference can silently reverse the sign of a rotation,
   alter a boundary control, and invalidate the physical initial
   conditions even when the assembled Jacobian is correct.

   **Official documentation:** `GetFEM — Generic Weak-Form
   Language, "Explicit matrices"
   <https://getfem.org/userdoc/gasm_high.html#explicit-matrices>`_.
   GetFEM explicitly describes the nested entries as column vectors
   and gives the equivalent notation using semicolons to separate rows.

SCRIMP runs in Python, but weak-form expressions given to
``S.Brick`` and ``gf_model.add_macro`` are **GWFL strings**,
not Python expressions. Consequently:

.. code-block:: python

   import numpy as np

   # Row-wise Python/NumPy matrix: clockwise rotation.
   R_numpy = np.array([[0, 1], [-1, 0]])

   # WRONG if interpreted as the same matrix:
   # GWFL treats the nested lists in this string as COLUMNS.
   R_wrong_gwfl = "[[0,1],[-1,0]]"

For the intended clockwise matrix

.. math::

   R = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix},

write the following in SCRIMP:

.. code-block:: python

   # Correct nested GWFL columns:
   dphs.gf_model.add_macro("Rot", "[[0,-1],[1,0]]")

   # Equivalent, less ambiguous GWFL row/semicolon syntax:
   dphs.gf_model.add_macro("Rot", "[0,1;-1,0]")

.. tip::

   **Prefer the row/semicolon form** ``[a,b;c,d]`` when transcribing
   a matrix from a Python/NumPy calculation or an article.
   The semicolon-separated rows follow the familiar visual order.
   When using the nested form, supply the **columns**, not the rows,
   of the mathematical matrix.

Sign-check a boundary tangent
-----------------------------

The GWFL keyword ``Normal`` is the outward unit normal
on an external boundary (see the `GetFEM reference for Normal
<https://getfem.org/userdoc/gasm_high.html#special-expressions-linked-to-the-current-position>`_).

With the corrected rotation:

.. code-block:: python

   dphs.gf_model.add_macro("Rot", "[0,1;-1,0]")
   dphs.gf_model.add_macro("Tangent", "(Rot*Normal)")

at the rightmost point of the unit circle :math:`(1,0)`,
the normal :math:`n=(1,0)` must be mapped to the
clockwise tangent :math:`Rn=(0,-1)` (pointing down).

A Python source test can check the literal matrix convention,
but a **native GetFEM** test should also check the actual
``Normal`` and boundary integration whenever geometric mappings,
regions or MPI assembly could influence the result.

Native #28 example: viscous Saint-Venant
----------------------------------------

On 20 November 2026, native GetFEM boundary integrals exposed
a wrong-sign rotating wall control in the viscous Saint-Venant
sandbox. The prescribed initial momentum was clockwise, but the
old nested GWFL rotation macro acted counterclockwise.

Correcting the macro from ``"[[0,1],[-1,0]]"`` to
``"[[0,-1],[1,0]]"`` restored the original intended
orientation. In the quadratic (P2 **geometric**) disk and on
the same prescribed initial state:

* The relative velocity/tangent trace mismatch fell from
  approximately **200%** to **0.00697%**.
* The fixed-differential-state algebraic LSQR residual fell
  from **0.02773364** to **4.1851e-8**
  (more than **660,000-fold**).
* Both sign and residual checks agree across MPI 1 and 2.

These results **do not** qualify a complete nonlinear DAE
initializer. The projected rotating residual is still above
the strict numerical acceptance tolerance
(:math:`2.61\times 10^{-9}`), and the resting control has a
small nonzero LSQR residual despite a known equilibrium.
Further DAE initialization work is therefore still needed;
see :doc:`initialization`.

The gyroscopic expression
``Gyro(p) = div(Rot*p)*Rot`` is unchanged under a
consistent sign reversal of both instances of ``Rot``:
the two signs cancel. **Do not assume** that this holds
for other formulas containing only one rotation operator.

For the complete native evidence and development status,
see ``docs/V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md``
and GitHub issue #28.
