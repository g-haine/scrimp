# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Ensure the step-by-step Sphinx code IS the qualified wave1d model.

The first tutorial must not silently drift from examples/wave_1d.py:
all public SCRIMP constructors, GWFL bricks, control values, initial
data, actual DAE/CN/LU options and energy/plot calls must match
the original function's Python AST. This is source-only, and cannot
substitute for the already completed five-backend native acceptance.

The pedagogical recipe uses the same names as the source, but is
split across consecutive RST blocks; we parse those blocks as
one complete Python module to prove they are copy/paste coherent.
"""

import ast
from collections import Counter
from pathlib import Path
import re
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "source" / "started.rst"
SOURCE = ROOT / "examples" / "wave_1d.py"
SCRIMP_CALLS = frozenset({
    "Domain", "State", "CoState", "Control_Port", "FEM",
    "Parameter", "Brick", "Term",
})
SYSTEM_CALLS = frozenset({
    "set_domain", "add_state", "add_costate", "add_control_port",
    "add_FEM", "add_parameter", "add_brick", "set_control",
    "set_initial_value", "set_time_scheme", "solve",
    "plot_Hamiltonian",
})
HAMILTONIAN_CALLS = frozenset({"set_name", "add_term"})


def code_blocks_from_coding_section():
    text = DOC.read_text(encoding="utf-8")
    start = text.index("Coding within SCRIMP\n")
    end = text.index(
        "The script performs the check explicitly from the accepted FE",
        start,
    )
    lines = text[start:end].splitlines()
    blocks = []
    index = 0
    while index < len(lines):
        if re.match(r"^\.\. code(?:-block)?:: python\s*$", lines[index]):
            index += 1
            # RST requires one empty separator before indented code.
            while index < len(lines) and not lines[index].strip():
                index += 1
            block = []
            while index < len(lines):
                line = lines[index]
                if line.startswith("    ") or not line.strip():
                    block.append(line)
                    index += 1
                else:
                    break
            blocks.append(textwrap.dedent("\n".join(block)).strip())
        else:
            index += 1
    if len(blocks) < 8:
        raise RuntimeError("Missing progressive Sphinx tutorial snippets")
    return blocks


def source_function_ast():
    source = ast.parse(SOURCE.read_text(encoding="utf-8"))
    return next(
        node for node in source.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "wave_1d"
    )


def relevant_calls(tree):
    """Keep the AST of actual constructors/calls, not loose substrings."""
    result = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(
            node.func, ast.Attribute
        ):
            continue
        function = node.func
        if isinstance(function.value, ast.Name):
            if function.value.id == "S" and function.attr in SCRIMP_CALLS:
                result.append(("S."+function.attr, ast.dump(node)))
            if (function.value.id == "system"
                    and function.attr in SYSTEM_CALLS):
                result.append(("system."+function.attr, ast.dump(node)))
        elif (
            isinstance(function.value, ast.Attribute)
            and isinstance(function.value.value, ast.Name)
            and function.value.value.id == "system"
            and function.value.attr == "hamiltonian"
            and function.attr in HAMILTONIAN_CALLS
        ):
            result.append(("system.hamiltonian."+function.attr, ast.dump(node)))
    return Counter(result)


class Wave1DDocCodeParity(unittest.TestCase):
    def test_every_python_snippet_parses_and_recipe_compiles(self):
        blocks = code_blocks_from_coding_section()
        for number, block in enumerate(blocks, 1):
            with self.subTest(block=number):
                ast.parse(block, filename=f"started.rst code block {number}")
        joined = "\n\n".join(blocks) + "\n"
        compile(joined, str(DOC), "exec")

    def test_all_scrimp_calls_match_qualified_python_exactly(self):
        source = source_function_ast()
        blocks = code_blocks_from_coding_section()
        recipe = ast.parse("\n\n".join(blocks), filename=str(DOC))
        native_calls = relevant_calls(source)
        doc_calls = relevant_calls(recipe)
        self.assertTrue(native_calls, "No real calls in source")
        self.assertEqual(
            doc_calls, native_calls,
            "Sphinx constructors, GWFL expressions, ports or solver calls "
            "differ from the native-qualified examples/wave_1d.py",
        )

    def test_source_default_parameters_match_copyable_tutorial(self):
        source = source_function_ast()
        defaults = {
            argument.arg: ast.literal_eval(value)
            for argument, value in zip(
                source.args.kwonlyargs, source.args.kw_defaults
            )
        }
        self.assertEqual(defaults, {
            "final_time": 0.04,
            "time_step": 0.004,
            "mesh_size": 0.05,
            "factor_solver": None,
            "plot": False,
        })
        blocks = code_blocks_from_coding_section()
        recipe = ast.parse("\n\n".join(blocks), filename=str(DOC))
        assignments = {}
        for statement in recipe.body:
            if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
                target = statement.targets[0]
                if isinstance(target, ast.Name) and target.id in defaults:
                    assignments[target.id] = ast.literal_eval(statement.value)
        self.assertEqual(
            assignments, defaults,
            "The RST recipe must declare the tested values of every "
            "wave_1d() user parameter",
        )

    def test_historical_1d_source_is_retired(self):
        self.assertFalse((ROOT/"sandbox"/"wave_1D.py").exists())
        self.assertIn(
            "examples/wave_1d.py", DOC.read_text(encoding="utf-8")
        )


if __name__ == "__main__":
    unittest.main()
