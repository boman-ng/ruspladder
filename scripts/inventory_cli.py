#!/usr/bin/env python3
"""Extract reference CLI declarations without importing or modifying SplAdder."""
import argparse
import ast
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    tree = ast.parse(args.source.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "parse_options")
    owners = {}
    result = {name: [] for name in ("prep", "build", "test")}
    for node in function.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            call = node.value
            if isinstance(call.func, ast.Attribute) and call.func.attr == "add_parser":
                owners[node.targets[0].id] = ast.literal_eval(call.args[0])
            elif isinstance(call.func, ast.Attribute) and call.func.attr == "add_argument_group":
                owners[node.targets[0].id] = owners.get(call.func.value.id)
        if not (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)):
            continue
        call = node.value
        if not (isinstance(call.func, ast.Attribute) and call.func.attr == "add_argument"):
            continue
        owner = owners.get(call.func.value.id)
        if owner not in result:
            continue
        entry = {"flags": [ast.literal_eval(a) for a in call.args], "source_line": node.lineno}
        for keyword in call.keywords:
            try: value = ast.literal_eval(keyword.value)
            except (ValueError, TypeError): value = ast.unparse(keyword.value)
            entry[keyword.arg] = value
        result[owner].append(entry)
    args.output.write_text(json.dumps({"reference": "SplAdder v3.1.1", "commands": result}, indent=2) + "\n")
    print({mode: len(entries) for mode, entries in result.items()})


if __name__ == "__main__":
    main()
