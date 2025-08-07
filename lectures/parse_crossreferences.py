import sys
import re

if len(sys.argv) < 2:
    print("Usage: python parse_crossreferences.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, "r", encoding="utf8") as f:
    content = f.read()
    # Regular expression to match cross-references
    matches = re.findall(":(fig:[a-zA-Z0-9_]+):", content)
    r = {f"{figure}": f"Figure {j+1}" for j, figure in enumerate(matches)}
    for reference, replacement in r.items():
        content = re.sub(f"[:@]{reference}", replacement, content)

with open(filename, "w", encoding="utf8") as f:
    f.write(content)
