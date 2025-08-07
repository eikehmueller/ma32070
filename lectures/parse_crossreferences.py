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
        content = re.sub(
            f":{reference}", f'<a name="reference">{replacement}</a>', content
        )
        content = re.sub(
            f"@{reference}", f'<a href="#reference">{replacement}</a>', content
        )

print(content)

"""
Take me to <a href="#pookie">pookie</a>

... 

<a name="pookie">this is pookie</a>
"""
