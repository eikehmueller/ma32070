"""Script for inserting html references into markdown

The script searches for all occurences of :fig:xxx and replaces them
with "Figure n", where n is the consecutive number of the figure, and
it inserts a html anchor at the same position. In addition, it
replaces all occurances of @fig:xxx by a html link to the anchor.

It also searches for all occurances of :eqn:xxx and replaces them
with "(n)", where n is the consecutive number of the equation, and it
inserts a html anchor just before the equation. In addition, it
replaces all occura
"""

import sys
import re

if len(sys.argv) < 2:
    print("Usage: python parse_crossreferences.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, "r", encoding="utf8") as f:
    content = f.read()
    # figures
    matches = re.findall(":(fig:[a-zA-Z0-9_]+)", content)
    r = {f"{figure}": f"Figure {j+1}" for j, figure in enumerate(matches)}
    for reference, replacement in r.items():
        pos = content.rfind("![", 0, content.find(f":{reference}"))
        content = content[:pos] + f'<a name="{reference}"></a>' + content[pos:]
        content = re.sub(f":{reference}", f"{replacement}", content)
        content = re.sub(
            f"@{reference}", f'<a href="#{reference}">{replacement}</a>', content
        )
    # equations
    matches = re.findall(":(eqn:[a-zA-Z0-9_]+)", content)
    r = {f"{equation}": f"({j+1})" for j, equation in enumerate(matches)}
    for reference, replacement in r.items():
        pos = content.rfind("$$", 0, content.find(f":{reference}"))
        content = content[:pos] + f'<a name="{reference}"></a>' + content[pos:]
        content = re.sub(
            f":{reference}",
            f"{replacement}",
            content,
        )
        content = re.sub(
            f"@{reference}", f'<a href="#{reference}">{replacement}</a>', content
        )

print(content)
