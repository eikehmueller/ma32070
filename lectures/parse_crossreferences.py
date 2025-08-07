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
        content = re.sub(
            f":{reference}", f'<a name="{reference}">{replacement}</a>', content
        )
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
