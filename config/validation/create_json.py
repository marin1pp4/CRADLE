import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "custom_coefficients.json"

COEFFICIENTS = ["a", "b", "c", "A", "B", "D"]

with open(SOURCE) as f:
    expressions = json.load(f)["11C"]

# Verify that all general expressions are available.
for coeff in COEFFICIENTS:
    if coeff not in expressions or not expressions[coeff].strip():
        raise ValueError(f"Missing expression for coefficient {coeff}")

    if expressions[coeff].strip() == "0":
        raise ValueError(
            f"Coefficient {coeff} is zero in the source JSON. "
            "Restore its general expression before generating the files."
        )

# One nonzero coefficient per benchmark.
for coeff in COEFFICIENTS:

    config = {
        "11C": {name: "0" for name in COEFFICIENTS}
    }

    config["11C"][coeff] = expressions[coeff]

    output = ROOT / f"coefficients_{coeff}.json"

    with open(output, "w") as f:
        json.dump(config, f, indent=4)

    print(f"Created: {output.name}")

# Reference sample for b = 0.
reference = {
    "11C": {name: "0" for name in COEFFICIENTS}
}

with open(ROOT / "coefficients_b_ref.json", "w") as f:
    json.dump(reference, f, indent=4)

print("Created: coefficients_b_ref.json")
