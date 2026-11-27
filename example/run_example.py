"""Run the ResistFlow example: estimate fitness parameters from example data.

   Usage (from the repository root):
       poetry run python example/run_example.py
"""

import json
from pathlib import Path

from resistflow import estimate_fitness_parameters, load_local_data

ROOT = Path(__file__).resolve().parent.parent
INPUT_CSV = ROOT / "example" / "input" / "example_data.csv"
OUTPUT_JSON = ROOT / "example" / "output" / "fitted_parameters.json"


def main():
    data = load_local_data(INPUT_CSV)

    params = estimate_fitness_parameters(
        p0=data["p0"],
        observed_years=data["observed_years"],
        observed_frequencies=data["observed_frequencies"],
        n_chromosomes=data["n_chromosomes"],
    )

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSON, "w", newline="\n") as f:
        json.dump(params, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"Wrote {OUTPUT_JSON}")


if __name__ == "__main__":
    main()

