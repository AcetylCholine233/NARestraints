
from __future__ import annotations

import argparse

from .base_pairs import read_base_pair_file
from .builder import build_phil_from_pdb
from .validator import validate_and_fix_pdb


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Phenix nucleic-acid restraints from a PDB and base-pair range file."
    )
    parser.add_argument("pdb", help="Input PDB file")
    parser.add_argument("pairs", help="Base-pair specification text file")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output .phil filename (default: <PDB stem>.phil)",
    )
    parser.add_argument("--no-parallels", action="store_true")
    parser.add_argument("--no-planes", action="store_true")
    parser.add_argument("--no-stacking", action="store_true")

    anomalous = parser.add_argument_group(
        "anomalous scatterers",
        "Detect heavy atoms in HETATM records (waters excluded) and write them as "
        "phenix.refine anomalous scatterer groups, one per element, with f' and f'' "
        "starting at 0 and refined.",
    )
    anomalous.add_argument(
        "--anomalous", action="store_true",
        help="detect heavy atoms (atomic number 21/Sc and above) in HETATM records",
    )
    anomalous.add_argument(
        "--anomalous-elements", metavar="LIST", default=None,
        help="comma-separated elements to use instead of the default cutoff, "
             "e.g. AG or AG,HG,BR (implies --anomalous)",
    )
    anomalous.add_argument(
        "--no-anomalous-strategy", action="store_true",
        help="do not add group_anomalous to the refinement strategy (use when the "
             "strategy is set elsewhere, e.g. by NASolve or on the command line)",
    )

    args = parser.parse_args()

    anomalous_elements = None
    if args.anomalous_elements is not None:
        from .anomalous import parse_element_list

        try:
            anomalous_elements = parse_element_list(args.anomalous_elements)
        except ValueError as exc:
            parser.error(f"--anomalous-elements: {exc}")

    output = args.output
    if output is None:
        output = args.pdb.rsplit(".", 1)[0] + ".phil"

    pdb = validate_and_fix_pdb(args.pdb)
    stretches = read_base_pair_file(args.pairs)

    build_phil_from_pdb(
        pdb,
        stretches,
        output,
        parallels=not args.no_parallels,
        planes=not args.no_planes,
        include_stacking=not args.no_stacking,
        anomalous=args.anomalous,
        anomalous_elements=anomalous_elements,
        anomalous_strategy=not args.no_anomalous_strategy,
    )

    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
