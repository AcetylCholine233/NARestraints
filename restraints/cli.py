
from __future__ import annotations

import argparse

from .base_pairs import read_base_pair_file
from .builder import build_phil_from_pdb
from .metal_pairs import DEFAULT_MAX_DISTANCE
from .validator import validate_and_fix_pdb


def add_metal_pair_arguments(parser: argparse.ArgumentParser) -> None:
    """Options shared by narestraints and narestraints-guesser."""
    group = parser.add_argument_group(
        "metal-mediated base pairs",
        "A metal ion (HETATM) within bonding distance of base N/O/S/Se atoms of two "
        "non-neighbouring nucleotides marks a metal-mediated pair. Its bases get no "
        "base-pair restraints and are never guessed as a canonical pair.",
    )
    group.add_argument(
        "--metal-distance", type=float, default=DEFAULT_MAX_DISTANCE, metavar="A",
        help=f"maximum metal-to-base-atom distance in Angstrom (default {DEFAULT_MAX_DISTANCE})",
    )
    group.add_argument(
        "--metal-elements", metavar="LIST", default=None,
        help="comma-separated metals to consider, e.g. AG,HG "
             "(default: metals with atomic number 21/Sc and above)",
    )
    group.add_argument(
        "--no-metal-pairs", action="store_true",
        help="do not look for metal-mediated pairs (previous behaviour)",
    )


def metal_elements_from_args(parser: argparse.ArgumentParser, args):
    if args.metal_distance <= 0:
        parser.error("--metal-distance must be positive")
    if args.metal_elements is None:
        return None
    from .anomalous import parse_element_list

    try:
        return parse_element_list(args.metal_elements)
    except ValueError as exc:
        parser.error(f"--metal-elements: {exc}")


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

    add_metal_pair_arguments(parser)
    parser.add_argument(
        "--keep-metal-pair-stacking", action="store_true",
        help="keep stacking restraints for bases in metal-mediated pairs "
             "(their base-pair restraints are still skipped)",
    )

    args = parser.parse_args()
    metal_elements = metal_elements_from_args(parser, args)

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
        metal_pairs=not args.no_metal_pairs,
        metal_distance=args.metal_distance,
        metal_elements=metal_elements,
        keep_metal_pair_stacking=args.keep_metal_pair_stacking,
    )

    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
