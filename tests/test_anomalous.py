from pathlib import Path

import pytest

from restraints.anomalous import (
    ELEMENT_COLUMN,
    INFERRED_FROM_ATOM_NAME,
    INFERRED_FROM_ION_RESIDUE,
    anomalous_phil_block,
    find_heavy_atoms,
    parse_element_list,
)
from restraints.phenix import write_phil


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def hetatm(
    serial, name, resname, chain, resseq, element="", *, altloc=" ", icode=" ",
    record="HETATM", occupancy=1.0,
):
    """Build one fixed-column PDB coordinate record."""
    return (
        f"{record:<6}{serial:>5} {name:<4}{altloc}{resname:>3} {chain}{resseq:>4}{icode}"
        f"   {0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{occupancy:>6.2f}{20.0:>6.2f}"
        f"          {element:>2}"
    )


def write_pdb(tmp_path, lines, name="model.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_record_helper_uses_pdb_columns():
    line = hetatm(1, "AG", "AG", "Z", 6, "AG")
    assert line[0:6] == "HETATM"
    assert line[12:16] == "AG  "
    assert line[17:20] == " AG"
    assert line[21] == "Z"
    assert line[22:26] == "   6"
    assert line[76:78] == "AG"


def test_detects_metal_from_element_column_and_ignores_atom_records(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "N1", "DG", "A", 1, "N", record="ATOM"),
        hetatm(2, "AG", "AG", "Z", 6, "AG"),
    ])
    sites, warnings = find_heavy_atoms(pdb)
    assert warnings == []
    assert [(s.element, s.chain, s.resid, s.atom_name) for s in sites] == [("AG", "Z", "6", "AG")]
    assert sites[0].element_source == ELEMENT_COLUMN
    assert sites[0].selection() == "chain Z and resid 6 and name AG"


def test_default_cutoff_skips_water_light_ions_and_ligand_atoms(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "O", "HOH", "A", 101, "O"),
        hetatm(2, "MG", "MG", "A", 102, "MG"),
        hetatm(3, "NA", "NA", "A", 103, "NA"),
        hetatm(4, "K", "K", "A", 104, "K"),
        hetatm(5, "CA", "CA", "A", 105, "CA"),
        hetatm(6, "P", "DF", "B", 4, "P"),
        hetatm(7, "S1", "DF", "B", 4, "S"),
        hetatm(8, "ZN", "ZN", "A", 106, "ZN"),
    ])
    sites, _ = find_heavy_atoms(pdb)
    assert [s.element for s in sites] == ["ZN"]


def test_heavy_atom_inside_modified_nucleotide_is_found(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "C5", "5IU", "B", 8, "C"),
        hetatm(2, "I5", "5IU", "B", 8, "I"),
        hetatm(3, "BR", "BRU", "B", 9, "BR"),
        hetatm(4, "SE", "SE4", "B", 10, "SE"),
    ])
    sites, _ = find_heavy_atoms(pdb)
    assert [(s.element, s.residue_name, s.atom_name) for s in sites] == [
        ("I", "5IU", "I5"),
        ("BR", "BRU", "BR"),
        ("SE", "SE4", "SE"),
    ]


def test_explicit_element_list_replaces_the_cutoff(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "MG", "MG", "A", 102, "MG"),
        hetatm(2, "ZN", "ZN", "A", 106, "ZN"),
        hetatm(3, "AG", "AG", "A", 107, "AG"),
    ])
    sites, _ = find_heavy_atoms(pdb, elements="mg, Ag")
    assert [s.element for s in sites] == ["MG", "AG"]


def test_parse_element_list_rejects_non_elements():
    assert parse_element_list("AG,hg;Br") == frozenset({"AG", "HG", "BR"})
    with pytest.raises(ValueError, match="Not an element"):
        parse_element_list("AG,XX")
    with pytest.raises(ValueError, match="No element"):
        parse_element_list(" , ")


def test_blank_element_column_inferred_from_ion_residue(tmp_path):
    # " CA " follows the one-letter convention and would read as carbon, but a
    # CA residue whose atom is CA is calcium.
    pdb = write_pdb(tmp_path, [
        hetatm(1, " CA", "CA", "A", 1, ""),
        hetatm(2, " AG", "AG", "A", 2, ""),
    ])
    sites, warnings = find_heavy_atoms(pdb, elements="CA,AG")
    assert warnings == []
    assert [(s.element, s.element_source) for s in sites] == [
        ("CA", INFERRED_FROM_ION_RESIDUE),
        ("AG", INFERRED_FROM_ION_RESIDUE),
    ]


def test_blank_element_column_inferred_from_atom_name(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, " C1", "LIG", "A", 1, ""),   # one-letter convention: carbon
        hetatm(2, "HG1", "LIG", "A", 1, ""),   # two-letter convention: mercury
        hetatm(3, " I5", "5IU", "B", 8, ""),
    ])
    sites, warnings = find_heavy_atoms(pdb)
    assert warnings == []
    assert [(s.element, s.atom_name, s.element_source) for s in sites] == [
        ("HG", "HG1", INFERRED_FROM_ATOM_NAME),
        ("I", "I5", INFERRED_FROM_ATOM_NAME),
    ]


def test_undeterminable_element_is_warned_and_skipped(tmp_path):
    pdb = write_pdb(tmp_path, [hetatm(1, "QQ", "UNK", "A", 1, "")])
    sites, warnings = find_heavy_atoms(pdb)
    assert sites == []
    assert len(warnings) == 1 and "could not determine the element" in warnings[0]


def test_only_first_model_is_read(tmp_path):
    pdb = write_pdb(tmp_path, [
        "MODEL        1",
        hetatm(1, "AG", "AG", "Z", 6, "AG"),
        "ENDMDL",
        "MODEL        2",
        hetatm(2, "HG", "HG", "Z", 7, "HG"),
        "ENDMDL",
    ])
    sites, _ = find_heavy_atoms(pdb)
    assert [s.element for s in sites] == ["AG"]


def test_selection_handles_insertion_code_and_blank_chain(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "AG", "AG", "A", 52, "AG", icode="A"),
        hetatm(2, "AG", "AG", " ", 60, "AG"),
    ])
    sites, _ = find_heavy_atoms(pdb)
    assert [s.selection() for s in sites] == [
        "chain A and resid 52A and name AG",
        "resid 60 and name AG",
    ]


def test_phil_block_groups_by_element_and_merges_alternate_locations(tmp_path):
    pdb = write_pdb(tmp_path, [
        hetatm(1, "AG", "AG", "Z", 6, "AG", altloc="A", occupancy=0.5),
        hetatm(2, "AG", "AG", "Z", 6, "AG", altloc="B", occupancy=0.5),
        hetatm(3, "AG", "AG", "Z", 7, "AG"),
        hetatm(4, "HG", "HG", "Z", 8, "HG"),
    ])
    sites, _ = find_heavy_atoms(pdb)
    block = anomalous_phil_block(sites)

    assert block.startswith("  refine {")
    assert (
        "    strategy = *individual_sites *individual_sites_real_space "
        "*individual_adp *occupancies *group_anomalous"
    ) in block
    assert block.count("group {") == 2
    assert (
        'selection = "(chain Z and resid 6 and name AG) or '
        '(chain Z and resid 7 and name AG)"'
    ) in block
    assert 'selection = "chain Z and resid 8 and name HG"' in block
    assert block.count("f_prime = 0") == 2
    assert block.count("f_double_prime = 0") == 2
    assert block.count("refine = *f_prime *f_double_prime") == 2
    assert block.count("{") == block.count("}")


def test_phil_block_without_strategy_and_with_no_sites(tmp_path):
    pdb = write_pdb(tmp_path, [hetatm(1, "AG", "AG", "Z", 6, "AG")])
    sites, _ = find_heavy_atoms(pdb)
    assert "strategy" not in anomalous_phil_block(sites, include_strategy=False)
    # No groups: write nothing, because group_anomalous without groups is an error.
    assert anomalous_phil_block([]) == ""


def test_write_phil_default_output_is_unchanged(tmp_path):
    before = tmp_path / "before.phil"
    write_phil(before, ["    bond {\n    }"], "        stacking_pair {\n        }")
    text = before.read_text()
    assert text.endswith("    bond {\n    }\n  }\n}\n")
    assert "anomalous" not in text and "refine {" not in text


def test_write_phil_places_anomalous_scope_inside_refinement(tmp_path):
    pdb = write_pdb(tmp_path, [hetatm(1, "AG", "AG", "Z", 6, "AG")])
    sites, _ = find_heavy_atoms(pdb)
    out = tmp_path / "out.phil"
    write_phil(out, [], "", anomalous_phil_block(sites))
    text = out.read_text()

    assert text.startswith("refinement {\n")
    assert text.endswith("  }\n}\n")
    edits = text.index("  geometry_restraints.edits {")
    refine = text.index("  refine {")
    assert edits < refine
    assert text.count("{") == text.count("}")


def test_bundled_examples():
    # dfPucker.pdb writes silver with left-justified names and element "Ag".
    sites, _ = find_heavy_atoms(EXAMPLES / "dfPucker.pdb")
    assert [(s.element, s.chain, s.resid) for s in sites] == [("AG", "Z", "6")]

    # 2j02.pdb: two Zn ions are heavy; its 638 Mg ions are not, unless asked for.
    sites, _ = find_heavy_atoms(EXAMPLES / "2j02.pdb")
    assert [s.element for s in sites] == ["ZN", "ZN"]
    sites, _ = find_heavy_atoms(EXAMPLES / "2j02.pdb", elements="MG")
    assert len(sites) == 638

    # 1BNA.pdb has only waters as HETATM.
    sites, _ = find_heavy_atoms(EXAMPLES / "1BNA.pdb")
    assert sites == []


def test_builder_writes_anomalous_groups_end_to_end(tmp_path):
    pytest.importorskip("Bio")
    pytest.importorskip("openpyxl")
    from restraints.builder import build_phil_from_pdb

    # No base pairs: this test is about wiring the toggle through the builder;
    # pair and stacking recipes are covered by their own tests.
    out = tmp_path / "dfPucker.phil"
    stretches = []
    sites = build_phil_from_pdb(
        EXAMPLES / "dfPucker-fixed.pdb", stretches, out, anomalous=True
    )
    text = out.read_text()

    assert [s.element for s in sites] == ["AG"]
    assert "geometry_restraints.edits" in text
    assert "*group_anomalous" in text
    assert 'selection = "chain Z and resid 6 and name AG"' in text
    assert text.count("{") == text.count("}")

    # Without the toggle, the output has no anomalous scope at all.
    plain = tmp_path / "plain.phil"
    assert build_phil_from_pdb(EXAMPLES / "dfPucker-fixed.pdb", stretches, plain) == []
    assert "anomalous" not in plain.read_text()
