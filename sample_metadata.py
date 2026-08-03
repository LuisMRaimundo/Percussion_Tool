"""
Filename / folder metadata parsing for validation grouping.

METADATA-ONLY: never estimates physical parameters from audio.
Circularity refusal is a hard rule of the validation pipeline.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from model import PlateInstrument

# Catalogue entries keyed by exact inch size (crash-type Chladni anchors).
_CYMBAL_INCH_CATALOGUE: Dict[int, PlateInstrument] = {
    16: PlateInstrument(
        "cymbal_16in_thin", 0.406, 0.0008, chladni=(10.8, 2.0, 1.81)
    ),
    18: PlateInstrument(
        "cymbal_18in_medium", 0.457, 0.0012, chladni=(13.4, 2.0, 1.65)
    ),
}

# internal_default thickness per class (metres) when no exact catalogue row.
_THICKNESS_INTERNAL_DEFAULT = {
    "cymbal": 0.0012,   # medium-plate default (matches 18in / 46cm catalogue)
    "tamtam": 0.0015,   # matches nontunperc tamtam_80cm_bronze
    "windgong": 0.0015,
}

INCH_M = 0.0254

# Scientific pitch tokens in library names (Thai gong packs).
_NOTE_RE = re.compile(
    r"(?:^|[.\-_])([A-Ga-g](?:bb|b|#)?)([0-9])(?:[.\-_]|$)"
)

_STROKE_PATTERNS: Tuple[Tuple[str, str], ...] = (
    (r"stick\.bell", "stick.bell"),
    (r"stick\.normal", "stick.normal"),
    (r"stick\.shoulder", "stick.shoulder"),
    (r"\bbow\b", "bow"),
    (r"\bmallet\b", "mallet"),
)

_DYNAMIC_RE = re.compile(r"(?:^|[.\-_])(pp|mf|ff)(?:[.\-_]|$)", re.I)

# Size + type glued tokens: 17crash, 20windgong, 22tamtam, 16chinese, 21ride
_SIZED_TYPE_RE = re.compile(
    r"(?:^|[.\-_])(\d{2})(crash|chinese|ride|tamtam|windgong)(?:[.\-_]|$)",
    re.I,
)
_SPLASH_RE = re.compile(r"(?:^|[.\-_])splash(?:[.\-_]|$)", re.I)
_THAIGONG_RE = re.compile(r"thaigong", re.I)


@dataclass
class SampleMeta:
    """Parsed metadata for one sample path (no audio features)."""

    path: Path
    status: str  # "ok" | "skip_pitched" | "unparseable"
    instrument_id: Optional[str] = None  # group instrument key
    plate_class: Optional[str] = None  # cymbal | tamtam | windgong
    subtype: Optional[str] = None  # crash | chinese | splash | ride | …
    diameter_in: Optional[float] = None
    stroke: Optional[str] = None
    dynamic: Optional[str] = None
    skip_reason: Optional[str] = None
    notes: List[str] = field(default_factory=list)

    @property
    def group_key(self) -> Optional[Tuple[str, str, str]]:
        if self.status != "ok":
            return None
        if not (self.instrument_id and self.stroke and self.dynamic):
            return None
        return (self.instrument_id, self.stroke, self.dynamic)

    @property
    def is_transfer_caution(self) -> bool:
        """Chinese / splash / ride: Chladni crash anchors may not transfer."""
        return self.subtype in {"chinese", "splash", "ride"}

    @property
    def contributes_to_aggregate(self) -> bool:
        """PRIMARY aggregate excludes transfer-caution subtypes."""
        return (
            self.status == "ok"
            and not self.is_transfer_caution
            and self.dynamic in {"pp", "mf"}
        )


@dataclass
class ModelMapping:
    """Resolved PlateInstrument + provenance for one instrument_id."""

    instrument_id: str
    instrument: PlateInstrument
    provenance: str
    plate_class: str
    subtype: str
    diameter_in: float
    diameter_m: float
    thickness_m: float
    catalogue_match: Optional[str]
    transfer_caution: bool


def _path_text(path: Path) -> str:
    """Lowercased stem + parent folder names joined by dots."""
    parts = [path.stem] + [p for p in path.parts[:-1]]
    return ".".join(parts).replace("\\", ".").replace("/", ".").lower()


def parse_sample_path(path: Path) -> SampleMeta:
    """Parse size/type/stroke/dynamic from filename and folders only."""
    path = Path(path)
    text = _path_text(path)
    meta = SampleMeta(path=path, status="ok")

    # --- pitched Thai / note-named gongs: skip (validity scope) --------
    if _THAIGONG_RE.search(text) or (
        _NOTE_RE.search(text)
        and "gong" in text
        and "windgong" not in text
        and "tamtam" not in text
    ):
        meta.status = "skip_pitched"
        meta.skip_reason = (
            "tuned percussion (note name / Thai gong) — outside "
            "NonTunPerc validity scope; not forced onto an unpitched model"
        )
        note = _NOTE_RE.search(text)
        if note:
            meta.notes.append(f"note_token={note.group(1).upper()}{note.group(2)}")
        return meta

    # --- type + size ---------------------------------------------------
    diameter_in: Optional[float] = None
    subtype: Optional[str] = None
    plate_class: Optional[str] = None

    m = _SIZED_TYPE_RE.search(text)
    if m:
        diameter_in = float(m.group(1))
        token = m.group(2).lower()
        if token == "tamtam":
            plate_class, subtype = "tamtam", "tamtam"
        elif token == "windgong":
            plate_class, subtype = "windgong", "windgong"
        else:
            plate_class, subtype = "cymbal", token
    elif _SPLASH_RE.search(text):
        plate_class, subtype = "cymbal", "splash"
        # splash packs in this library often omit diameter → unparseable
    elif "tamtam" in text or "tam-tam" in text or "tam_tam" in text:
        plate_class, subtype = "tamtam", "tamtam"
    elif "windgong" in text or "wind-gong" in text:
        plate_class, subtype = "windgong", "windgong"
    elif "cymbal" in text or "crash" in text:
        plate_class = "cymbal"
        subtype = "crash" if "crash" in text else "cymbal"

    # Folder class hints when type token missing size glue
    if plate_class is None:
        if any(p.lower() in {"tam-tam", "tamtam"} for p in path.parts):
            plate_class, subtype = "tamtam", "tamtam"
        elif any(p.lower() == "cymbals" for p in path.parts):
            plate_class = "cymbal"
            subtype = subtype or "cymbal"

    if diameter_in is None:
        # e.g. folder 17crash already consumed; try bare NN before type words
        m2 = re.search(
            r"(?:^|[.\-_])(\d{2})(?:in)?(?:[.\-_]|$)", text
        )
        # Prefer sized-type; only accept bare size if class known and no note
        if m2 and plate_class in {"tamtam", "windgong", "cymbal"}:
            # Avoid grabbing dynamics-adjacent junk; require plausible inches
            val = int(m2.group(1))
            if 6 <= val <= 60:
                diameter_in = float(val)

    # --- stroke --------------------------------------------------------
    stroke: Optional[str] = None
    for pat, label in _STROKE_PATTERNS:
        if re.search(pat, text, re.I):
            stroke = label
            break
    if stroke is None:
        # Parsed absence — not a physical guess; groups stay separate.
        stroke = "unmarked"

    # --- dynamic -------------------------------------------------------
    dyn_m = _DYNAMIC_RE.search(text)
    dynamic = dyn_m.group(1).lower() if dyn_m else None

    # --- validate completeness -----------------------------------------
    missing = []
    if plate_class is None or subtype is None:
        missing.append("type")
    if diameter_in is None:
        missing.append("size")
    if dynamic is None:
        missing.append("dynamic")

    if missing:
        meta.status = "unparseable"
        meta.skip_reason = (
            "unparseable metadata — missing "
            + ", ".join(missing)
            + "; no guess made"
        )
        meta.plate_class = plate_class
        meta.subtype = subtype
        meta.diameter_in = diameter_in
        meta.stroke = stroke
        meta.dynamic = dynamic
        return meta

    instrument_id = f"{plate_class}_{int(diameter_in)}in_{subtype}"
    meta.instrument_id = instrument_id
    meta.plate_class = plate_class
    meta.subtype = subtype
    meta.diameter_in = diameter_in
    meta.stroke = stroke
    meta.dynamic = dynamic
    if subtype == "windgong":
        meta.notes.append("subtype=wind_gong (distinct plate sub-type)")
    if meta.is_transfer_caution:
        meta.notes.append(
            "Chladni crash-type anchors may not transfer; "
            "band-profile metrics only; excluded from aggregate pass/fail"
        )
    return meta


def resolve_model(meta: SampleMeta) -> ModelMapping:
    """Map parsed metadata → PlateInstrument (catalogue or parametric)."""
    if meta.status != "ok" or meta.diameter_in is None or meta.plate_class is None:
        raise ValueError("resolve_model requires status=ok metadata")

    inches = int(meta.diameter_in)
    subtype = meta.subtype or meta.plate_class
    transfer = meta.is_transfer_caution

    if meta.plate_class == "cymbal" and inches in _CYMBAL_INCH_CATALOGUE:
        base = _CYMBAL_INCH_CATALOGUE[inches]
        instr = PlateInstrument(
            name=meta.instrument_id or base.name,
            diameter=base.diameter,
            thickness=base.thickness,
            material=base.material,
            chladni=base.chladni,
            decay_tau_1k=base.decay_tau_1k,
            decay_alpha=base.decay_alpha,
        )
        return ModelMapping(
            instrument_id=meta.instrument_id or base.name,
            instrument=instr,
            provenance="filename_metadata + exact_catalogue_inch_match",
            plate_class="cymbal",
            subtype=subtype,
            diameter_in=float(inches),
            diameter_m=base.diameter,
            thickness_m=base.thickness,
            catalogue_match=base.name,
            transfer_caution=transfer,
        )

    thickness = _THICKNESS_INTERNAL_DEFAULT[meta.plate_class]
    diameter_m = inches * INCH_M
    material = "bronze_B20"
    name = meta.instrument_id or f"{meta.plate_class}_{inches}in"
    instr = PlateInstrument(
        name=name,
        diameter=diameter_m,
        thickness=thickness,
        material=material,
        chladni=None,  # scaled internal_default via PlateInstrument.low_modes
    )
    return ModelMapping(
        instrument_id=name,
        instrument=instr,
        provenance="filename_metadata + internal_default_thickness",
        plate_class=meta.plate_class,
        subtype=subtype,
        diameter_in=float(inches),
        diameter_m=diameter_m,
        thickness_m=thickness,
        catalogue_match=None,
        transfer_caution=transfer,
    )


def classify_paths(paths: List[Path]) -> Tuple[List[SampleMeta], List[SampleMeta], List[SampleMeta]]:
    """Return (ok, skip_pitched, unparseable) lists."""
    ok: List[SampleMeta] = []
    pitched: List[SampleMeta] = []
    bad: List[SampleMeta] = []
    for p in paths:
        m = parse_sample_path(p)
        if m.status == "ok":
            ok.append(m)
        elif m.status == "skip_pitched":
            pitched.append(m)
        else:
            bad.append(m)
    return ok, pitched, bad
