"""Typed, standalone QM-log -> thermochemistry -> geometric-learning bridge.

Default execution is a SYNTHETIC software self-test, never a DFT calculation.
Dependencies: numpy, scipy, matplotlib, torch. No external geometry framework.
See qm_autobridge/README.md for the supported log dialects and evidence limits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import norm, t as student_t
import torch
from torch import Tensor, nn

HARTREE_KCAL = 627.5095
HARTREE_EV = 27.211386245988
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
ELEMENTS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
            "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
            "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
            "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


class ParseError(ValueError):
    """Missing, incomplete or ambiguous quantum output; never substitute zeros."""


def number(value: str) -> float:
    result = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(result):
        raise ParseError('Nonfinite numeric field')
    return result


@dataclass
class QuantumRecord:
    program: str
    source_sha256: str
    electronic_hartree: float
    atomic_numbers: list[int]
    coordinates_angstrom: list[list[float]]
    charge: int
    multiplicity: int
    optimization_converged: bool
    temperature_kelvin: float | None = None
    zpve_hartree: float | None = None
    enthalpy_hartree: float | None = None
    gibbs_hartree: float | None = None
    enthalpy_correction_hartree: float | None = None
    gibbs_correction_hartree: float | None = None
    frequencies_cm1: list[float] = field(default_factory=list)
    normal_modes: dict[int, list[list[float]]] = field(default_factory=dict)
    partial_charges: dict[str, list[float]] = field(default_factory=dict)
    frontier_ev: dict[str, dict[str, float]] = field(default_factory=dict)
    synthetic: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def imaginary_indices(self) -> list[int]:
        return [i for i, frequency in enumerate(self.frequencies_cm1) if frequency < 0]

    def validate_stationary_point(self, transition_state: bool = False) -> None:
        """Necessary harmonic checks only; a mode count cannot replace an IRC."""
        n = len(self.atomic_numbers)
        rank = np.linalg.matrix_rank(np.asarray(self.coordinates_angstrom) - np.mean(self.coordinates_angstrom, axis=0), tol=1e-5)
        vibrational_count = 3*n - (5 if rank == 1 else 6)
        if not self.optimization_converged or n < 2:
            raise ValueError('An explicitly converged molecular optimization is required')
        expected_counts = (3*n,) if self.program == 'ORCA' else (vibrational_count,)
        if len(self.frequencies_cm1) not in expected_counts:
            raise ValueError('Incomplete frequency spectrum; cannot classify stationary point')
        expected = 1 if transition_state else 0
        if len(self.imaginary_indices) != expected:
            raise ValueError(f'Expected {expected} strictly negative frequencies')
        if transition_state:
            mode = self.normal_modes.get(self.imaginary_indices[0])
            if mode is None or np.asarray(mode).shape != (n, 3) or not np.linalg.norm(mode) > 0:
                raise ValueError('Complete nonzero imaginary-mode displacements required')


class QuantumLogParser:
    """Conservative single-job Gaussian 16 / ORCA 5/6 text parser.

    Concatenated jobs and repeated frequency/thermochemistry analyses are rejected
    rather than mixing geometries, orbital sets and thermochemical corrections.
    Standard low-precision Gaussian displacement tables and ORCA row matrices are
    supported; unsupported optional sections stay absent with explicit warnings.
    """

    def parse_file(self, path: str | Path) -> QuantumRecord:
        return self.parse(Path(path).read_text(encoding='utf-8', errors='strict'))

    @staticmethod
    def _last(text: str, pattern: str, required: bool = False) -> float | None:
        matches = re.findall(pattern, text, re.M | re.I)
        if not matches:
            if required:
                raise ParseError(f'Missing required field: {pattern}')
            return None
        return number(matches[-1])

    def parse(self, text: str) -> QuantumRecord:
        gaussian = bool(re.search(r'Gaussian\s+(?:16|\(R\))|SCF Done:', text))
        orca = bool(re.search(r'O\s+R\s+C\s+A|ORCA|FINAL SINGLE POINT ENERGY', text))
        if gaussian == orca:
            raise ParseError('Unknown or mixed program output')
        program = 'Gaussian' if gaussian else 'ORCA'
        termination = r'Normal termination of Gaussian[^\n]*' if gaussian else r'\*?\s*ORCA TERMINATED NORMALLY\s*\*?'
        ends = list(re.finditer(termination, text))
        if len(ends) != 1:
            raise ParseError('Require exactly one normally terminated job; split multi-job logs explicitly')
        if re.search(r'SCF Done:|FINAL SINGLE POINT ENERGY|Entering Gaussian|Program Version|Error termination|aborting', text[ends[0].end():], re.I):
            raise ParseError('New/failed calculation follows normal termination')
        if re.search(r'Error termination|SCF NOT CONVERGED|ORCA finished by error termination', text, re.I):
            raise ParseError('Failure marker in job')
        if len(re.findall(r'Zero-point correction\s*=|Zero point energy\s*\.\.\.', text, re.I)) > 1:
            raise ParseError('Multiple thermochemistry blocks; select one job/temperature explicitly')
        energy_pattern = (rf'SCF Done:\s+E\([^\n]*?\)\s*=\s*({NUMBER})' if gaussian
                          else rf'FINAL SINGLE POINT ENERGY\s+({NUMBER})')
        energies = list(re.finditer(energy_pattern, text))
        if not energies:
            raise ParseError('Missing converged SCF energy')
        thermo_start = re.search(r'Zero-point correction\s*=|Zero point energy\s*\.\.\.', text, re.I)
        if thermo_start and energies[-1].start() > thermo_start.start():
            raise ParseError('SCF after thermal block; stale thermal correction is unsafe')
        if thermo_start and re.search(r'Standard orientation:|Input orientation:|CARTESIAN COORDINATES \(ANGSTR', text[thermo_start.start():]):
            raise ParseError('Geometry after thermochemistry; separate calculation stages explicitly')
        z, xyz = self._coordinates(text, gaussian)
        if gaussian:
            cm = re.findall(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)', text)
        else:
            charges = re.findall(r'Total Charge\s+[^\n]*?(-?\d+)\s*$', text, re.M)
            mults = re.findall(r'Multiplicity\s+[^\n]*?(\d+)\s*$', text, re.M)
            cm = list(zip(charges, mults))
        if not cm or len(set(cm)) != 1:
            raise ParseError('Missing or changing charge/multiplicity')
        charge, multiplicity = map(int, cm[-1])
        electrons = sum(z)-charge
        if multiplicity < 1 or electrons < multiplicity-1 or (electrons-(multiplicity-1)) % 2:
            raise ParseError('Inconsistent electron-count / multiplicity parity')
        record = QuantumRecord(program, hashlib.sha256(text.encode()).hexdigest(), number(energies[-1][1]), z, xyz,
                               charge, multiplicity, bool(re.search(r'Optimization completed|THE OPTIMIZATION HAS CONVERGED', text)),
                               synthetic='SYNTHETIC_QM_AUTOBRIDGE' in text)
        record.temperature_kelvin = self._last(text, rf'Temperature\s+(?:\.\.\.\s*)?({NUMBER})\s*(?:Kelvin|K)')
        record.zpve_hartree = self._last(text, rf'(?:Zero-point correction\s*=|Zero point energy\s*\.\.\.)\s*({NUMBER})')
        record.enthalpy_hartree = self._last(text, rf'(?:Sum of electronic and thermal Enthalpies\s*=|Total Enthalpy\s*\.\.\.)\s*({NUMBER})')
        record.gibbs_hartree = self._last(text, rf'(?:Sum of electronic and thermal Free Energies\s*=|Final Gibbs free energy\s*\.\.\.)\s*({NUMBER})')
        for kind, label in [('enthalpy', 'Enthalpy'), ('gibbs', 'Gibbs Free Energy')]:
            correction = self._last(text, rf'Thermal correction to {label}\s*=\s*({NUMBER})')
            total = getattr(record, kind+'_hartree')
            if total is not None and correction is not None and abs(total-record.electronic_hartree-correction) > 3e-6:
                raise ParseError('Thermal totals and corrections disagree')
            if total is None and correction is not None:
                total = record.electronic_hartree+correction
                setattr(record, kind+'_hartree', total)
            setattr(record, kind+'_correction_hartree', None if total is None else total-record.electronic_hartree)
        record.frequencies_cm1, record.normal_modes = self._frequencies(text, gaussian, z)
        record.partial_charges = self._charges(text, gaussian, z)
        for scheme, charges in record.partial_charges.items():
            if abs(sum(charges)-record.charge) > max(.02, len(z)*1e-5):
                raise ParseError(f'{scheme} charges do not sum to molecular charge')
        record.frontier_ev = self._orbitals(text, gaussian)
        for attribute in ['temperature_kelvin', 'zpve_hartree', 'enthalpy_hartree', 'gibbs_hartree', 'frequencies_cm1', 'partial_charges', 'frontier_ev']:
            value = getattr(record, attribute)
            if value is None or value == [] or value == {}:
                record.warnings.append('Missing optional field: '+attribute)
        if record.temperature_kelvin is not None and abs(record.temperature_kelvin-298.15) > .01:
            record.warnings.append('Thermochemistry is not at 298.15 K; no temperature conversion performed')
        if record.imaginary_indices:
            record.warnings.append('Negative modes require inspection; exactly one is necessary but does not establish reaction connectivity (IRC needed)')
        record.warnings.append('Method, solvent and standard-state compatibility require caller review')
        return record

    @staticmethod
    def _coordinates(text: str, gaussian: bool) -> tuple[list[int], list[list[float]]]:
        lines = text.splitlines()
        headers = [i for i, line in enumerate(lines) if ('orientation:' in line.lower() and ('Standard' in line or 'Input' in line))] if gaussian else [i for i, line in enumerate(lines) if 'CARTESIAN COORDINATES (ANGSTROEM)' in line or 'CARTESIAN COORDINATES (ANGSTROM)' in line]
        if not headers:
            raise ParseError('Missing supported Cartesian coordinate block in angstrom')
        z: list[int] = []
        xyz: list[list[float]] = []
        for line in lines[headers[-1]+1:]:
            fields = line.split()
            valid = (len(fields) == 6 and fields[0].isdigit() and fields[1].isdigit() and fields[2].lstrip('-').isdigit()) if gaussian else (len(fields) == 4 and fields[0] in ELEMENTS[1:])
            if valid:
                if gaussian and int(fields[0]) != len(z)+1:
                    raise ParseError('Coordinate indices not contiguous')
                atomic_number = int(fields[1]) if gaussian else ELEMENTS.index(fields[0])
                if not 1 <= atomic_number <= 118:
                    raise ParseError('Ghost/dummy/unknown atom unsupported')
                z.append(atomic_number)
                xyz.append([number(x) for x in fields[-3:]])
            elif z:
                if (gaussian and not re.fullmatch(r'\s*-+\s*', line)) or (not gaussian and line.strip()):
                    raise ParseError('Malformed or unterminated Cartesian coordinate table')
                break
        if not z:
            raise ParseError('Empty Cartesian block')
        return z, xyz

    @staticmethod
    def _frequencies(text: str, gaussian: bool, z: list[int]) -> tuple[list[float], dict[int, list[list[float]]]]:
        lines = text.splitlines()
        frequencies: list[float] = []
        modes: dict[int, list[list[float]]] = {}
        if gaussian:
            if len(re.findall(r'Harmonic frequencies', text, re.I)) > 1:
                raise ParseError('Repeated frequency analyses')
            for i, line in enumerate(lines):
                if 'Frequencies --' not in line:
                    continue
                vals = [number(x) for x in line.split('--', 1)[1].split()]
                offset = len(frequencies)
                frequencies.extend(vals)
                for j in range(i+1, min(i+12, len(lines))):
                    if re.match(r'\s*Atom\s+AN\s+', lines[j]):
                        block = [row.split() for row in lines[j+1:j+1+len(z)]]
                        if any(len(row) != 2+3*len(vals) for row in block) or len(block) != len(z):
                            raise ParseError('Incomplete Gaussian displacement table')
                        for atom, row in enumerate(block):
                            if int(row[0]) != atom+1 or int(row[1]) != z[atom]:
                                raise ParseError('Normal mode atom mapping mismatch')
                        for k in range(len(vals)):
                            modes[offset+k] = [[number(x) for x in row[2+3*k:5+3*k]] for row in block]
                        break
        else:
            if text.count('VIBRATIONAL FREQUENCIES') > 1:
                raise ParseError('Repeated frequency analyses')
            match = re.search(r'VIBRATIONAL FREQUENCIES(.*?)(?:NORMAL MODES|THERMOCHEMISTRY|\Z)', text, re.S)
            if match:
                rows = re.findall(rf'^\s*(\d+):\s*({NUMBER})\s*cm', match[1], re.M)
                if [int(row[0]) for row in rows] != list(range(len(rows))):
                    raise ParseError('Noncontiguous ORCA frequency indices')
                frequencies = [number(row[1]) for row in rows]
            normal = re.search(r'NORMAL MODES(.*?)(?:IR SPECTRUM|RAMAN SPECTRUM|THERMOCHEMISTRY|\Z)', text, re.S)
            if normal:
                columns: list[int] = []
                matrix: dict[int, dict[int, float]] = {}
                for line in normal[1].splitlines():
                    fields = line.split()
                    if fields and all(re.fullmatch(r'\d+', f) for f in fields):
                        columns = list(map(int, fields))
                    elif columns and len(fields) == len(columns)+1 and fields[0].isdigit():
                        row = int(fields[0])
                        for column, value in zip(columns, fields[1:]):
                            if column >= len(frequencies) or row >= 3*len(z) or row in matrix.setdefault(column, {}):
                                raise ParseError('Invalid or duplicate normal mode matrix index')
                            matrix[column][row] = number(value)
                for column, rows in matrix.items():
                    if set(rows) != set(range(3*len(z))):
                        raise ParseError('Incomplete ORCA normal mode')
                    modes[column] = np.array([rows[i] for i in range(3*len(z))]).reshape(-1, 3).tolist()
        return frequencies, modes

    @staticmethod
    def _charges(text: str, gaussian: bool, z: list[int]) -> dict[str, list[float]]:
        lines = text.splitlines()
        output: dict[str, list[float]] = {}
        for scheme, pattern in [('mulliken', r'^\s*Mulliken (?:atomic )?charges(?: and spin (?:densities|populations))?\s*:?\s*$'),
                                ('hirshfeld', r'^\s*(?:Hirshfeld charges.*|HIRSHFELD ANALYSIS)\s*$')]:
            headers = [i for i, line in enumerate(lines) if re.match(pattern, line, re.I)]
            if not headers:
                continue
            values: list[float] = []
            for line in lines[headers[-1]+1:]:
                row = re.match(rf'^\s*(\d+)\s+([A-Z][a-z]?)\s*:?\s*({NUMBER})(?:\s|$)', line)
                if row:
                    if int(row[1]) != len(values)+(1 if gaussian else 0) or row[2] != ELEMENTS[z[len(values)]]:
                        raise ParseError('Atomic charge mapping mismatch')
                    values.append(number(row[3]))
                    if len(values) == len(z):
                        break
                elif values:
                    break
            if len(values) != len(z):
                raise ParseError('Incomplete atomic charge table')
            output[scheme] = values
        return output

    @staticmethod
    def _orbitals(text: str, gaussian: bool) -> dict[str, dict[str, float]]:
        occupied: dict[str, list[float]] = {}
        virtual: dict[str, list[float]] = {}
        if gaussian:
            # Orbital output from the last SCF cycle only; do not concatenate old cycles.
            tail = text[text.rfind('SCF Done:'):]
            for spin, kind, values in re.findall(r'(Alpha|Beta)\s+(occ\.|virt\.) eigenvalues\s*--([^\n]+)', tail):
                target = occupied if kind == 'occ.' else virtual
                target.setdefault(spin.lower(), []).extend(number(x)*HARTREE_EV for x in re.findall(NUMBER, values))
        else:
            tail = text[text.rfind('ORBITAL ENERGIES'):] if 'ORBITAL ENERGIES' in text else ''
            spin = 'restricted'
            for line in tail.splitlines()[1:]:
                if 'SPIN UP ORBITALS' in line:
                    spin = 'alpha'
                elif 'SPIN DOWN ORBITALS' in line:
                    spin = 'beta'
                row = re.match(rf'^\s*\d+\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*$', line)
                if row:
                    occupation = number(row[1])
                    if occupation < 0 or occupation > 2:
                        raise ParseError('Invalid orbital occupation')
                    if min(abs(occupation-i) for i in (0, 1, 2)) > 1e-6:
                        raise ParseError('Fractional orbital occupations unsupported for frontier assignment')
                    (occupied if occupation > 0 else virtual).setdefault(spin, []).append(number(row[2])*HARTREE_EV)
                elif occupied and line.strip() and not re.match(r'^\s*(?:[-]+|NO\s|SPIN)', line):
                    break
        result: dict[str, dict[str, float]] = {}
        for spin in occupied:
            if spin in virtual:
                homo, lumo = max(occupied[spin]), min(virtual[spin])
                result[spin] = dict(homo=homo, lumo=lumo, orbital_gap=lumo-homo)
        return result


@dataclass
class BEPFit:
    alpha: float
    beta: float
    covariance: list[list[float]]
    degrees_of_freedom: int
    absolute_sigma: bool
    residual_sum_squares_weighted: float
    x: list[float]
    y: list[float]
    sigma: list[float] | None

    def confidence_band(self, x: Sequence[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        values = np.asarray(x, dtype=float)
        if values.ndim != 1 or not np.isfinite(values).all():
            raise ValueError('Finite one-dimensional prediction coordinates required')
        design = np.column_stack([values, np.ones_like(values)])
        mean = design @ np.array([self.alpha, self.beta])
        sd = np.sqrt(np.maximum(0, np.einsum('ij,jk,ik->i', design, self.covariance, design)))
        critical = norm.ppf(.975) if self.absolute_sigma else student_t.ppf(.975, self.degrees_of_freedom)
        return mean, mean-critical*sd, mean+critical*sd


class ReactionThermochemistryEngine:
    """Closed-composition, equal-charge pathway arithmetic and weighted BEP fit."""

    @staticmethod
    def pathway(states: Mapping[str, QuantumRecord], *, compatible_protocol_confirmed: bool = False) -> dict[str, Any]:
        labels = ['R', 'PC', 'TS', 'INT', 'P']
        if set(states) != set(labels):
            raise ValueError('Supply exactly R, PC, TS, INT, P')
        if not compatible_protocol_confirmed:
            raise ValueError('Caller must review method, solvent, standard state and reference compatibility')
        reference = states['R']
        for name, record in states.items():
            if Counter(record.atomic_numbers) != Counter(reference.atomic_numbers) or (record.charge, record.multiplicity) != (reference.charge, reference.multiplicity):
                raise ValueError('Unbalanced composition, charge or spin; electron/proton reservoirs are not implicit')
            if record.gibbs_hartree is None or not math.isfinite(record.gibbs_hartree) or record.temperature_kelvin is None or abs(record.temperature_kelvin-298.15) > .01:
                raise ValueError('Finite Gibbs energies at 298.15 K required')
            if record.synthetic != reference.synthetic:
                raise ValueError('Do not mix synthetic and external evidence')
            record.validate_stationary_point(name == 'TS')
        g = {name: (states[name].gibbs_hartree-reference.gibbs_hartree)*HARTREE_KCAL for name in labels}
        return dict(relative_g_kcal_mol=g, activation_from_R_kcal_mol=g['TS'],
                    activation_from_PC_kcal_mol=g['TS']-g['PC'], reaction_kcal_mol=g['P'],
                    synthetic=reference.synthetic, temperature_kelvin=298.15,
                    scope='Closed-composition illustrative pathway; no electrochemical reference conversion or IRC validation')

    @staticmethod
    def fit_bep(x: Sequence[float], y: Sequence[float], sigma: Sequence[float] | None = None) -> BEPFit:
        xs, ys = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        if xs.ndim != 1 or ys.shape != xs.shape or len(xs) < 3 or not np.isfinite(xs).all() or not np.isfinite(ys).all():
            raise ValueError('At least three finite paired BEP observations required')
        sig = np.ones_like(xs) if sigma is None else np.asarray(sigma, dtype=float)
        if sig.shape != xs.shape or not np.isfinite(sig).all() or np.any(sig <= 0):
            raise ValueError('Uncertainties must be finite positive standard deviations')
        design = np.column_stack([xs, np.ones_like(xs)])/sig[:, None]
        u, singular, vt = np.linalg.svd(design, full_matrices=False)
        if singular[-1] <= singular[0]*1e-12:
            raise ValueError('Rank-deficient or ill-conditioned BEP design')
        coefficients = vt.T @ ((u.T @ (ys/sig))/singular)
        covariance = (vt.T / singular**2) @ vt
        residual = float(np.sum(((ys-np.column_stack([xs, np.ones_like(xs)]) @ coefficients)/sig)**2))
        if sigma is None:
            covariance *= residual/(len(xs)-2)
        return BEPFit(float(coefficients[0]), float(coefficients[1]), covariance.tolist(), len(xs)-2,
                      sigma is not None, residual, xs.tolist(), ys.tolist(), None if sigma is None else sig.tolist())


class GaussianRBF(nn.Module):
    def __init__(self, bins: int = 32, cutoff: float = 5.0) -> None:
        super().__init__()
        if bins < 2 or not math.isfinite(cutoff) or cutoff <= 0:
            raise ValueError('RBF requires >=2 bins and a positive cutoff')
        self.cutoff = cutoff
        self.register_buffer('centers', torch.linspace(0, cutoff, bins))
        self.gamma = .5 / (cutoff/(bins-1))**2

    def envelope(self, distances: Tensor) -> Tensor:
        return torch.where(distances < self.cutoff, .5*(torch.cos(math.pi*distances/self.cutoff)+1), torch.zeros_like(distances))

    def forward(self, distances: Tensor) -> Tensor:
        return torch.exp(-self.gamma*(distances[..., None]-self.centers)**2)


class ContinuousFilterBlock(nn.Module):
    def __init__(self, width: int, bins: int) -> None:
        super().__init__()
        self.filters = nn.Sequential(nn.Linear(bins, width), nn.SiLU(), nn.Linear(width, width))
        self.project = nn.Linear(width, width, bias=False)
        self.update = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, width))

    def forward(self, h: Tensor, rbf: Tensor, envelope: Tensor, edge_mask: Tensor) -> Tensor:
        # Apply cutoff after the filter MLP, so its biases cannot leak beyond r_cut.
        weights = self.filters(rbf)*envelope[..., None]*edge_mask[..., None]
        messages = torch.sum(weights*self.project(h).unsqueeze(0), dim=1)
        return h+self.update(messages)


class SchNetE3InvariantEngine(nn.Module):
    """Distance-only E(3)-invariant scalars; atom outputs permute with atoms.

    Single molecule per forward call; O(N^2) pairs with a configurable atom cap.
    Two scalar heads are independently trainable AIP/electrode-potential targets;
    neither is a conversion from the other. Atomic barriers are latent additive
    contributions, not uniquely identifiable physical atomic activation energies.
    """
    def __init__(self, width: int = 48, blocks: int = 3, bins: int = 32, cutoff: float = 5.0, max_atoms: int = 512) -> None:
        super().__init__()
        if min(width, blocks, max_atoms) < 1:
            raise ValueError('Positive model dimensions required')
        self.max_atoms = max_atoms
        self.embedding = nn.Embedding(119, width)
        self.state_embedding = nn.Sequential(nn.Linear(2, width), nn.SiLU(), nn.Linear(width, width))
        self.rbf = GaussianRBF(bins, cutoff)
        self.blocks = nn.ModuleList([ContinuousFilterBlock(width, bins) for _ in range(blocks)])
        self.property_head = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, 2))
        self.barrier_head = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, 1))

    def forward(self, atomic_numbers: Tensor, positions: Tensor, charge: int = 0, multiplicity: int = 1) -> dict[str, Tensor]:
        if atomic_numbers.ndim != 1 or not 0 < len(atomic_numbers) <= self.max_atoms or positions.shape != (len(atomic_numbers), 3):
            raise ValueError('Expected Z[N], positions[N,3] within configured atom cap')
        if atomic_numbers.dtype != torch.long or not torch.all((atomic_numbers >= 1) & (atomic_numbers <= 118)):
            raise ValueError('Atomic numbers must be int64 in 1..118')
        if not positions.is_floating_point() or not torch.isfinite(positions).all():
            raise ValueError('Finite floating-point coordinates required')
        if positions.dtype != self.embedding.weight.dtype or positions.device != self.embedding.weight.device or atomic_numbers.device != positions.device:
            raise ValueError('Move model, coordinates and atom numbers to compatible dtype/device')
        if type(charge) is not int or type(multiplicity) is not int or multiplicity < 1:
            raise ValueError('Integer charge and positive integer multiplicity required')
        displacement = positions[:, None, :] - positions[None, :, :]
        # Nonzero diagonal avoids undefined sqrt derivative at self-distance zero.
        diagonal = torch.eye(len(positions), device=positions.device, dtype=positions.dtype)
        squared = displacement.square().sum(-1)
        if torch.any((squared < 1e-16) & ~diagonal.bool()):
            raise ValueError('Coincident atoms are unsupported')
        distances = torch.sqrt(squared+diagonal)
        edge_mask = ~diagonal.bool()
        state = positions.new_tensor([charge, multiplicity-1])
        h = self.embedding(atomic_numbers)+self.state_embedding(state)
        for block in self.blocks:
            h = block(h, self.rbf(distances), self.rbf.envelope(distances), edge_mask)
        properties = self.property_head(h).mean(0)
        atomic_barriers = self.barrier_head(h).squeeze(-1)
        return {'aip': properties[0], 'eox': properties[1], 'atomic_barriers': atomic_barriers,
                'barrier_sum': atomic_barriers.sum()}


# One fixed RDKit ETKDG embedding (seed 42), not optimized or quantum calculated.
# SMILES COc1ccc2[nH]c(-c3ccccc3)cc2c1; formula C15H13NO.
GEOMETRY = '''C -5.303304 .334088 .679589
O -4.776903 -.075180 -.549602
C -3.441805 -.222368 -.850297
C -2.968624 -.630031 -2.081067
C -1.616495 -.763556 -2.336433
C -.678544 -.492747 -1.364278
N .678128 -.531237 -1.318617
C 1.100754 -.158691 -.088392
C 2.479752 -.069567 .363539
C 3.503185 -.383481 -.502836
C 4.825005 -.294021 -.056915
C 5.080768 .111381 1.256541
C 4.046196 .423035 2.114410
C 2.735422 .330117 1.658905
C -.036596 .116492 .641500
C -1.148441 -.080625 -.119014
C -2.511896 .052179 .132675
H -4.824113 1.196486 1.144677
H -6.384913 .572581 .517368
H -5.305236 -.488116 1.440666
H -3.689212 -.841523 -2.838733
H -1.317772 -1.086001 -3.319336
H 1.247660 -.815742 -2.148028
H 3.348596 -.700594 -1.524610
H 5.660892 -.538925 -.728202
H 6.100087 .184007 1.611976
H 4.198432 .739751 3.134673
H 1.908441 .576180 2.334974
H -.054491 .443049 1.671949
H -2.854973 .378717 1.122922'''


class SyntheticDFTSandbox:
    """Format-realistic numerical fixtures; not physical DFT/Hessian results."""
    relative_g = {'R': 0., 'PC': -3., 'TS': 18., 'INT': -5., 'P': -9.}

    @staticmethod
    def generate(program: str, state: str) -> str:
        if program not in ('Gaussian', 'ORCA') or state not in ('substrate', 'cation', 'R', 'PC', 'TS', 'INT', 'P'):
            raise ValueError('Unknown fixture program or state')
        atoms = [line.split() for line in GEOMETRY.splitlines()]
        n = len(atoms)
        charge, mult = (0, 1) if state == 'substrate' else (1, 2)
        g = -650.+SyntheticDFTSandbox.relative_g.get(state, 0.)/HARTREE_KCAL
        if state == 'substrate':
            g -= .2
        energy, h, zpve = g-.12, g+.03, .10
        vib = [(-342.5 if state == 'TS' else 80.)]+[100.+i*35 for i in range(3*n-7)]
        # Arbitrary displacement fixture emphasizes the C3-H pair (atoms 15,29).
        displacement = np.zeros((n, 3)); displacement[14, 0] = -.1; displacement[28, 0] = 1.
        q = [charge/n]*n
        electrons = sum(ELEMENTS.index(a[0]) for a in atoms)-charge
        n_alpha, n_beta = (electrons+mult-1)//2, (electrons-mult+1)//2
        alpha = np.linspace(-20., -.21, n_alpha).tolist()
        beta = np.linspace(-20.01, -.24, n_beta).tolist()
        rows = ['SYNTHETIC_QM_AUTOBRIDGE: invented energies, frequencies, charges and modes; SOFTWARE TEST ONLY']
        if program == 'Gaussian':
            rows += ['Gaussian 16', f'Charge = {charge} Multiplicity = {mult}',
                     f'SCF Done: E(UB3LYP) = {energy:.12f} A.U. after 10 cycles', 'Standard orientation:',
                     '-----', 'Center Atomic Atomic Coordinates (Angstroms)', 'Number Number Type X Y Z', '-----']
            rows += [f'{i+1} {ELEMENTS.index(a[0])} 0 '+ ' '.join(a[1:]) for i, a in enumerate(atoms)]
            rows += ['-----', 'Optimization completed.']
            for spin, orbitals, virtuals in [('Alpha', alpha, [.03,.1])]+([('Beta', beta, [-.02,.11])] if mult == 2 else []):
                for start in range(0, len(orbitals), 5):
                    rows.append(f'{spin}  occ. eigenvalues -- '+' '.join(f'{v:.6f}' for v in orbitals[start:start+5]))
                rows.append(f'{spin} virt. eigenvalues -- '+' '.join(f'{v:.6f}' for v in virtuals))
            for name in ['Mulliken charges:', 'Hirshfeld charges, spin densities, dipoles, and CM5 charges']:
                rows.append(name)
                rows += [f'{i+1} {a[0]} {q[i]:.9f} 0.0000' for i, a in enumerate(atoms)]
                rows.append('End charge table')
            rows.append('Harmonic frequencies (cm**-1)')
            for start in range(0, len(vib), 3):
                vals = vib[start:start+3]
                rows += ['Frequencies -- '+' '.join(f'{v:.6f}' for v in vals), 'Red. masses -- 1.0 1.0 1.0', 'Atom AN '+ ' X Y Z'*len(vals)]
                rows += [f'{i+1} {ELEMENTS.index(a[0])} '+' '.join(f'{v:.6f}' for _ in vals for v in displacement[i]) for i, a in enumerate(atoms)]
            rows += ['Temperature 298.150 Kelvin. Pressure 1.00000 Atm.', f'Zero-point correction= {zpve:.12f}',
                     f'Thermal correction to Enthalpy= {h-energy:.12f}', f'Thermal correction to Gibbs Free Energy= {g-energy:.12f}',
                     f'Sum of electronic and thermal Enthalpies= {h:.12f}', f'Sum of electronic and thermal Free Energies= {g:.12f}',
                     'Normal termination of Gaussian 16']
        else:
            rows += ['O   R   C   A', 'Program Version 6.0.0', f'Total Charge .... {charge}', f'Multiplicity .... {mult}',
                     f'FINAL SINGLE POINT ENERGY {energy:.12f}', 'CARTESIAN COORDINATES (ANGSTROEM)', '-----', GEOMETRY, '',
                     'THE OPTIMIZATION HAS CONVERGED', 'ORBITAL ENERGIES']
            channels = [('SPIN UP ORBITALS', alpha, [.03,.1]), ('SPIN DOWN ORBITALS', beta, [-.02,.11])] if mult == 2 else [('',alpha,[.03,.1])]
            for label, orbitals, virtuals in channels:
                rows += [label, 'NO OCC E(Eh) E(eV)']
                for i, eigenvalue in enumerate(orbitals+virtuals):
                    occupation = (1 if mult == 2 else 2) if i < len(orbitals) else 0
                    rows.append(f'{i} {occupation:.4f} {eigenvalue:.6f} {eigenvalue*HARTREE_EV:.6f}')
            rows.append('')
            for name in ['MULLIKEN ATOMIC CHARGES', 'HIRSHFELD ANALYSIS']:
                rows.append(name)
                rows += [f'{i} {a[0]} : {q[i]:.9f}' for i, a in enumerate(atoms)]
                rows += ['', 'End charge table']
            freqs = [0.]*6+vib
            rows += ['VIBRATIONAL FREQUENCIES', '-----']+[f'{i}: {v:.6f} cm**-1' for i, v in enumerate(freqs)]
            rows += ['NORMAL MODES', 'These modes are synthetic, not Hessian eigenvectors.']
            for start in range(0, len(freqs), 6):
                columns = list(range(start, min(start+6, len(freqs))))
                rows.append(' '.join(map(str, columns)))
                rows += [f'{i} '+' '.join(f'{displacement.ravel()[i]:.6f}' for _ in columns) for i in range(3*n)]
            rows += ['IR SPECTRUM', 'THERMOCHEMISTRY AT 298.15K', 'Temperature ... 298.150 K',
                     f'Zero point energy ... {zpve:.12f} Eh', f'Total Enthalpy ... {h:.12f} Eh',
                     f'Final Gibbs free energy ... {g:.12f} Eh', '****ORCA TERMINATED NORMALLY****']
        return '\n'.join(rows)+'\n'

    def write_and_verify(self, directory: Path) -> dict[str, dict[str, QuantumRecord]]:
        directory.mkdir(parents=True, exist_ok=False)
        records: dict[str, dict[str, QuantumRecord]] = {}
        parser = QuantumLogParser()
        for program in ('Gaussian', 'ORCA'):
            records[program] = {}
            for state in ('substrate', 'cation', 'R', 'PC', 'TS', 'INT', 'P'):
                text = self.generate(program, state)
                path = directory/f'{program.lower()}_{state}.out'
                path.write_text(text, encoding='utf-8', newline='\n')
                record = parser.parse_file(path)
                record.validate_stationary_point(state == 'TS')
                expected = -650.+self.relative_g.get(state, 0.)/HARTREE_KCAL-(.2 if state == 'substrate' else 0.)
                if abs(record.gibbs_hartree-expected) > 1e-10 or abs(record.zpve_hartree-.1) > 1e-10:
                    raise AssertionError('Synthetic thermochemistry target mismatch')
                if state == 'TS' and record.frequencies_cm1[record.imaginary_indices[0]] != -342.5:
                    raise AssertionError('Synthetic imaginary frequency mismatch')
                if set(record.partial_charges) != {'mulliken', 'hirshfeld'} or not record.frontier_ev:
                    raise AssertionError('Missing synthetic charges/orbitals')
                records[program][state] = record
        for state in records['Gaussian']:
            a, b = records['Gaussian'][state], records['ORCA'][state]
            if a.electronic_hartree != b.electronic_hartree or a.gibbs_hartree != b.gibbs_hartree or a.coordinates_angstrom != b.coordinates_angstrom:
                raise AssertionError('Cross-format target mismatch')
        return records


class PublicationArtifactGenerator:
    """ACS-style booktabs table and 300-dpi figure with mandatory provenance labels."""

    @staticmethod
    def generate(directory: Path, states: Mapping[str, QuantumRecord], pathway: dict[str, Any], fit: BEPFit,
                 model: SchNetE3InvariantEngine, database: dict[str, Any]) -> None:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        lines = [r'\begin{table}[htbp]', r'\centering',
                 r'\caption{SYNTHETIC software fixtures, not DFT results. Energies in hartree; activation energies relative to R in kcal mol$^{-1}$.}',
                 r'\label{tab:qm-autobridge}', r'\begin{tabular}{lrrrrrr}', r'\toprule',
                 r'State & $E_{\mathrm{SCF}}$ & ZPVE & $H$ & $G$ & $\Delta G^\ddagger$ & $\nu_{\mathrm{imag}}$/cm$^{-1}$ \\', r'\midrule']
        for name, record in states.items():
            barrier = f"{pathway['activation_from_R_kcal_mol']:.2f}" if name == 'TS' else '--'
            imag = f'{record.frequencies_cm1[record.imaginary_indices[0]]:.1f}' if record.imaginary_indices else '--'
            lines.append(f'{name} & {record.electronic_hartree:.8f} & {record.zpve_hartree:.8f} & {record.enthalpy_hartree:.8f} & {record.gibbs_hartree:.8f} & {barrier} & {imag} '+r'\\')
        lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}',
                  '% Requires booktabs. H/G are TOTAL energies; do not add ZPVE again.',
                  '% All states in this toy pathway have identical composition, +1 charge and doublet spin.']
        (directory/'qm_autobridge_table.tex').write_text('\n'.join(lines)+'\n', encoding='utf-8')
        standalone = [r'\documentclass[10pt,landscape]{article}', r'\usepackage[margin=1.5cm]{geometry}',
                      r'\usepackage{booktabs}', r'\begin{document}', *lines,
                      r'\noindent Software demonstration only. One imaginary mode does not prove reaction connectivity. The neural network is untrained.',
                      r'\end{document}']
        (directory/'qm_autobridge_report.tex').write_text('\n'.join(standalone)+'\n', encoding='utf-8')
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True)
        labels, values = list(pathway['relative_g_kcal_mol']), list(pathway['relative_g_kcal_mol'].values())
        for i, value in enumerate(values):
            axes[0].plot([i-.22, i+.22], [value, value], color='#175b89', lw=2.5)
            if i:
                axes[0].plot([i-1+.22, i-.22], [values[i-1], value], '--', color='gray')
        axes[0].set(xticks=range(5), xticklabels=labels, ylabel='Relative Gibbs energy (kcal/mol)', title='A  Synthetic pathway (not DFT)')
        axes[0].annotate('', xy=(2, values[2]), xytext=(2, 0), arrowprops={'arrowstyle':'<->', 'color':'#b64322'})
        axes[0].text(2.14, 8, f"{pathway['activation_from_R_kcal_mol']:.1f}")
        grid = np.linspace(min(fit.x), max(fit.x), 150)
        mean, low, high = fit.confidence_band(grid)
        axes[1].errorbar(fit.x, fit.y, yerr=fit.sigma, fmt='o', color='#175b89')
        axes[1].plot(grid, mean, color='#b64322', label=f'alpha={fit.alpha:.3f}, beta={fit.beta:.3f}')
        axes[1].fill_between(grid, low, high, alpha=.2, color='#b64322', label='95% mean confidence band')
        axes[1].set(xlabel='Reaction free energy (kcal/mol)', ylabel='Activation free energy (kcal/mol)', title='B  Synthetic BEP calibration')
        axes[1].legend(fontsize=8)
        distance = torch.linspace(0, 6, 240, dtype=next(model.parameters()).dtype)
        with torch.no_grad():
            envelope = model.rbf.envelope(distance)
            radial = model.rbf(distance)*envelope[:, None]
            filters = model.blocks[0].filters(model.rbf(distance))*envelope[:, None]
        for k in (0, 8, 16, 24, 31):
            axes[2].plot(distance.numpy(), radial[:, k].numpy(), alpha=.6, lw=1)
        axes[2].plot(distance.numpy(), envelope.numpy(), 'k--', label='Cosine envelope')
        axes[2].plot(distance.numpy(), filters[:, 0].numpy(), color='#b64322', label='Untrained filter channel 0')
        axes[2].set(xlabel='Distance (angstrom)', ylabel='Basis / filter value (arbitrary)', title='C  RBF and cutoff filters')
        axes[2].legend(fontsize=8)
        fig.suptitle('QM AutoBridge: synthetic integration test; neural network untrained', fontsize=12)
        fig.savefig(directory/'qm_autobridge_composite.png', dpi=300)
        fig.savefig(directory/'qm_autobridge_composite.pdf')
        plt.close(fig)
        database['figure_data'] = dict(distance_angstrom=distance.tolist(), cutoff=envelope.tolist(),
                                       rbf=radial.tolist(), untrained_filter_channel_0=filters[:, 0].tolist(),
                                       bep_x=grid.tolist(), bep_mean=mean.tolist(), bep_lower=low.tolist(), bep_upper=high.tolist())
        (directory/'qm_autobridge_database.json').write_text(json.dumps(database, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def run_pipeline(output: Path, seed: int = 42) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False)
    torch.manual_seed(seed)
    torch.set_num_threads(2)
    records = SyntheticDFTSandbox().write_and_verify(output/'synthetic_logs')
    states = {name: records['Gaussian'][name] for name in SyntheticDFTSandbox.relative_g}
    pathway = ReactionThermochemistryEngine.pathway(states, compatible_protocol_confirmed=True)
    x = [-20., -15., -10., -5., 0., 5.]
    y = [.45*v+18.+noise for v, noise in zip(x, [.3, -.4, .2, -.1, .35, -.2])]
    fit = ReactionThermochemistryEngine.fit_bep(x, y, [.6]*len(x))
    model = SchNetE3InvariantEngine().double()
    z = torch.tensor(states['R'].atomic_numbers, dtype=torch.long)
    xyz = torch.tensor(states['R'].coordinates_angstrom, dtype=torch.float64, requires_grad=True)
    prediction = model(z, xyz, charge=1, multiplicity=2)
    loss = sum(value.square().mean() for value in prediction.values())
    loss.backward()
    if not torch.isfinite(xyz.grad).all() or not all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters()):
        raise AssertionError('Nonfinite or missing network gradients')
    with torch.no_grad():
        rotation, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64))
        moved = model(z, xyz @ rotation+torch.tensor([2., -3., 4.]), charge=1, multiplicity=2)
        permutation = torch.arange(len(z)-1, -1, -1)
        permuted = model(z[permutation], xyz[permutation], charge=1, multiplicity=2)
        invariance = max(abs(float(prediction[key]-moved[key])) for key in ('aip', 'eox', 'barrier_sum'))
        permutation_error = float(torch.max(torch.abs(prediction['atomic_barriers'][permutation]-permuted['atomic_barriers'])))
    if max(invariance, permutation_error) > 1e-9:
        raise AssertionError('Geometric symmetry test failed')
    diagnostics = dict(parsed_logs=14, synthetic=True, actual_dft_jobs=0, model_trained=False,
                       invariance_max_abs=invariance, permutation_max_abs=permutation_error,
                       finite_gradients=True, seed=seed)
    database = dict(schema_version=1, provenance={'role':'synthetic_software_self_test', 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    'torch':torch.__version__, 'numpy':np.__version__, 'structure':'5-methoxy-2-phenyl-1H-indole',
                    'geometry':'Fixed RDKit ETKDG embedding, seed 42, not optimized',
                    'limitations':'No experimental data, physical Hessian, IRC, trained predictor, AIP or electrode potential validation.'},
                    records={program:{state:asdict(record) for state, record in group.items()} for program, group in records.items()},
                    pathway=pathway, bep=asdict(fit), bep_synthetic_center_labels=['Cu-N4','Cu-N3C1','Cu-N2O2','Cu-N2S2','Cu-N3O1','Cu-N3S1'], diagnostics=diagnostics,
                    network_untrained_arbitrary_units={key:value.detach().tolist() for key, value in prediction.items()})
    PublicationArtifactGenerator.generate(output, states, pathway, fit, model, database)
    return diagnostics


def validate_saved(directory: Path) -> dict[str, Any]:
    """Reparse archived logs and independently recompute the saved arithmetic."""
    database = json.loads((directory/'qm_autobridge_database.json').read_text('utf-8'))
    checks: dict[str, bool] = {}
    checks['source_hash'] = database['provenance']['source_sha256'] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    parsed: dict[str, dict[str, QuantumRecord]] = {}
    for program, records in database['records'].items():
        parsed[program] = {}
        for name, saved in records.items():
            record = QuantumLogParser().parse_file(directory/'synthetic_logs'/f'{program.lower()}_{name}.out')
            record.validate_stationary_point(name == 'TS')
            # JSON serializes integer dictionary keys as strings.
            checks[f'{program}_{name}'] = json.loads(json.dumps(asdict(record))) == saved
            parsed[program][name] = record
    states = {name:parsed['Gaussian'][name] for name in SyntheticDFTSandbox.relative_g}
    pathway = ReactionThermochemistryEngine.pathway(states, compatible_protocol_confirmed=True)
    checks['pathway_recomputed'] = pathway == database['pathway']
    fit = ReactionThermochemistryEngine.fit_bep(database['bep']['x'],database['bep']['y'],database['bep']['sigma'])
    checks['bep_recomputed'] = asdict(fit) == database['bep']
    checks['evidence_labels'] = database['provenance']['role'] == 'synthetic_software_self_test' and database['diagnostics']['actual_dft_jobs'] == 0 and database['diagnostics']['model_trained'] is False
    torch.manual_seed(database['diagnostics']['seed'])
    model = SchNetE3InvariantEngine().double()
    with torch.no_grad():
        result = model(torch.tensor(states['R'].atomic_numbers),torch.tensor(states['R'].coordinates_angstrom,dtype=torch.float64),charge=1,multiplicity=2)
    checks['untrained_inference_reproduced'] = all(np.allclose(value.tolist(),database['network_untrained_arbitrary_units'][key],rtol=1e-10,atol=1e-10) for key,value in result.items())
    mean, low, high = fit.confidence_band(database['figure_data']['bep_x'])
    checks['figure_confidence_band'] = all(np.allclose(actual,database['figure_data'][key],rtol=1e-12,atol=1e-12) for actual,key in [(mean,'bep_mean'),(low,'bep_lower'),(high,'bep_upper')])
    png = (directory/'qm_autobridge_composite.png').read_bytes()
    checks['png_dimensions'] = png[:8] == b'\x89PNG\r\n\x1a\n' and int.from_bytes(png[16:20],'big') == 4200 and int.from_bytes(png[20:24],'big') == 1380
    checks['png_300dpi'] = b'pHYs'+(11811).to_bytes(4,'big')*2+b'\x01' in png
    checks['pdf_export'] = (directory/'qm_autobridge_composite.pdf').read_bytes().startswith(b'%PDF')
    checks['table_synthetic_label'] = 'SYNTHETIC' in (directory/'qm_autobridge_table.tex').read_text('utf-8')
    result_summary = {'passed':all(checks.values()),'checks':checks,'scope':'Software and saved-artifact integrity, not physical validation'}
    if not result_summary['passed']:
        raise ValueError('Saved release validation failed: '+str([key for key,value in checks.items() if not value]))
    return result_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('qm_autobridge_output'), help='New directory; existing directories are never overwritten')
    parser.add_argument('--parse', type=Path, help='Parse one external log instead of running the synthetic sandbox')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--validate-saved', type=Path, help='Validate an existing sandbox release without overwriting it')
    args = parser.parse_args()
    try:
        if args.parse and args.validate_saved:
            raise ValueError('--parse and --validate-saved are mutually exclusive')
        if args.validate_saved:
            print(json.dumps(validate_saved(args.validate_saved),indent=2))
        elif args.parse:
            record = QuantumLogParser().parse_file(args.parse)
            args.output.mkdir(parents=True, exist_ok=False)
            (args.output/'parsed_record.json').write_text(json.dumps(asdict(record), indent=2, allow_nan=False)+'\n', encoding='utf-8')
            print(json.dumps({'parsed':True, 'synthetic':record.synthetic, 'warnings':record.warnings}))
        else:
            print(json.dumps(run_pipeline(args.output, args.seed), indent=2))
    except (ValueError, OSError, RuntimeError) as error:
        parser.exit(2, f'QM AutoBridge failed: {error}\n')


if __name__ == '__main__':
    main()
